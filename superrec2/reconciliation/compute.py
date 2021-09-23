"""Compute reconciliations with a arbitrary event costs."""
from collections import defaultdict
from itertools import product
from typing import DefaultDict, List, NamedTuple, Tuple
from ete3 import PhyloTree, PhyloNode
from ..utils.trees import LowestCommonAncestor
from ..utils.min_sequence import MinSequence
from ..model.tree_mapping import TreeMapping
from .tools import CostType, CostVector, ExtendedIntegral


class MappingInfo(NamedTuple):
    """Information about an assignment of a species to a gene."""

    # Assigned species
    species: PhyloNode

    # Species assigned to the left child
    left: PhyloNode

    # Species assigned to the right child
    right: PhyloNode


def reconcile_lca(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
) -> TreeMapping:
    """
    Compute the minimum-cost reconciliation between species tree and a gene
    tree, allowing for duplications and losses and where a duplication
    costs the same as loss.

    :param gene_tree: gene tree to reconcile
    :param species_lca: species ancestry information
    :returns: the reconciliation result
    """
    rec = {}

    for gene in gene_tree.traverse("postorder"):
        if gene.is_leaf():
            species = species_lca.tree.get_leaves_by_name(gene.species)[0]
            rec[gene] = species
        else:
            left, right = gene.children
            rec[gene] = species_lca(rec[left], rec[right])

    return rec


THLTable = DefaultDict[Tuple[PhyloNode, PhyloNode], MinSequence[MappingInfo]]


def _compute_thl_try_speciation(  # pylint:disable=too-many-locals
    species_lca: LowestCommonAncestor,
    root_species: PhyloNode,
    root_gene: PhyloNode,
    table: THLTable,
    costs: CostVector,
) -> None:
    loss_cost = costs[CostType.FULL_LOSS]

    options = table[(root_gene, root_species)]
    left_species, right_species = root_species.children
    left_gene, right_gene = root_gene.children

    # Optimal costs obtained by mapping the left or right gene below
    # the left or right species
    min_left_to_left: MinSequence[PhyloNode] = MinSequence()
    min_right_to_left: MinSequence[PhyloNode] = MinSequence()
    min_left_to_right: MinSequence[PhyloNode] = MinSequence()
    min_right_to_right: MinSequence[PhyloNode] = MinSequence()

    for left_child in left_species.traverse():
        min_left_to_left.update(
            (table[(left_gene, left_child)].min, left_child)
        )
        min_right_to_left.update(
            (table[(right_gene, left_child)].min, left_child)
        )

    for right_child in right_species.traverse():
        min_left_to_right.update(
            (table[(left_gene, right_child)].min, right_child)
        )
        min_right_to_right.update(
            (table[(right_gene, right_child)].min, right_child)
        )

    # Map left gene below left species and right gene below right
    for left_child, right_child in product(
        min_left_to_left, min_right_to_right
    ):
        options.update(
            (
                min_left_to_left.min
                + min_right_to_right.min
                + loss_cost
                * (
                    species_lca.distance(root_species, left_child)
                    + species_lca.distance(root_species, right_child)
                    - 2
                ),
                MappingInfo(
                    species=root_species,
                    left=left_child,
                    right=right_child,
                ),
            )
        )

    # Map right gene below left species and left gene below right
    for left_child, right_child in product(
        min_right_to_left, min_left_to_right
    ):
        options.update(
            (
                min_right_to_left.min
                + min_left_to_right.min
                + loss_cost
                * (
                    species_lca.distance(root_species, left_child)
                    + species_lca.distance(root_species, right_child)
                    - 2
                ),
                MappingInfo(
                    species=root_species,
                    left=right_child,
                    right=left_child,
                ),
            )
        )


def _compute_thl_try_duplication_transfer(  # pylint:disable=too-many-locals
    species_lca: LowestCommonAncestor,
    root_species: PhyloNode,
    root_gene: PhyloNode,
    table: THLTable,
    costs: CostVector,
) -> None:
    dup_cost = costs[CostType.DUPLICATION]
    hgt_cost = costs[CostType.HORIZONTAL_GENE_TRANSFER]
    loss_cost = costs[CostType.FULL_LOSS]

    options = table[(root_gene, root_species)]
    left_gene, right_gene = root_gene.children

    # Optimal costs obtained by mapping the left or right gene inside
    # root_species’ subtree or outside of it
    min_left_to_child: MinSequence[PhyloNode] = MinSequence()
    min_left_to_sep: MinSequence[PhyloNode] = MinSequence()
    min_right_to_child: MinSequence[PhyloNode] = MinSequence()
    min_right_to_sep: MinSequence[PhyloNode] = MinSequence()

    for other_species in species_lca.tree.traverse():
        if species_lca.is_ancestor_of(root_species, other_species):
            min_left_to_child.update(
                (table[(left_gene, other_species)].min, other_species)
            )
            min_right_to_child.update(
                (table[(right_gene, other_species)].min, other_species)
            )
        elif not species_lca.is_ancestor_of(other_species, root_species):
            min_left_to_sep.update(
                (table[(left_gene, other_species)].min, other_species)
            )
            min_right_to_sep.update(
                (table[(right_gene, other_species)].min, other_species)
            )

    # Try mapping as a duplication
    for left_species, right_species in product(
        min_left_to_child, min_right_to_child
    ):
        options.update(
            (
                dup_cost
                + min_left_to_child.min
                + min_right_to_child.min
                + loss_cost
                * (
                    species_lca.distance(root_species, left_species)
                    + species_lca.distance(root_species, right_species)
                ),
                MappingInfo(
                    species=root_species,
                    left=left_species,
                    right=right_species,
                ),
            )
        )

    # Try mapping as a horizontal gene transfer
    # Map left gene as a transfer and right gene as conserved
    for left_species, right_species in product(
        min_left_to_sep, min_right_to_child
    ):
        options.update(
            (
                hgt_cost
                + min_left_to_sep.min
                + min_right_to_child.min
                + loss_cost * species_lca.distance(root_species, right_species),
                MappingInfo(
                    species=root_species,
                    left=left_species,
                    right=right_species,
                ),
            )
        )

    # Map left gene as conserved and right gene as a transfer
    for left_species, right_species in product(
        min_left_to_child, min_right_to_sep
    ):
        options.update(
            (
                hgt_cost
                + min_left_to_child.min
                + min_right_to_sep.min
                + loss_cost * species_lca.distance(root_species, left_species),
                MappingInfo(
                    species=root_species,
                    left=left_species,
                    right=right_species,
                ),
            )
        )


def _compute_thl_table(
    gene_tree: PhyloTree, species_lca: LowestCommonAncestor, costs: CostVector
) -> THLTable:
    table: THLTable = defaultdict(MinSequence)

    for root_gene in gene_tree.traverse("postorder"):
        if root_gene.is_leaf():
            root_species = species_lca.tree.get_leaves_by_name(
                root_gene.species
            )[0]
            table[(root_gene, root_species)].update(
                (0, MappingInfo(species=root_species, left=None, right=None))
            )
        else:
            for root_species in species_lca.tree.traverse("postorder"):
                if not root_species.is_leaf():
                    _compute_thl_try_speciation(
                        species_lca, root_species, root_gene, table, costs
                    )

                _compute_thl_try_duplication_transfer(
                    species_lca, root_species, root_gene, table, costs
                )

    return table


def _decode_thl_table(root, solutions, table):
    results = []

    for solution in solutions:
        if solution.left is None:
            assert solution.right is None
            results.append({root: solution.species})
        else:
            left, right = root.children
            mappings = product(
                _decode_thl_table(
                    left,
                    table[(left, solution.left)],
                    table,
                ),
                _decode_thl_table(
                    right,
                    table[(right, solution.right)],
                    table,
                ),
            )

            for map_left, map_right in mappings:
                results.append(
                    {
                        root: solution.species,
                        **map_left,
                        **map_right,
                    }
                )

    return results


def reconcile_thl(
    gene_tree: PhyloTree, species_lca: LowestCommonAncestor, costs: CostVector
) -> Tuple[ExtendedIntegral, List[TreeMapping]]:
    """
    Find all minimum-cost reconciliations between a species tree and a gene
    tree, allowing for duplications, horizontal gene transfers and losses.

    If you wish to disallow horizontal gene transfers, using
    `reconcile_lca` is faster than setting `costs[HorizontalGeneTransfer]`
    to infinity.

    :param gene_tree: gene tree to reconcile
    :param species_lca: species ancestry information
    :param costs: cost of each evolutionary event
    :returns: a tuple containing the minimum cost of a reconciliation and a
        list of all reconciliations with such a cost
    """
    table = _compute_thl_table(gene_tree, species_lca, costs or {})
    solutions: MinSequence[MappingInfo] = MinSequence()

    for root_species in species_lca.tree.traverse("postorder"):
        for option in table[(gene_tree, root_species)]:
            solutions.update((table[(gene_tree, root_species)].min, option))

    return solutions.min, _decode_thl_table(gene_tree, solutions, table)
