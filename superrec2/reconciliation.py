from collections import defaultdict, namedtuple
from enum import Enum, auto
from itertools import product
from typing import List, Mapping, Tuple, Union
from numbers import Integral
from ete3 import PhyloTree, PhyloNode
from infinity import Infinity
from .lowest_common_ancestor import LowestCommonAncestor
from .min_sequence import MinSequence


ExtendedIntegral = Union[Integral, Infinity]
Reconciliation = Mapping[PhyloNode, PhyloNode]


class Event(Enum):
    """Evolutionary events affecting genes."""
    # Sentinel for extant (current) genes
    Leaf = auto()

    # Sentinel for scenarios that are invalid wrt our evolutionary model
    Invalid = auto()

    # Transmission of the parent gene to both children species
    Speciation = auto()

    # Duplication of the parent gene in the same genome
    Duplication = auto()

    # Transfer of the parent gene to a foreign genome
    HorizontalGeneTransfer = auto()


def get_event(
    gene_node: PhyloNode,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
) -> Event:
    """
    Find the event associated to a gene by a reconciliation.

    :param gene_node: gene to query
    :param species_lca: species ancestry information
    :param rec: reconciliation to use
    :returns: event associated to the given gene node
    """
    if gene_node.is_leaf():
        return (
            Event.Leaf if rec[gene_node].name == gene_node.species
            else Event.Invalid
        )

    vl, vr = gene_node.children

    if (
        species_lca.is_strict_ancestor_of(rec[vl], rec[gene_node])
        or species_lca.is_strict_ancestor_of(rec[vr], rec[gene_node])
    ):
        return Event.Invalid

    if (
        species_lca.is_ancestor_of(rec[gene_node], rec[vl])
        and species_lca.is_ancestor_of(rec[gene_node], rec[vr])
    ):
        return Event.Speciation if (
            rec[gene_node] == species_lca(rec[vl], rec[vr])
            and not species_lca.is_comparable(rec[vl], rec[vr])
        ) else Event.Duplication

    if (
        species_lca.is_ancestor_of(rec[gene_node], rec[vl])
        or species_lca.is_ancestor_of(rec[gene_node], rec[vr])
    ):
        return Event.HorizontalGeneTransfer

    return Event.Invalid


class CostType(Enum):
    """Evolutionary events to which a cost can be assigned."""
    Duplication = auto()
    HorizontalGeneTransfer = auto()
    Loss = auto()


# Cost values for each possible evolutionary event
CostVector = Mapping[CostType, ExtendedIntegral]


def get_cost(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
    costs: CostVector
) -> ExtendedIntegral:
    """
    Compute the total cost of a reconciliation.

    :param gene_tree: reconciled gene tree
    :param species_lca: ancestry information about the reconciled species tree
    :param rec: reconciliation to use
    :param costs: costs for each evolutionary event
    """
    event = get_event(gene_tree, species_lca, rec)

    if event == Event.Invalid:
        return float('inf')

    if event == Event.Leaf:
        return 0

    vl, vr = gene_tree.children
    cost_vl = get_cost(vl, species_lca, rec, costs)
    dist_vl = species_lca.distance(rec[gene_tree], rec[vl])
    cost_vr = get_cost(vr, species_lca, rec, costs)
    dist_vr = species_lca.distance(rec[gene_tree], rec[vr])

    if event == Event.Speciation:
        return (
            cost_vl + cost_vr
            + costs[CostType.Loss] * (dist_vl + dist_vr - 2)
        )

    if event == Event.Duplication:
        return (
            costs[CostType.Duplication]
            + cost_vl + cost_vr
            + costs[CostType.Loss] * (dist_vl + dist_vr)
        )

    assert event == Event.HorizontalGeneTransfer

    dist_conserved = (
        dist_vl if species_lca.is_ancestor_of(rec[gene_tree], rec[vl])
        else dist_vr
    )
    return (
        costs[CostType.HorizontalGeneTransfer]
        + cost_vl + cost_vr
        + costs[CostType.Loss] * dist_conserved
    )


MappingInfo = namedtuple('MappingInfo', ['species', 'left', 'right'])


def _compute_thl_table(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
    costs: CostVector
) -> MinSequence[MappingInfo]:
    dup_cost = costs[CostType.Duplication]
    hgt_cost = costs[CostType.HorizontalGeneTransfer]
    loss_cost = costs[CostType.Loss]

    min_costs = defaultdict(MinSequence)
    species_lca = LowestCommonAncestor(species_tree)

    for v in gene_tree.traverse("postorder"):
        if v.is_leaf():
            s = species_tree.get_leaves_by_name(v.species)[0]
            min_costs[(v, s)].update((0, MappingInfo(
                species=s,
                left=None,
                right=None
            )))
        else:
            vl, vr = v.children

            for u in species_tree.traverse("postorder"):
                options = min_costs[(v, u)]

                # Try mapping as a speciation
                if not u.is_leaf():
                    ul, ur = u.children
                    min_vl_ul_ch = MinSequence()
                    min_vr_ul_ch = MinSequence()
                    min_vl_ur_ch = MinSequence()
                    min_vr_ur_ch = MinSequence()

                    for ul_ch in ul.traverse():
                        min_vl_ul_ch.update((min_costs[(vl, ul_ch)].min, ul_ch))
                        min_vr_ul_ch.update((min_costs[(vr, ul_ch)].min, ul_ch))

                    for ur_ch in ur.traverse():
                        min_vl_ur_ch.update((min_costs[(vl, ur_ch)].min, ur_ch))
                        min_vr_ur_ch.update((min_costs[(vr, ur_ch)].min, ur_ch))

                    # Map left gene to left species and right gene to right
                    for ul_ch, ur_ch in product(min_vl_ul_ch, min_vr_ur_ch):
                        options.update((
                            min_vl_ul_ch.min + min_vr_ur_ch.min
                            + loss_cost * (
                                species_lca.distance(u, ul_ch)
                                + species_lca.distance(u, ur_ch) - 2
                            ),
                            MappingInfo(
                                species=u,
                                left=ul_ch,
                                right=ur_ch,
                            )
                        ))

                    # Map right gene to left species and left gene to right
                    for ul_ch, ur_ch in product(min_vr_ul_ch, min_vl_ur_ch):
                        options.update((
                            min_vr_ul_ch.min + min_vl_ur_ch.min
                            + loss_cost * (
                                species_lca.distance(u, ul_ch)
                                + species_lca.distance(u, ur_ch) - 2
                            ),
                            MappingInfo(
                                species=u,
                                left=ur_ch,
                                right=ul_ch,
                            )
                        ))

                min_vl_u_ch = MinSequence()
                min_vl_u_inc = MinSequence()
                min_vr_u_ch = MinSequence()
                min_vr_u_inc = MinSequence()

                for w in species_tree.traverse():
                    if species_lca.is_ancestor_of(u, w):
                        min_vl_u_ch.update((min_costs[(vl, w)].min, w))
                        min_vr_u_ch.update((min_costs[(vr, w)].min, w))
                    elif not species_lca.is_ancestor_of(w, u):
                        min_vl_u_inc.update((min_costs[(vl, w)].min, w))
                        min_vr_u_inc.update((min_costs[(vr, w)].min, w))

                # Try mapping as a duplication
                for l_u_ch, r_u_ch in product(min_vl_u_ch, min_vr_u_ch):
                    options.update((
                        dup_cost + min_vl_u_ch.min + min_vr_u_ch.min
                        + loss_cost * (
                            species_lca.distance(u, l_u_ch)
                            + species_lca.distance(u, r_u_ch)
                        ),
                        MappingInfo(
                            species=u,
                            left=l_u_ch,
                            right=r_u_ch,
                        )
                    ))

                # Try mapping as a horizontal gene transfer
                # Map left gene as a transfer and right gene as conserved
                for l_u_inc, r_u_ch in product(min_vl_u_inc, min_vr_u_ch):
                    options.update((
                        hgt_cost + min_vl_u_inc.min + min_vr_u_ch.min
                        + loss_cost * species_lca.distance(u, r_u_ch),
                        MappingInfo(
                            species=u,
                            left=l_u_inc,
                            right=r_u_ch,
                        )
                    ))

                # Map left gene as conserved and right gene as a transfer
                for l_u_ch, r_u_inc in product(min_vl_u_ch, min_vr_u_inc):
                    options.update((
                        hgt_cost + min_vl_u_ch.min + min_vr_u_inc.min
                        + loss_cost * species_lca.distance(u, l_u_ch),
                        MappingInfo(
                            species=u,
                            left=l_u_ch,
                            right=r_u_inc,
                        )
                    ))

    return min_costs


def _decode_thl_table(root, solutions, min_costs):
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
                    min_costs[(left, solution.left)],
                    min_costs,
                ),
                _decode_thl_table(
                    right,
                    min_costs[(right, solution.right)],
                    min_costs,
                ),
            )

            for map_left, map_right in mappings:
                results.append({
                    root: solution.species,
                    **map_left, **map_right,
                })

    return results


def reconcile_lca(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
) -> Reconciliation:
    """
    Compute the minimum-cost reconciliation between species tree and a gene
    tree, allowing for duplications and losses and where a duplication
    costs the same as loss.

    :param gene_tree: gene tree to reconcile
    :param species_tree: species tree to reconcile
    :returns: the reconciliation result
    """
    rec = {}
    species_lca = LowestCommonAncestor(species_tree)

    for gene in gene_tree.traverse("postorder"):
        if gene.is_leaf():
            species = species_tree.get_leaves_by_name(gene.species)[0]
            rec[gene] = species
        else:
            left, right = gene.children
            rec[gene] = species_lca(rec[left], rec[right])

    return rec


def reconcile_thl(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
    costs: CostVector
) -> Tuple[ExtendedIntegral, List[Reconciliation]]:
    """
    Find all minimum-cost reconciliations between a species tree and a gene
    tree, allowing for duplications, horizontal gene transfers and losses.

    If you wish to disallow horizontal gene transfers, using
    `reconcile_lca` is faster than setting `costs[HorizontalGeneTransfer]`
    to infinity.

    :param gene_tree: gene tree to reconcile
    :param species_tree: species tree to reconcile
    :param costs: cost of each evolutionary event
    :returns: minimum cost of each reconciliation and list of reconciliations
        with such costs
    """
    min_costs = _compute_thl_table(gene_tree, species_tree, costs or {})
    solutions = MinSequence()

    for u in species_tree.traverse("postorder"):
        for option in min_costs[(gene_tree, u)]:
            solutions.update((min_costs[(gene_tree, u)].min, option))

    return solutions.min, _decode_thl_table(gene_tree, solutions, min_costs)
