"""Reconciliation and super-reconciliation representation."""
from enum import Enum, auto
from itertools import product
from typing import Any, Dict, Generator, List, Mapping, NamedTuple, Sequence
from ete3 import PhyloTree, PhyloNode
from infinity import inf
from ..utils.lowest_common_ancestor import LowestCommonAncestor
from ..utils.min_sequence import ExtendedIntegral
from ..utils.subsequences import (
    subseq_complete,
    mask_from_subseq,
    subseq_segment_dist,
)


Reconciliation = Mapping[PhyloNode, PhyloNode]
Labeling = Mapping[PhyloNode, Sequence[Any]]


class SuperReconciliationInput(NamedTuple):
    """Input for the super-reconciliation problem."""

    # Tree of syntenies to map onto the species tree
    synteny_tree: PhyloTree

    # Ancestry information of the species to reconcile
    species_lca: LowestCommonAncestor

    # Extant syntenies at the leaves of the synteny tree
    leaf_labeling: Labeling

    def __repr__(self):
        synteny_tree = self.synteny_tree.write(format=8, format_root_node=True)
        species_tree = self.species_lca.tree.write(
            format=8, format_root_node=True
        )
        leaf_labeling = serialize_labeling(self.leaf_labeling)
        return f"""SuperReconciliationInput({", ".join([
            f'synteny_tree="{synteny_tree}", '
            f'species_tree="{species_tree}", '
            f'leaf_labeling="{leaf_labeling}"'
        ])})"""


class SuperReconciliation(NamedTuple):
    """Super-reconciliation of a synteny tree with a species tree."""

    # Input to the super-reconciliation problem
    input: SuperReconciliationInput

    # Mapping of the synteny tree onto the species tree
    reconciliation: Reconciliation

    # Labeling of ancestral syntenies
    labeling: Labeling

    def __repr__(self):
        reconciliation = serialize_reconciliation(self.reconciliation)
        labeling = serialize_labeling(self.labeling)
        return f"""SuperReconciliation({", ".join([
            f"input={repr(self.input)}",
            f'reconciliation="{reconciliation}"',
            f'labeling="{labeling}"',
        ])})"""


def get_species_name(gene_name):
    """
    Extract the species name out of a gene name.

    >>> get_species_name("x_21")
    "X"

    :param gene_name: full gene name
    :returns: first part of the gene name in uppercase
    """
    return gene_name.split("_")[0].upper()


class Event(Enum):
    """Evolutionary events affecting genes."""

    # Sentinel for extant (current) genes
    LEAF = auto()

    # Sentinel for scenarios that are invalid wrt our evolutionary model
    INVALID = auto()

    # Transmission of the parent gene to both children species
    SPECIATION = auto()

    # Duplication of the parent gene in the same genome
    DUPLICATION = auto()

    # Transfer of the parent gene to a foreign genome
    HORIZONTAL_GENE_TRANSFER = auto()


def get_event(
    root_gene: PhyloNode,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
) -> Event:
    """
    Find the event associated to a gene by a reconciliation.

    :param root_gene: gene to query
    :param species_lca: species ancestry information
    :param rec: reconciliation to use
    :returns: event associated to the given gene node
    """
    if root_gene.is_leaf():
        return (
            Event.LEAF
            if rec[root_gene].name == root_gene.species
            else Event.INVALID
        )

    left_gene, right_gene = root_gene.children

    if species_lca.is_strict_ancestor_of(
        rec[left_gene], rec[root_gene]
    ) or species_lca.is_strict_ancestor_of(rec[right_gene], rec[root_gene]):
        return Event.INVALID

    if species_lca.is_ancestor_of(
        rec[root_gene], rec[left_gene]
    ) and species_lca.is_ancestor_of(rec[root_gene], rec[right_gene]):
        return (
            Event.SPECIATION
            if (
                rec[root_gene] == species_lca(rec[left_gene], rec[right_gene])
                and not species_lca.is_comparable(
                    rec[left_gene], rec[right_gene]
                )
            )
            else Event.DUPLICATION
        )

    if species_lca.is_ancestor_of(
        rec[root_gene], rec[left_gene]
    ) or species_lca.is_ancestor_of(rec[root_gene], rec[right_gene]):
        return Event.HORIZONTAL_GENE_TRANSFER

    return Event.INVALID


class CostType(Enum):
    """Evolutionary events to which a cost can be assigned."""

    DUPLICATION = auto()
    HORIZONTAL_GENE_TRANSFER = auto()
    FULL_LOSS = auto()
    SEGMENTAL_LOSS = auto()


# Cost values for each possible evolutionary event
CostVector = Mapping[CostType, ExtendedIntegral]


def get_reconciliation_cost(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
    costs: CostVector,
) -> ExtendedIntegral:
    """
    Compute the total cost of a reconciliation.

    :param gene_tree: reconciled gene tree
    :param species_lca: ancestry information about the reconciled species tree
    :param rec: reconciliation to use
    :param costs: costs for each evolutionary event
    """
    event = get_event(gene_tree, species_lca, rec)

    if event == Event.INVALID:
        return inf

    if event == Event.LEAF:
        return 0

    left_gene, right_gene = gene_tree.children
    cost_vl = get_reconciliation_cost(left_gene, species_lca, rec, costs)
    dist_vl = species_lca.distance(rec[gene_tree], rec[left_gene])
    cost_vr = get_reconciliation_cost(right_gene, species_lca, rec, costs)
    dist_vr = species_lca.distance(rec[gene_tree], rec[right_gene])

    if event == Event.SPECIATION:
        return (
            cost_vl
            + cost_vr
            + costs[CostType.FULL_LOSS] * (dist_vl + dist_vr - 2)
        )

    if event == Event.DUPLICATION:
        return (
            costs[CostType.DUPLICATION]
            + cost_vl
            + cost_vr
            + costs[CostType.FULL_LOSS] * (dist_vl + dist_vr)
        )

    assert event == Event.HORIZONTAL_GENE_TRANSFER

    dist_conserved = (
        dist_vl
        if species_lca.is_ancestor_of(rec[gene_tree], rec[left_gene])
        else dist_vr
    )
    return (
        costs[CostType.HORIZONTAL_GENE_TRANSFER]
        + cost_vl
        + cost_vr
        + costs[CostType.FULL_LOSS] * dist_conserved
    )


def get_labeling_cost(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
    labeling: Labeling,
) -> int:
    """
    Compute the cost of a synteny labeling.

    :param gene_tree: reconciled gene tree
    :param species_lca: ancestry information about the reconciled species tree
    :param rec: reconciliation to use
    :param labeling: synteny labeling
    """
    total_cost = 0
    root_syn = labeling[gene_tree]
    masks = {gene_tree: subseq_complete(root_syn)}

    for sub_gene in gene_tree.traverse("preorder"):
        if not sub_gene.is_leaf():
            event = get_event(sub_gene, species_lca, rec)
            sub_mask = masks[sub_gene]
            left_gene, right_gene = sub_gene.children

            left_mask = masks[left_gene] = mask_from_subseq(
                labeling[left_gene], root_syn
            )

            right_mask = masks[right_gene] = mask_from_subseq(
                labeling[right_gene], root_syn
            )

            if event == Event.SPECIATION:
                total_cost += subseq_segment_dist(
                    left_mask, sub_mask, True
                ) + subseq_segment_dist(right_mask, sub_mask, True)
            elif event == Event.DUPLICATION:
                total_cost += min(
                    (
                        subseq_segment_dist(left_mask, sub_mask, True)
                        + subseq_segment_dist(right_mask, sub_mask, False)
                    ),
                    (
                        subseq_segment_dist(left_mask, sub_mask, False)
                        + subseq_segment_dist(right_mask, sub_mask, True)
                    ),
                )
            elif event == Event.HORIZONTAL_GENE_TRANSFER:
                keep_left = species_lca.is_comparable(
                    rec[sub_gene], rec[left_gene]
                )
                total_cost += subseq_segment_dist(
                    left_mask, sub_mask, keep_left
                ) + subseq_segment_dist(right_mask, sub_mask, not keep_left)
            else:
                raise ValueError("Invalid event")

    return total_cost


def reconcile_leaves(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
) -> Reconciliation:
    """
    Compute a partial reconciliation containing forced leaves assignments.

    :param gene_tree: gene tree to reconcile
    :param species_tree: species tree to reconcile
    :returns: mapping of the leaf genes onto leaf species
    """
    return {
        gene_leaf: species_tree & gene_leaf.species
        for gene_leaf in gene_tree.get_leaves()
    }


def reconcile_all(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
) -> Generator[Reconciliation, None, None]:
    """
    Generate all valid reconciliations.

    :param gene_tree: gene tree to reconcile
    :param species_lca: species ancestry information
    :returns: all valid mappings of the gene tree onto the species tree
    """
    if gene_tree.is_leaf():
        yield {gene_tree: species_lca.tree & gene_tree.species}
        return

    left_gene, right_gene = gene_tree.children

    for map_left, map_right in product(
        reconcile_all(left_gene, species_lca),
        reconcile_all(right_gene, species_lca),
    ):
        left_species = map_left[left_gene]
        right_species = map_right[right_gene]
        lca = species_lca(left_species, right_species)

        parent_species = lca
        while parent_species is not None:
            yield {
                gene_tree: parent_species,
                **map_left,
                **map_right,
            }
            parent_species = parent_species.up

        for (transfer_target, other_target) in (
            (left_species, right_species),
            (right_species, left_species),
        ):
            if species_lca.is_ancestor_of(other_target, transfer_target):
                continue

            transfer_species = transfer_target

            while transfer_species != lca:
                yield {
                    gene_tree: transfer_species,
                    **map_left,
                    **map_right,
                }
                transfer_species = transfer_species.up


def parse_reconciliation(
    gene_tree: PhyloTree, species_tree: PhyloTree, source: str
) -> Reconciliation:
    """
    Parse a string representation of a reconciliation.

    :param gene_tree: gene tree of the reconciliation
    :param species_tree: species tree of the reconciliation
    :param source: string to parse
    :returns: parsed reconciliation
    """
    result: Dict[PhyloNode, PhyloNode] = {}

    for pair in source.split(","):
        if pair.strip():
            gene, species = pair.split(":")
            result[gene_tree & gene.strip()] = species_tree & species.strip()

    return result


def serialize_reconciliation(reconciliation: Reconciliation) -> str:
    """
    Serialize a reconciliation.

    :param reconciliation: reconciliation to serialize
    :returns: serialized representation
    """
    return ",".join(
        f"{gene.name}:{species.name}"
        for gene, species in reconciliation.items()
    )


def parse_labeling(gene_tree: PhyloTree, source: str) -> Labeling:
    """
    Parse a string representation of a synteny labeling.

    :param gene_tree: labeled gene tree
    :param source: string to parse
    :returns: parsed labeling
    """
    result: Dict[PhyloNode, List[str]] = {}

    for pair in source.split(","):
        if pair.strip():
            node, synteny = pair.split(":")
            result[gene_tree & node.strip()] = list(synteny.strip())

    return result


def serialize_labeling(labeling: Labeling) -> str:
    """
    Serialize a synteny labeling.

    :param labeling: labeling to serialize
    :returns: serialized representation
    """
    return ",".join(
        f"{node.name}:{''.join(synteny)}" for node, synteny in labeling.items()
    )
