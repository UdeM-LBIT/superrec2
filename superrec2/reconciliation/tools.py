from enum import Enum, auto
from typing import Mapping
from ete3 import PhyloTree, PhyloNode
from infinity import inf
from ..utils.lowest_common_ancestor import LowestCommonAncestor
from ..utils.min_sequence import ExtendedIntegral


Reconciliation = Mapping[PhyloNode, PhyloNode]


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
        return inf

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


def parse_reconciliation(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
    source: str
) -> Reconciliation:
    """
    Parse a string representation of a reconciliation.

    :param gene_tree: gene tree of the reconciliation
    :param species_tree: species tree of the reconciliation
    :param source: string to parse
    :returns: parsed reconciliation
    """
    result: Reconciliation = {}

    for pair in source.split(","):
        if pair.strip():
            gene, species = pair.split(":")
            result[gene_tree & gene.strip()] = species_tree & species.strip()

    return result
