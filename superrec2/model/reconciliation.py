from typing import Mapping, Union
from dataclasses import dataclass, field
from enum import Enum, auto
from ete3 import Tree, TreeNode
from textwrap import indent
from infinity import inf, Infinity
from .tree_mapping import TreeMapping, serialize_tree_mapping
from .synteny import SyntenyMapping, serialize_synteny_mapping
from ..utils.subsequences import (
    subseq_complete,
    mask_from_subseq,
    subseq_segment_dist,
)
from ..utils.trees import LowestCommonAncestor


class NodeEvent(Enum):
    """Evolutionary events that can happen at a node of the object tree."""
    # Sentinel for extant (current) object
    LEAF = auto()

    # Sentinel for scenarios that are invalid wrt the evolutionary model
    INVALID = auto()

    # Transmission of the parent object to both children species
    SPECIATION = auto()

    # Duplication of the parent object in the same genome
    DUPLICATION = auto()

    # Transfer of the parent object to a foreign genome
    HORIZONTAL_TRANSFER = auto()


class EdgeEvent(Enum):
    """Evolutionary events that can happen on an edge of the object tree."""
    # Complete loss of the object in a parallel lineage
    FULL_LOSS = auto()

    # Loss of a segment of the synteny
    SEGMENTAL_LOSS = auto()


# Any kind of evolutionary event
Event = Union[NodeEvent, EdgeEvent]


# Cost values for evolutionary events
CostValues = Mapping[Event, Union[int, Infinity]]


def get_default_cost() -> CostValues:
    """Get the default event cost vector."""
    return {
        NodeEvent.SPECIATION: 0,
        NodeEvent.DUPLICATION: 1,
        NodeEvent.HORIZONTAL_TRANSFER: 1,
        EdgeEvent.FULL_LOSS: 1,
        EdgeEvent.SEGMENTAL_LOSS: 1,
    }


@dataclass(frozen=True)
class ReconciliationInput:
    """Input to the reconciliation problem."""

    # Tree of objects to map onto the species tree
    object_tree: Tree

    # Ancestry information of the species to reconcile
    species_lca: LowestCommonAncestor

    # Mapping of the leaf objects onto the species tree
    leaf_object_species: TreeMapping

    # Costs of evolutionary events
    costs: CostValues = field(default_factory=get_default_cost)

    def _repr(self):
        object_tree = self.object_tree.write(format=8, format_root_node=True)
        species_tree = self.species_lca.tree.write(
            format=8, format_root_node=True
        )
        leaf_object_species = serialize_tree_mapping(self.leaf_object_species)
        return ",\n".join((
            f'object_tree="{object_tree}"',
            f'species_tree="{species_tree}"',
            f'leaf_object_species="{leaf_object_species}"',
            f'costs={self.costs}',
        ))

    def __repr__(self):
        return f'{self.__class__.__name__}(\n{indent(self._repr(), "    ")}\n)'

    def __hash__(self):
        return hash((
            self.object_tree,
            self.species_lca,
            serialize_tree_mapping(self.leaf_object_species),
            tuple(self.costs.items()),
        ))


@dataclass(frozen=True)
class ReconciliationOutput:
    """Output for the reconciliation problem."""

    # Problem input
    input: ReconciliationInput

    # Mapping of the object tree onto the species tree
    object_species: TreeMapping

    def _repr(self):
        object_species = serialize_tree_mapping(self.object_species)
        return ",\n".join((
            f'input={repr(self.input)}',
            f'object_species="{object_species}"',
        ))

    def __repr__(self):
        return f'{self.__class__.__name__}(\n{indent(self._repr(), "    ")}\n)'

    def node_event(self, node: TreeNode) -> NodeEvent:
        """
        Find the event associated to a node.

        :param node: node to query
        :returns: event associated to the node
        """
        species_lca = self.input.species_lca
        rec = self.object_species

        if node.is_leaf():
            return (
                NodeEvent.LEAF
                if rec[node] == self.input.leaf_object_species[node]
                else NodeEvent.INVALID
            )

        left_node, right_node = node.children

        if species_lca.is_strict_ancestor_of(
            rec[left_node], rec[node]
        ) or species_lca.is_strict_ancestor_of(rec[right_node], rec[node]):
            return NodeEvent.INVALID

        if species_lca.is_ancestor_of(
            rec[node], rec[left_node]
        ) and species_lca.is_ancestor_of(rec[node], rec[right_node]):
            return (
                NodeEvent.SPECIATION
                if (
                    rec[node] == species_lca(rec[left_node], rec[right_node])
                    and not species_lca.is_comparable(
                        rec[left_node],
                        rec[right_node],
                    )
                )
                else NodeEvent.DUPLICATION
            )

        if species_lca.is_ancestor_of(
            rec[node], rec[left_node]
        ) or species_lca.is_ancestor_of(rec[node], rec[right_node]):
            return NodeEvent.HORIZONTAL_TRANSFER

        return NodeEvent.INVALID

    def _cost_rec(self, node: TreeNode = None) -> Union[int, Infinity]:
        event = self.node_event(node)
        species_lca = self.input.species_lca
        costs = self.input.costs
        rec = self.object_species

        if event == NodeEvent.INVALID:
            return inf

        if event == NodeEvent.LEAF:
            return 0

        left_node, right_node = node.children
        left_cost = self._cost_rec(left_node)
        left_dist = species_lca.distance(rec[node], rec[left_node])
        right_cost = self._cost_rec(right_node)
        right_dist = species_lca.distance(rec[node], rec[right_node])

        if event == NodeEvent.SPECIATION:
            return (
                left_cost
                + right_cost
                + costs[EdgeEvent.FULL_LOSS] * (left_dist + right_dist - 2)
            )

        if event == NodeEvent.DUPLICATION:
            return (
                costs[NodeEvent.DUPLICATION]
                + left_cost
                + right_cost
                + costs[EdgeEvent.FULL_LOSS] * (left_dist + right_dist)
            )

        assert event == NodeEvent.HORIZONTAL_TRANSFER

        dist_conserved = (
            left_dist
            if species_lca.is_ancestor_of(rec[node], rec[left_node])
            else right_dist
        )
        return (
            costs[NodeEvent.HORIZONTAL_TRANSFER]
            + left_cost
            + right_cost
            + costs[EdgeEvent.FULL_LOSS] * dist_conserved
        )

    def cost(self) -> Union[int, Infinity]:
        """Compute the total cost of this reconciliation."""
        return self._cost_rec(self.input.object_tree)

    def __hash__(self):
        return hash((
            self.input,
            tuple(self.object_species.items()),
        ))


@dataclass(frozen=True)
class SuperReconciliationInput(ReconciliationInput):
    """Input to the super-reconciliation problem."""

    # Extant syntenies at the leaves of the synteny tree
    leaf_syntenies: SyntenyMapping = field(default_factory = dict)

    def _repr(self):
        leaf_syntenies = serialize_synteny_mapping(self.leaf_syntenies)
        return "\n".join((
            super()._repr(),
            f'leaf_syntenies="{leaf_syntenies}"'
        ))

    def __repr__(self):
        return super().__repr__()

    def __hash__(self):
        return hash((
            super().__hash__(),
            serialize_synteny_mapping(self.leaf_syntenies),
        ))


@dataclass(frozen=True)
class SuperReconciliationOutput(ReconciliationOutput):
    """Output for the reconciliation problem."""

    # Problem input
    input: SuperReconciliationInput

    # Mapping of internal nodes to syntenies
    syntenies: SyntenyMapping

    def _repr(self):
        syntenies = serialize_synteny_mapping(self.syntenies)
        return "\n".join((
            super()._repr(),
            f'syntenies="{syntenies}"'
        ))

    def __repr__(self):
        return super().__repr__()

    def reconciliation_cost(self):
        """Compute the cost of the reconciliation part."""
        return super().cost()

    def labeling_cost(self):
        """Compute the cost of the labeling part."""
        tree = self.input.object_tree
        rec = self.object_species

        total_cost = 0
        root_syn = self.syntenies[tree]
        masks = {tree: subseq_complete(root_syn)}

        for node in tree.traverse("preorder"):
            if not node.is_leaf():
                event = self.node_event(node)
                sub_mask = masks[node]
                left_node, right_node = node.children

                left_mask = masks[left_node] = mask_from_subseq(
                    self.syntenies[left_node], root_syn
                )

                right_mask = masks[right_node] = mask_from_subseq(
                    self.syntenies[right_node], root_syn
                )

                if event == NodeEvent.SPECIATION:
                    total_cost += subseq_segment_dist(
                        left_mask, sub_mask, True
                    ) + subseq_segment_dist(right_mask, sub_mask, True)
                elif event == NodeEvent.DUPLICATION:
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
                else:
                    assert event == NodeEvent.HORIZONTAL_TRANSFER
                    keep_left = self.input.species_lca.is_comparable(
                        rec[node], rec[left_node]
                    )
                    total_cost += subseq_segment_dist(
                        left_mask, sub_mask, keep_left
                    ) + subseq_segment_dist(right_mask, sub_mask, not keep_left)

        return total_cost


    def cost(self):
        """Compute the total cost of this super-reconciliation."""
        return self.reconciliation_cost() + self.labeling_cost()

    def __hash__(self):
        return hash((
            super().__hash__(),
            serialize_synteny_mapping(self.syntenies),
        ))
