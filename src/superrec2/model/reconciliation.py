"""Represent and parse reconciliation problems and results."""
from dataclasses import dataclass, field
from enum import Enum, auto
from itertools import chain, product
from textwrap import indent
from typing import Generator, Mapping, Type, TypeVar, Union
from ete3 import Tree, TreeNode
from infinity import inf, Infinity
from .tree_mapping import (
    TreeMapping,
    get_species_mapping,
    parse_tree_mapping,
    serialize_tree_mapping,
)
from .synteny import (
    SyntenyMapping,
    parse_synteny_mapping,
    serialize_synteny_mapping,
)
from ..utils.subsequences import (
    subseq_complete,
    mask_from_subseq,
    subseq_segment_dist,
)
from ..utils.trees import LowestCommonAncestor, is_binary, binarize


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


Self = TypeVar("Self", bound="ReconciliationInput")


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

    def to_dict(self):
        """Convert the problem to a plain dictionary."""
        return {
            "object_tree": self.object_tree.write(
                format=8,
                format_root_node=True,
                features=["color"],
            ),
            "species_tree": self.species_lca.tree.write(
                format=8,
                format_root_node=True,
                features=["color"],
            ),
            "leaf_object_species": serialize_tree_mapping(self.leaf_object_species),
            "costs": dict(((event.name, value) for event, value in self.costs.items())),
        }

    def __repr__(self):
        params = ",\n".join(
            f"{key}={repr(value)}" for key, value in self.to_dict().items()
        )
        return f'{self.__class__.__name__}(\n{indent(params, " " * 2)}\n)'

    @classmethod
    def _from_dict(cls, data):
        object_tree = Tree(data["object_tree"], format=1)
        species_tree = Tree(data["species_tree"], format=1)

        if "leaf_object_species" in data:
            leaf_object_species = parse_tree_mapping(
                object_tree, species_tree, data["leaf_object_species"]
            )
        else:
            leaf_object_species = get_species_mapping(
                object_tree,
                species_tree,
            )

        if "costs" in data:
            costs = {}

            for event, value in data["costs"].items():
                if isinstance(event, (NodeEvent, EdgeEvent)):
                    event_enum = event
                else:
                    if hasattr(NodeEvent, event):
                        event_enum = getattr(NodeEvent, event)
                    else:
                        event_enum = getattr(EdgeEvent, event)

                costs[event_enum] = value
        else:
            costs = get_default_cost()

        return {
            "object_tree": object_tree,
            "species_lca": LowestCommonAncestor(species_tree),
            "leaf_object_species": leaf_object_species,
            "costs": costs,
        }

    @classmethod
    def from_dict(cls: Type[Self], data) -> Self:
        """Reconstruct a problem from its plain dictionary representation."""
        return cls(**cls._from_dict(data))

    def binarize(self: Self) -> Generator[Self, None, None]:
        """
        Generate all possible inputs that can be derived from this one
        by resolving polytomies in the object and species trees. If the
        input is already binary, it is returned unchanged.
        """
        if is_binary(self.object_tree) and is_binary(self.species_lca.tree):
            yield self
            return

        for object_tree, species_tree in product(
            binarize(self.object_tree),
            binarize(self.species_lca.tree),
        ):
            result = self.__class__.from_dict(
                {
                    **self.to_dict(),
                    "object_tree": object_tree.write(
                        format=8,
                        format_root_node=True,
                        features=["color"],
                    ),
                    "species_tree": species_tree.write(
                        format=8,
                        format_root_node=True,
                        features=["color"],
                    ),
                }
            )
            yield result

    def label_internal(self) -> None:
        """
        Ensure that all internal nodes have a label by automatically
        generating new ones for unlabeled nodes.
        """
        next_object = 0
        next_species = 0

        for node_object in self.object_tree.traverse("preorder"):
            if not node_object.name or node_object.name == "NoName":
                while f"O{next_object}" in self.object_tree:
                    next_object += 1

                node_object.name = f"O{next_object}"

        for node_species in self.species_lca.tree.traverse("preorder"):
            if not node_species.name or node_species.name == "NoName":
                while f"S{next_species}" in self.species_lca.tree:
                    next_species += 1

                node_species.name = f"S{next_species}"

    def __hash__(self):
        return hash(
            (
                self.object_tree,
                self.species_lca,
                tuple(
                    sorted(serialize_tree_mapping(self.leaf_object_species).items()),
                ),
                tuple(
                    (event, self.costs.get(event))
                    for event in chain(
                        NodeEvent.__members__.keys(),
                        EdgeEvent.__members__.keys(),
                    )
                ),
            )
        )


@dataclass(frozen=True)
class ReconciliationOutput:
    """Output for the reconciliation problem."""

    # Problem input
    input: ReconciliationInput

    # Mapping of the object tree onto the species tree
    object_species: TreeMapping

    def to_dict(self):
        """Convert the output to a plain dictionary."""
        return {
            "input": self.input.to_dict(),
            "object_species": serialize_tree_mapping(self.object_species),
        }

    def __repr__(self):
        keyvals = self.to_dict()
        keyvals["input"] = self.input
        params = ",\n".join(f"{key}={repr(value)}" for key, value in keyvals.items())
        return f'{self.__class__.__name__}(\n{indent(params, " " * 2)}\n)'

    @classmethod
    def _from_dict(cls, data):
        input_problem = ReconciliationInput.from_dict(data["input"])
        return {
            "input": input_problem,
            "object_species": parse_tree_mapping(
                input_problem.object_tree,
                input_problem.species_lca.tree,
                data["object_species"],
            ),
        }

    @classmethod
    def from_dict(cls, data):
        """Reconstruct an output from its plain dictionary representation."""
        return cls(**cls._from_dict(data))

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
                costs[NodeEvent.SPECIATION]
                + left_cost
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
        return hash(
            (
                self.input,
                tuple(sorted(serialize_tree_mapping(self.object_species).items())),
            )
        )


@dataclass(frozen=True, repr=False)
class SuperReconciliationInput(ReconciliationInput):
    """Input to the super-reconciliation problem."""

    # Extant syntenies at the leaves of the synteny tree
    leaf_syntenies: SyntenyMapping = field(default_factory=dict)

    def to_dict(self):
        return {
            **super().to_dict(),
            "leaf_syntenies": serialize_synteny_mapping(self.leaf_syntenies),
        }

    @classmethod
    def _from_dict(cls, data):
        parent = super()._from_dict(data)
        return {
            **parent,
            "leaf_syntenies": parse_synteny_mapping(
                parent["object_tree"],
                data["leaf_syntenies"],
            ),
        }

    def __hash__(self):
        return hash(
            (
                super().__hash__(),
                tuple(
                    sorted(
                        (
                            (node, tuple(synteny))
                            for node, synteny in serialize_synteny_mapping(
                                self.leaf_syntenies
                            ).items()
                        )
                    )
                ),
            )
        )


@dataclass(frozen=True, repr=False)
class SuperReconciliationOutput(ReconciliationOutput):
    """Output for the reconciliation problem."""

    # Problem input
    input: SuperReconciliationInput

    # Mapping of internal nodes to syntenies
    syntenies: SyntenyMapping

    # Whether the syntenies are to be considered as ordered sequences
    # or unordered sets
    ordered: bool

    def to_dict(self):
        return {
            **super().to_dict(),
            "syntenies": serialize_synteny_mapping(self.syntenies),
            "ordered": self.ordered,
        }

    @classmethod
    def _from_dict(cls, data):
        parent = super()._from_dict(data)
        return {
            **parent,
            "syntenies": parse_synteny_mapping(
                parent["input"].object_tree,
                data["syntenies"],
            ),
            "ordered": data.get("ordered", True),
        }

    def reconciliation_cost(self):
        """Compute the cost of the reconciliation part."""
        return super().cost()

    def _ordered_labeling_cost(self):
        """Compute the ordered segmental loss cost of the labeling."""
        tree = self.input.object_tree
        rec = self.object_species

        total_cost = 0
        root_syn = self.syntenies[tree]
        masks = {tree: subseq_complete(root_syn)}

        sloss_cost = self.input.costs[EdgeEvent.SEGMENTAL_LOSS]

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
                    total_cost += (
                        subseq_segment_dist(left_mask, sub_mask, True)
                        + subseq_segment_dist(right_mask, sub_mask, True)
                    ) * sloss_cost
                elif event == NodeEvent.DUPLICATION:
                    total_cost += (
                        min(
                            subseq_segment_dist(left_mask, sub_mask, True)
                            + subseq_segment_dist(right_mask, sub_mask, False),
                            subseq_segment_dist(left_mask, sub_mask, False)
                            + subseq_segment_dist(right_mask, sub_mask, True),
                        )
                        * sloss_cost
                    )
                else:
                    assert event == NodeEvent.HORIZONTAL_TRANSFER
                    keep_left = self.input.species_lca.is_comparable(
                        rec[node], rec[left_node]
                    )
                    total_cost += (
                        subseq_segment_dist(left_mask, sub_mask, keep_left)
                        + subseq_segment_dist(right_mask, sub_mask, not keep_left)
                    ) * sloss_cost

        return total_cost

    def _unordered_labeling_cost(self):
        """Compute the unordered segmental loss cost of the labeling."""
        tree = self.input.object_tree
        rec = self.object_species

        total_cost = 0
        sloss_cost = self.input.costs[EdgeEvent.SEGMENTAL_LOSS]

        for node in tree.traverse("preorder"):
            if not node.is_leaf():
                event = self.node_event(node)
                left_node, right_node = node.children

                node_set = set(self.syntenies[node])
                left_cost = (
                    0 if node_set <= set(self.syntenies[left_node]) else sloss_cost
                )
                right_cost = (
                    0 if node_set <= set(self.syntenies[right_node]) else sloss_cost
                )

                if event == NodeEvent.SPECIATION:
                    total_cost += left_cost + right_cost
                elif event == NodeEvent.DUPLICATION:
                    total_cost += min(left_cost, right_cost)
                else:
                    assert event == NodeEvent.HORIZONTAL_TRANSFER
                    if self.input.species_lca.is_comparable(rec[node], rec[left_node]):
                        total_cost += left_cost
                    else:
                        total_cost += right_cost

        return total_cost

    def labeling_cost(self):
        """Compute the segmental loss cost of the labeling."""
        if self.ordered:
            return self._ordered_labeling_cost()

        return self._unordered_labeling_cost()

    def cost(self):
        """Compute the cost of this super-reconciliation."""
        return self.reconciliation_cost() + self.labeling_cost()

    def __hash__(self):
        return hash(
            (
                super().__hash__(),
                tuple(
                    sorted(
                        (
                            (node, tuple(synteny))
                            for node, synteny in serialize_synteny_mapping(
                                self.syntenies
                            ).items()
                        )
                    )
                ),
                self.ordered,
            )
        )
