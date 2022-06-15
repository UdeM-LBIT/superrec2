"""Represent and parse reconciliation problems and results."""
from dataclasses import dataclass, field
from enum import Enum, auto
from itertools import chain, product
from textwrap import indent
from typing import (
    Generator,
    List,
    Mapping,
    NamedTuple,
    Optional,
    Type,
    TypeVar,
    Union,
)
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
from ..utils.dynamic_programming import (
    Entry,
    Candidate,
    RetentionPolicy,
    MergePolicy,
)


class NodeEvent(Enum):
    """Evolutionary events that can happen at a node of the object tree."""

    NONE = auto()

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

    # Transfer to a different lineage and loss of the original copy
    TRANSFER_LOSS = auto()

    # Speciation to an unsampled lineage followed by a return transfer
    UNSAMPLED_TRANSFER = auto()


# Any kind of evolutionary event
Event = Union[NodeEvent, EdgeEvent]


class EventsAtNode(NamedTuple):
    """Events at a node and on its outgoing edges."""

    # Event at the node
    node_event: NodeEvent

    # When the node event is a duplication or a transfer, which child
    # is the target of that event
    target_child: Optional[TreeNode]

    # Event counts on the outgoing edges
    edge_events: List[Mapping[EdgeEvent, int]]

    def __hash__(self):
        return hash(
            (
                self.node_event,
                self.target_child,
                tuple(map(
                    lambda events: tuple(sorted(events.items())),
                    self.edge_events,
                )),
            )
        )


# Any number, or infinity
ExtendedInt = Union[int, Infinity]


class EventCosts:
    """Costs for the different evolutionary events."""

    def __init__(self, values: Mapping[Union[str, Event], ExtendedInt] = {}):
        self.values = {
            NodeEvent.SPECIATION: 0,
            NodeEvent.DUPLICATION: 1,
            NodeEvent.HORIZONTAL_TRANSFER: 1,
            EdgeEvent.FULL_LOSS: 1,
            EdgeEvent.SEGMENTAL_LOSS: 1,
            EdgeEvent.TRANSFER_LOSS: 2,
            EdgeEvent.UNSAMPLED_TRANSFER: 2,
        }

        for event, value in values.items():
            if isinstance(event, (NodeEvent, EdgeEvent)):
                event_enum = event
            else:
                if hasattr(NodeEvent, event):
                    event_enum = getattr(NodeEvent, event)
                else:
                    event_enum = getattr(EdgeEvent, event)

            self.values[event_enum] = value

    def to_dict(self):
        """Convert the set of costs to a dictionary."""
        return self.values

    def get(self, event: Event, count: ExtendedInt = 1) -> ExtendedInt:
        """Get the cost for the given number of events of the given type."""
        if count == 0:
            # No events incur no cost, even if the cost is infinite
            return 0

        return self.values[event] * count

    def __hash__(self):
        return hash(tuple(self.values.items()))


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
    costs: EventCosts = EventCosts()

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
            "leaf_object_species": serialize_tree_mapping(
                self.leaf_object_species
            ),
            "costs": self.costs.to_dict(),
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

        return {
            "object_tree": object_tree,
            "species_lca": LowestCommonAncestor(species_tree),
            "leaf_object_species": leaf_object_species,
            "costs": EventCosts(data.get("costs", {})),
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
                    sorted(
                        serialize_tree_mapping(self.leaf_object_species).items()
                    ),
                ),
                self.costs,
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
        params = ",\n".join(
            f"{key}={repr(value)}" for key, value in keyvals.items()
        )
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

    def events_at_node(self, node: TreeNode) -> EventsAtNode:
        """
        Find the events associated to a node and its outgoing edges.

        :param node: node to query
        :raises RuntimeError: if the reconciliation mapping is invalid
            around the given node
        :returns: events associated to the node
        """
        lca = self.input.species_lca
        species_node = self.object_species[node]

        if node.is_leaf():
            if species_node != self.input.leaf_object_species[node]:
                raise RuntimeError(
                    f"Invalid reconciliation: Leaf {node.name} is not "
                    "mapped to its host species"
                )

            return EventsAtNode(NodeEvent.NONE, None, {})

        for name, child in zip(("left", "right"), node.children):
            species_child = self.object_species[child]

            if lca.is_strict_ancestor_of(species_child, species_node):
                raise RuntimeError(
                    f"Invalid reconciliation: the {name} child of "
                    f"{node.name} is mapped to an ancestor of its species"
                )

        result: Entry[int, EventsAtNode] = Entry(
            MergePolicy.MIN,
            RetentionPolicy.ANY,
        )
        costs = self.input.costs

        def make_candidate(events: Mapping[Event, int]) -> Candidate:
            return Candidate(
                sum(costs.get(event, count) for event, count in events.items()),
                tuple(sorted(
                    (event, count)
                    for event, count in events.items()
                    if count != 0
                )),
            )

        for node_event in (
            NodeEvent.SPECIATION,
            NodeEvent.DUPLICATION,
            NodeEvent.HORIZONTAL_TRANSFER,
        ):
            targets = (
                (None,)
                if node_event == NodeEvent.SPECIATION
                else node.children
            )

            for target in targets:
                edges_events = []
                edges_cost = 0

                for child in node.children:
                    if (
                        species_node == species_child
                        and node_event == NodeEvent.SPECIATION
                    ):
                        continue

                    edge_entry = Entry(MergePolicy.MIN, RetentionPolicy.ANY)

                    species_child = self.object_species[child]
                    dist = lca.distance(species_node, species_child)

                    if lca.is_ancestor_of(species_node, species_child):
                        # Descending child
                        edge_entry.update(make_candidate({
                            EdgeEvent.UNSAMPLED_TRANSFER: 1,
                        }))

                        if node_event == NodeEvent.SPECIATION:
                            edge_entry.update(make_candidate({
                                EdgeEvent.FULL_LOSS: dist - 1,
                            }))
                        elif (
                            node_event == NodeEvent.DUPLICATION
                            or child != target
                        ):
                            edge_entry.update(make_candidate({
                                EdgeEvent.FULL_LOSS: dist,
                            }))
                    else:
                        # Separated child
                        if (
                            node_event == NodeEvent.HORIZONTAL_TRANSFER
                            and child == target
                        ):
                            edge_entry.update(make_candidate({}))
                        else:
                            edge_entry.update(make_candidate({
                                EdgeEvent.TRANSFER_LOSS: 1,
                            }))

                    edges_events.append(dict(edge_entry.info()))
                    edges_cost += edge_entry.value()

                if len(edges_events) == len(node.children):
                    result.update(Candidate(
                        costs.get(node_event, 1) + edges_cost,
                        EventsAtNode(
                            node_event,
                            target,
                            edges_events
                        )
                    ))

        return result.info()

    def _cost_rec(self, node: TreeNode = None) -> ExtendedInt:
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
                costs.get(NodeEvent.SPECIATION)
                + left_cost
                + right_cost
                + costs.get(EdgeEvent.FULL_LOSS, left_dist + right_dist - 2)
            )

        if event == NodeEvent.DUPLICATION:
            return (
                costs.get(NodeEvent.DUPLICATION)
                + left_cost
                + right_cost
                + costs.get(EdgeEvent.FULL_LOSS, left_dist + right_dist)
            )

        assert event == NodeEvent.HORIZONTAL_TRANSFER

        dist_conserved = (
            left_dist
            if species_lca.is_ancestor_of(rec[node], rec[left_node])
            else right_dist
        )
        return (
            costs.get(NodeEvent.HORIZONTAL_TRANSFER)
            + left_cost
            + right_cost
            + costs.get(EdgeEvent.FULL_LOSS, dist_conserved)
        )

    def cost(self) -> ExtendedInt:
        """Compute the total cost of this reconciliation."""
        return self._cost_rec(self.input.object_tree)

    def __hash__(self):
        return hash(
            (
                self.input,
                tuple(
                    sorted(serialize_tree_mapping(self.object_species).items())
                ),
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
        costs = self.input.costs
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
                    total_cost += costs.get(
                        EdgeEvent.SEGMENTAL_LOSS,
                        subseq_segment_dist(left_mask, sub_mask, True)
                        + subseq_segment_dist(right_mask, sub_mask, True)
                    )
                elif event == NodeEvent.DUPLICATION:
                    total_cost += costs.get(
                        EdgeEvent.SEGMENTAL_LOSS,
                        min(
                            subseq_segment_dist(left_mask, sub_mask, True)
                            + subseq_segment_dist(right_mask, sub_mask, False),
                            subseq_segment_dist(left_mask, sub_mask, False)
                            + subseq_segment_dist(right_mask, sub_mask, True),
                        )
                    )
                else:
                    assert event == NodeEvent.HORIZONTAL_TRANSFER
                    keep_left = self.input.species_lca.is_comparable(
                        rec[node], rec[left_node]
                    )
                    total_cost += costs.get(
                        EdgeEvent.SEGMENTAL_LOSS,
                        subseq_segment_dist(left_mask, sub_mask, keep_left)
                        + subseq_segment_dist(
                            right_mask, sub_mask, not keep_left
                        )
                    )

        return total_cost

    def _unordered_labeling_cost(self):
        """Compute the unordered segmental loss cost of the labeling."""
        tree = self.input.object_tree
        costs = self.input.costs
        rec = self.object_species

        total_cost = 0

        for node in tree.traverse("preorder"):
            if not node.is_leaf():
                event = self.node_event(node)
                left_node, right_node = node.children

                node_set = set(self.syntenies[node])
                left_count = (
                    0
                    if node_set <= set(self.syntenies[left_node])
                    else 1
                )
                right_count = (
                    0
                    if node_set <= set(self.syntenies[right_node])
                    else 1
                )

                if event == NodeEvent.SPECIATION:
                    total_cost += costs.get(
                        EdgeEvent.SEGMENTAL_LOSS,
                        left_count + right_count
                    )
                elif event == NodeEvent.DUPLICATION:
                    total_cost += costs.get(
                        EdgeEvent.SEGMENTAL_LOSS,
                        min(left_count, right_count)
                    )
                else:
                    assert event == NodeEvent.HORIZONTAL_TRANSFER
                    left_conserved = self.input.species_lca.is_comparable(
                        rec[node], rec[left_node]
                    )
                    total_cost += costs.get(
                        EdgeEvent.SEGMENTAL_LOSS,
                        left_count if left_conserved else right_count
                    )

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
