"""Reconciliation and evolution model."""

from typing import Iterable, Self, TypeVar
from abc import ABC, abstractmethod
from ast import literal_eval
from itertools import product
from collections.abc import Mapping
from dataclasses import dataclass, field, fields, asdict, replace
from immutables import Map
from ..utils.graph import Edge, shortest_paths, CycleError
from sowing.comb.binary import binarize
from sowing.node import Node
from sowing.zipper import Zipper
from sowing.indexed import index_trees, IndexedTree
from sowing import traversal
from sowing.repr import newick
from sowing.util.dataclasses import repr_default


T = TypeVar("T")


def parse_tree(datatype: type[T], data: str) -> Node[T, None]:
    """Read a Newick-formatted tree and cast its nodes to the given datatype."""
    from_mapping = getattr(datatype, "from_mapping")
    return traversal.map(
        lambda node, edge, *_: (from_mapping(node or {}), None),
        traversal.depth(newick.parse(data)),
    )


def write_tree(tree: Node[T, None]) -> str:
    """Write a tree to a Newick-formatted string."""
    return newick.write(
        traversal.map(
            lambda node, edge, *_: (Map(node.to_mapping()), None),
            traversal.depth(tree),
        )
    )


# Ordered or unordered contents of an associate
Contents = tuple[str, ...] | frozenset[str]


# Segment or subset of the existing contents of an associate
Segment = tuple[int, int] | frozenset[str]


# Subsegment (with position) or subset of contents gained in an associate
GainedContents = tuple[int, tuple[str, ...]] | frozenset[str]


def _contents_from_str(data: str) -> Contents:
    result = literal_eval(data)

    if isinstance(result, set):
        return frozenset(result)

    return result


def _contents_to_str(contents: Contents) -> str:
    if isinstance(contents, frozenset):
        if not contents:
            return "set()"

        return f"{{{str(sorted(contents))[1:-1]}}}"

    return str(contents)


def _bool_from_str(data: str) -> bool:
    if isinstance(data, (bool, int)):
        return bool(data)

    if data.lower() in ("true", "yes", "1"):
        return True

    if data.lower() in ("false", "no", "0"):
        return False

    raise ValueError


@dataclass(frozen=True, slots=True)
class Associate:
    """Phylogenetic entity attached to a host."""

    # Name of the associate
    name: str | None = None

    # Associate host (if known)
    host: str | None = None

    # Associate contents (if applicable)
    contents: Contents | None = None

    def is_complete(self) -> bool:
        """Check whether this associate has complete information."""
        return self.name is not None and self.host is not None

    def gain(self, gain: GainedContents) -> Self:
        """
        Add an ordered or unordered gain inside this associate’s contents.

        :param gain: gained segment, with insertion index if ordered
        :returns: resulting associate
        """
        if isinstance(self.contents, tuple) and isinstance(gain, tuple):
            position, segment = gain
            return Associate(
                name=self.name,
                host=self.host,
                contents=self.contents[:position] + segment + self.contents[position:],
            )

        if isinstance(self.contents, frozenset) and isinstance(gain, frozenset):
            return Associate(
                name=self.name,
                host=self.host,
                contents=self.contents | gain,
            )

        raise TypeError(f"cannot add contents of type {type(gain)}")

    def split(self, segment: Segment) -> tuple[Self, Self]:
        """
        Split this associate into two associates with disjoint contents.

        :param gain: segment to extract
        :returns: extracted and remaining associates
        """
        if isinstance(self.contents, tuple) and isinstance(segment, tuple):
            start, end = segment
            return (
                Associate(
                    name=self.name, host=self.host, contents=self.contents[start:end]
                ),
                Associate(
                    name=self.name,
                    host=self.host,
                    contents=self.contents[:start] + self.contents[end:],
                ),
            )

        if isinstance(self.contents, frozenset) and isinstance(segment, frozenset):
            if not (segment <= self.contents):
                raise ValueError(
                    f"split argument {self._repr_contents(segment)} is not a subset"
                    f" of existing contents {self._repr_contents(self.contents)}"
                )

            return (
                Associate(name=self.name, host=self.host, contents=segment),
                Associate(
                    name=self.name, host=self.host, contents=self.contents - segment
                ),
            )

        raise TypeError(f"cannot split on contents of type {type(segment)}")

    def switch(self, host: str) -> Self:
        """
        Switch this associate to a different host.

        :param host: new host name
        :returns: resulting associate
        """
        return Associate(name=self.name, host=host, contents=self.contents)

    @staticmethod
    def _repr_contents(contents: Contents) -> str:
        if isinstance(contents, frozenset):
            return f"{{{', '.join(repr(i) for i in sorted(contents))}}}"

        return repr(contents)

    def __repr__(self):
        args = []

        for item in fields(self):
            value = getattr(self, item.name)

            if value != item.default and item.repr:
                if isinstance(value, frozenset):
                    value_str = self._repr_contents(value)
                else:
                    value_str = repr(value)

                args.append(f"{item.name}={value_str}")

        classname = self.__class__.__qualname__
        return f"{classname}({', '.join(args)})"

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        data = dict(data)

        if "contents" in data:
            data["contents"] = _contents_from_str(data["contents"])

        return Associate(**data)

    def to_mapping(self) -> dict[str, str]:
        result = {}

        if self.name is not None:
            result["name"] = self.name

        if self.host is not None:
            result["host"] = self.host

        if self.contents is not None:
            result["contents"] = _contents_to_str(self.contents)

        return result


@repr_default
@dataclass(frozen=True, slots=True)
class Host:
    """Phylogenetic entity on which associates depend."""

    # Name of the host
    name: str | None = None

    # If False, represents a part of the host tree that has not been sampled
    sampled: bool = True

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        data = dict(data)

        if "sampled" in data:
            data["sampled"] = _bool_from_str(data["sampled"])

        return Host(**data)

    def to_mapping(self) -> dict[str, str]:
        result = {}

        if self.name is not None:
            result["name"] = self.name

        if not self.sampled:
            result["sampled"] = "False"

        return result


def graft_unsampled_hosts(host_tree: Node[Host, None]) -> Node[Host, None]:
    """Extend a host tree to add all possible unsampled species."""

    def graft(cursor: Zipper[Host, None]) -> Zipper[Host, None]:
        node = cursor.node
        host = node.data
        return cursor.replace(
            node=(
                Node(Host(name=f"{host.name}[P]"))
                .add(Node(Host(name=f"{host.name}[U]", sampled=False)))
                .add(node)
            ),
        )

    return traversal.fold(graft, traversal.depth(host_tree))


class InvalidReconciliation(Exception):
    def __init__(self, message, node=None):
        if node is not None:
            message += f" (at node {node!r})"

        super().__init__(message)
        self.node = node


@dataclass(frozen=True, slots=True)
@index_trees
class Reconciliation:
    """Mapping of an associate phylogeny onto an host phylogeny."""

    # Host phylogeny
    host_tree: Node[Host, None]
    host_index = field(metadata={"index_from_tree": "host_tree"})

    # Associate phylogeny, partially or completely mapped onto the host phylogeny
    associate_tree: Node[Associate, None]

    def binarize(self) -> Iterable[Self]:
        """Generate all possible binarizations of this reconciliation."""
        host_trees = binarize(self.host_tree, default=Node(Host()))
        associate_trees = binarize(self.associate_tree, default=Node(Associate()))

        for host_tree, associate_tree in product(host_trees, associate_trees):
            yield Reconciliation(host_tree=host_tree, associate_tree=associate_tree)

    def is_complete(self) -> bool:
        """Check that all ancestral and terminal nodes have complete information."""
        return all(
            cursor.node.data.is_complete()
            for cursor in traversal.depth(self.associate_tree)
        )

    def validate(self) -> None:
        """
        Check that this reconciliation is valid.

        A reconciliation is valid if all leaves have associate information,
        and if all nodes with associate information link to existing hosts.

        :raises InvalidReconciliation: if any node is invalid
        """
        for cursor in traversal.depth(self.associate_tree):
            node = cursor.node
            associate = node.data
            is_leaf = len(node.edges) == 0

            if is_leaf and not associate.is_complete():
                raise InvalidReconciliation(
                    "leaf associates must have complete information",
                    node,
                )

            if associate.host is not None:
                if associate.host not in self.host_index:
                    raise InvalidReconciliation(
                        f"associate host {associate.host!r} does not "
                        "exist in host tree",
                        node,
                    )

                if is_leaf and not self.host_index[associate.host].is_leaf():
                    raise InvalidReconciliation(
                        f"leaf associate host {associate.host!r} is not terminal",
                        node,
                    )

    def erase(self) -> Self:
        """
        Reduce this reconciliation to its original input by removing all
        ancestral associate mappings.
        """

        def erase_associate(cursor: Zipper[Associate, None]) -> Zipper[Associate, None]:
            node = cursor.node
            associate = node.data

            if len(node.edges) == 0:
                return cursor

            node = node.replace(data=Associate(name=associate.name))
            return cursor.replace(node=node)

        return replace(
            self,
            associate_tree=traversal.fold(
                erase_associate,
                traversal.depth(self.associate_tree),
            ),
        )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        return Reconciliation(
            host_tree=parse_tree(Host, data["host_tree"]),
            associate_tree=parse_tree(Associate, data["associate_tree"]),
        )

    def to_mapping(self) -> dict[str, str]:
        return {
            "host_tree": write_tree(self.host_tree),
            "associate_tree": write_tree(self.associate_tree),
        }


class InvalidEvent(Exception):
    def __init__(self, message, event=None):
        if event is not None:
            message += f" (at event {event!r})"

        super().__init__(message)
        self.event = event


@dataclass(frozen=True, slots=True, repr=False)
class Event(Associate, ABC):
    """Event in a cophylogeny history."""

    # Whether this event node will be visible as a node of the associate tree
    # after compressing the history
    apparent: bool = False

    def associate(self) -> Associate:
        """Get the associate corresponding to this event."""
        return Associate(name=self.name, host=self.host, contents=self.contents)

    def anon_associate(self) -> Associate:
        """Get the unnamed associate corresponding to this event."""
        return Associate(host=self.host, contents=self.contents)

    @abstractmethod
    def outcomes(self, host_index: IndexedTree[Host, None]) -> tuple[Associate, ...]:
        """Expected children resulting from this event."""

    @abstractmethod
    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        """
        Check that this event is valid.

        :param host_index: host tree, indexed by name
        :param children: children associates, if applicable
        :raises InvalidEvent: if the current event is invalid
        """
        # Check that associated host exists
        if self.host is None or self.host not in host_index:
            raise InvalidEvent(f"undefined event host {self.host!r}", self)

        # Check that out-degree matches expected arity
        arity = len(self.outcomes(host_index))
        arity_text = ("a leaf", "unary", "binary")[arity]

        if len(children) != arity:
            kind = self.__class__.__name__.lower()
            raise InvalidEvent(
                f"{kind} event must be {arity_text}, found"
                f" {len(children)} child(ren) instead",
                self,
            )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        data = dict(data)
        kind = data.pop("kind", "extant")

        match kind:
            case "extant":
                return Extant.from_mapping(data)

            case "codiverge":
                return Codiverge.from_mapping(data)

            case "diverge":
                return Diverge.from_mapping(data)

            case "gain":
                return Gain.from_mapping(data)

            case "loss":
                return Loss.from_mapping(data)

            case _:
                raise ValueError(f"unknown event kind {kind!r}")

    @staticmethod
    def _decode_mapping(data: Mapping) -> dict:
        data = dict(data)

        if "apparent" in data:
            data["apparent"] = _bool_from_str(data["apparent"])

        associate_data = {}

        for attr in ("name", "host", "contents"):
            if attr in data:
                associate_data[attr] = data.pop(attr)

        data.update(asdict(Associate.from_mapping(associate_data)))
        return data

    @staticmethod
    def to_mapping(self) -> dict[str, str]:
        result = super(Event, self).to_mapping()

        if self.apparent:
            result["apparent"] = "True"

        return result


@dataclass(frozen=True, slots=True, repr=False)
class Extant(Event):
    """Terminal event in an history."""

    def outcomes(self, host_index: IndexedTree[Host, None]) -> tuple[Associate, ...]:
        return ()

    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        Event.validate(self, host_index, children)

        if not host_index[self.host].is_leaf():
            raise InvalidEvent(
                f"extant host {self.host!r} is not terminal",
                self,
            )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        return Extant(**Event._decode_mapping(data))

    def to_mapping(self) -> dict[str, str]:
        result = super(Extant, self).to_mapping(self)
        result["kind"] = "extant"
        return result


@dataclass(frozen=True, slots=True, repr=False)
class Codiverge(Event):
    """Event where an associate follows a divergence of its host."""

    def outcomes(self, host_index: IndexedTree[Host, None]) -> tuple[Associate, ...]:
        host_node = host_index[self.host]
        left_host = host_node.down(0).node.data.name
        right_host = host_node.down(1).node.data.name

        assoc = self.anon_associate()
        return assoc.switch(left_host), assoc.switch(right_host)

    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        Event.validate(self, host_index, children)
        outcomes = self.outcomes(host_index)

        if children != outcomes and children != outcomes[::-1]:
            raise InvalidEvent(
                f"codivergence event children are {children}, expected {outcomes}",
                self,
            )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        return Codiverge(**Event._decode_mapping(data))

    def to_mapping(self) -> dict[str, str]:
        result = super(Codiverge, self).to_mapping(self)
        result["kind"] = "codiverge"
        return result


# Name of the placeholder host for transfer outcome child
TRANSFER_OUTCOME = "__X__"


@dataclass(frozen=True, slots=True, repr=False)
class Diverge(Event):
    """Event where an associate diverges inside its host."""

    # Segment of the associate contents targeted by the event
    segment: Segment = ()

    # Whether the target segment is cut out (True) or duplicated (False)
    cut: bool = False

    # Whether the divergence results in a separate host
    transfer: bool = False

    # Index of the child node that results from the divergence
    result: int = 0

    def outcomes(self, host_index: IndexedTree[Host, None]) -> tuple[Associate, ...]:
        assoc = self.anon_associate()
        result_host = TRANSFER_OUTCOME if self.transfer else self.host

        if self.contents is None:
            if self.cut:
                outcomes = (assoc.switch(result_host),)
            else:
                outcomes = assoc.switch(result_host), assoc
        else:
            target, remainder = assoc.split(self.segment)

            if self.cut:
                if not remainder.contents:
                    outcomes = (assoc.switch(result_host),)
                else:
                    outcomes = target.switch(result_host), remainder
            else:
                outcomes = target.switch(result_host), assoc

        if self.result == 0:
            return outcomes
        else:
            return outcomes[::-1]

    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        Event.validate(self, host_index, children)
        outcomes = self.outcomes(host_index)

        if self.result not in range(len(outcomes)):
            raise InvalidEvent(
                f"divergence result index {self.result} is out of"
                f" bounds (should be in range [0..{len(outcomes) - 1}])",
                self,
            )

        result = children[self.result]
        outcome_result = outcomes[self.result]

        if self.transfer and host_index.is_comparable(result.host, self.host):
            raise InvalidEvent(
                f"transfer target host {result.host!r}"
                f" is comparable to its origin host {self.host!r}",
                self,
            )

        if outcome_result.host == TRANSFER_OUTCOME:
            outcome_result = outcome_result.switch(result.host)

        if result != outcome_result:
            raise InvalidEvent(
                f"divergence result child is {result}, expected {outcome_result}",
                self,
            )

        if len(outcomes) > 1:
            conserved = children[1 - self.result]
            outcome_conserved = outcomes[1 - self.result]

            if conserved != outcome_conserved:
                raise InvalidEvent(
                    f"divergence conserved child is {conserved}, "
                    f"expected {outcome_conserved}",
                    self,
                )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        data = Event._decode_mapping(data)

        if "segment" in data:
            data["segment"] = _contents_from_str(data["segment"])

        if "cut" in data:
            data["cut"] = _bool_from_str(data["cut"])

        if "transfer" in data:
            data["transfer"] = _bool_from_str(data["transfer"])

        if "result" in data:
            data["result"] = int(data["result"])

        return Diverge(**data)

    def to_mapping(self) -> dict[str, str]:
        result = super(Diverge, self).to_mapping(self)
        result["kind"] = "diverge"

        if self.segment != ():
            result["segment"] = _contents_to_str(self.segment)

        if self.cut:
            result["cut"] = "True"

        if self.transfer:
            result["transfer"] = "True"

        if self.result != 0:
            result["result"] = str(self.result)

        return result


@dataclass(frozen=True, slots=True, repr=False)
class Gain(Event):
    """Event where an associate gains new contents."""

    # Added contents
    gained: GainedContents = ()

    def outcomes(self, host_index: IndexedTree[Host, None]) -> tuple[Associate, ...]:
        return (self.anon_associate().gain(self.gained),)

    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        Event.validate(self, host_index, children)
        child = children[0]
        outcome_child = self.outcomes(host_index)[0]

        if child != outcome_child:
            raise InvalidEvent(
                f"gain event child is {child}, expected {outcome_child}",
                self,
            )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        data = Event._decode_mapping(data)

        if "gained" in data:
            data["gained"] = _contents_from_str(data["gained"])

        return Gain(**data)

    def to_mapping(self) -> dict[str, str]:
        result = super(Gain, self).to_mapping(self)
        result["kind"] = "gain"

        if self.gained != ():
            result["gained"] = _contents_to_str(self.gained)

        return result


@dataclass(frozen=True, slots=True, repr=False)
class Loss(Event):
    """Event where an associate loses part or all of its contents."""

    # Lost segment
    segment: Segment = ()

    def outcomes(self, host_index: IndexedTree[Host, None]) -> tuple[Associate, ...]:
        lost, remainder = self.anon_associate().split(self.segment)

        if remainder.contents:
            return (remainder,)
        else:
            return ()

    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        Event.validate(self, host_index, children)
        outcomes = self.outcomes(host_index)

        if outcomes:
            child = children[0]
            outcome_child = outcomes[0]

            if child != outcome_child:
                raise InvalidEvent(
                    f"loss event child is {child}, expected {outcome_child}",
                    self,
                )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        data = Event._decode_mapping(data)

        if "segment" in data:
            data["segment"] = _contents_from_str(data["segment"])

        return Loss(**data)

    def to_mapping(self) -> dict[str, str]:
        result = super(Loss, self).to_mapping(self)
        result["kind"] = "loss"

        if self.segment != ():
            result["segment"] = _contents_to_str(self.segment)

        return result


@dataclass(frozen=True, slots=True)
class Epochs:
    hosts: dict[Zipper[Host, None], tuple[int, int]]
    events: dict[Zipper[Event, None], int]

    def range(self):
        start = min(before for before, _ in self.hosts.values())
        end = max(after for _, after in self.hosts.values())
        return range(start, end + 1)

    def hosts_at(
        self, start: int | None = None, end: int | None = None
    ) -> Iterable[Zipper[Host, None]]:
        return (
            host
            for host, (host_start, host_end) in self.hosts.items()
            if (start is None or start == host_start)
            and (end is None or end == host_end)
        )


class InfeasibleEpochs(Exception):
    def __init__(
        self, cycle: list[tuple[Zipper[Host, None], str] | Zipper[Event, None]]
    ):
        transfers = []
        cycle_text = []

        for item in cycle:
            if isinstance(item, tuple):
                cycle_text.append(f"{item[1]} of {item[0].node.data.name}")
            elif isinstance(item.node.data, Diverge) and item.node.data.transfer:
                source = item.node.data.host
                target = item.down(item.node.data.result).node.data.host
                transfers.append(item)
                cycle_text.append(f"transfer from {source} to {target}")

        super().__init__(
            "infeasible history because of epochs cycle: "
            + ", before ".join(cycle_text)
        )

        self.cycle = cycle
        self.transfers = transfers
        self.cycle_text = cycle_text


@dataclass(frozen=True, slots=True)
@index_trees
class History:
    """Evolutionary history of an associate phylogeny inside an host phylogeny."""

    # Host phylogeny
    host_tree: Node[Host, None]
    host_index = field(metadata={"index_from_tree": "host_tree"})

    # History tree, with event and associate information at each node
    event_tree: Node[Event, None]

    def validate(self) -> None:
        """
        Check that this history is valid.

        A history is valid if each node is labeled with an event and
        if this labeling follows the constraints for that event.

        :raises InvalidReconciliation: if any node is invalid
        """
        for cursor in traversal.depth(self.event_tree, preorder=False):
            node = cursor.node
            event = node.data

            children = tuple(edge.node.data.anon_associate() for edge in node.edges)
            event.validate(self.host_index, children)

    def epochs(self, ignore: Iterable[Zipper[Event, None]] = frozenset()) -> Epochs:
        """
        Compute minimum feasible dates for each host and event of this history,
        taking into account codivergence and horizontal transfer relations.

        :param ignore: transfer events to ignore when computing the dates
        :returns: minimal dates for each host and event
        :raises InfeasibleEpochs: if the history has no feasible datation
        """

        def host_start(host: Zipper[Host, None]) -> tuple[Zipper[Host, None], str]:
            return (host, "start")

        def host_end(host: Zipper[Host, None]) -> tuple[Zipper[Host, None], str]:
            return (host, "end")

        root = host_start(self.host_tree.unzip())
        extant_sink = object()
        nodes = [extant_sink]
        edges = []

        for host in traversal.depth(self.host_tree):
            nodes.append(host_start(host))
            nodes.append(host_end(host))

            host_data = host.node.data

            # Host intervals must not start after they end
            edges.append(Edge(start=host_start(host), end=host_end(host), weight=0))

            # Host intervals must end strictly before their descendants
            for i in range(len(host.node.edges)):
                edges.append(
                    Edge(start=host_end(host), end=host_start(host.down(i)), weight=-1)
                )

            # Sampled terminal hosts must be contemporaneous
            if host.is_leaf() and host_data.sampled:
                edges.append(Edge(start=extant_sink, end=host_end(host), weight=0))
                edges.append(Edge(start=host_end(host), end=extant_sink, weight=0))

        for event in traversal.depth(self.event_tree):
            nodes.append(event)

            event_data = event.node.data
            host = self.host_index[event_data.host]
            host_data = host.node.data

            # Sampled leaves must be contemporaneous
            if isinstance(event_data, Extant) and host_data.sampled:
                edges.append(Edge(start=extant_sink, end=event, weight=0))
                edges.append(Edge(start=event, end=extant_sink, weight=0))

            # Host intervals must enclose all their events
            edges.append(Edge(start=host_start(host), end=event, weight=0))
            edges.append(Edge(start=event, end=host_end(host), weight=0))

            if event not in ignore:
                # Parents must not come after their children
                for i in range(len(event.node.edges)):
                    edges.append(Edge(start=event, end=event.down(i), weight=0))

                # Transfers must go towards coexisting hosts
                if isinstance(event_data, Diverge) and event_data.transfer:
                    result = event.down(event_data.result)
                    target = self.host_index[result.node.data.host]
                    edges.append(Edge(start=host_start(target), end=event, weight=0))
                    edges.append(Edge(start=event, end=host_end(target), weight=0))

        # Assign minimum feasible epochs, if possible, using shortest paths
        try:
            epochs, _ = shortest_paths(root, nodes, edges)
        except CycleError as err:
            raise InfeasibleEpochs(err.args[1])

        epochs = {key: -value for key, value in epochs.items()}
        hosts_epochs = {
            host: (epochs[host_start(host)], epochs[host_end(host)])
            for host in traversal.depth(self.host_tree)
        }
        events_epochs = {
            event: epochs[event] for event in traversal.depth(self.event_tree)
        }

        return Epochs(hosts=hosts_epochs, events=events_epochs)

    def prune_unsampled(self) -> Self:
        """Remove unsampled species containing no non-extant events from the history."""

        # Collect events by the host they belong to
        events_by_host = {
            cursor.node.data.name: [] for cursor in traversal.depth(self.host_tree)
        }

        for cursor in traversal.depth(self.event_tree):
            event = cursor.node.data
            events_by_host[event.host].append(event)

        # Remove and regraft unsampled hosts that only contain extant events
        pruned_leaves = set()
        pruned_internal = set()

        def prune_hosts(cursor: Zipper[Host, None]) -> Zipper[Host, None]:
            host = cursor.node.data

            if not host.sampled and all(
                (
                    isinstance(event, Extant)
                    or (
                        isinstance(event, Loss)
                        and len(event.outcomes(self.host_index)) == 0
                    )
                )
                for event in events_by_host[host.name]
            ):
                pruned_leaves.add(host.name)
                return cursor.replace(node=None)

            if len(cursor.node.edges) == 1:
                pruned_internal.add(host.name)
                return cursor.replace(node=cursor.down().node)

            return cursor

        host_tree = traversal.fold(
            prune_hosts,
            traversal.depth(self.host_tree, preorder=False),
        )

        # Remove and regraft events happening in pruned hosts
        def prune_events(cursor: Zipper[Event, None]) -> Zipper[Event, None]:
            event = cursor.node.data
            actual_arity = len(cursor.node.edges)

            if event.host in pruned_leaves:
                return cursor.replace(node=None)

            if actual_arity < len(event.outcomes(self.host_index)):
                if actual_arity == 0:
                    return cursor.replace(node=None)

                if actual_arity == 1:
                    return cursor.replace(node=cursor.down().node)

            if event.host in pruned_internal:
                host = next(
                    edge.node.data.name
                    for edge in self.host_index[event.host].node.edges
                    if edge.node.data.name not in pruned_leaves
                )
                return cursor.replace(
                    node=cursor.node.replace(data=replace(event, host=host))
                )

            return cursor

        event_tree = traversal.fold(
            prune_events,
            traversal.depth(self.event_tree, preorder=False),
        )
        return History(host_tree=host_tree, event_tree=event_tree)

    def compress(self) -> Reconciliation:
        """
        Reduce this history to a binary associate phylogeny mapped onto
        its host phylogeny.

        :raises InvalidReconciliation: if apparent nodes are compressed
            or non-apparent nodes are not compressed
        """

        def compress_event(cursor: Zipper[Event, None]) -> Zipper[Associate, None]:
            node = cursor.node
            event = node.data

            sampled = self.host_index[event.host].node.data.sampled
            new_node = node.replace(data=event.associate())

            match len(node.edges):
                case 0:
                    if sampled and isinstance(event, Extant):
                        if not event.apparent:
                            raise InvalidEvent(
                                "sampled extant leaf should be apparent",
                                event,
                            )

                        return cursor.replace(node=new_node)

                    if event.apparent:
                        raise InvalidEvent(
                            "empty or unsampled leaf should not be apparent",
                            event,
                        )

                    return cursor.replace(node=None)

                case 1:
                    if event.apparent:
                        raise InvalidEvent("unary event should not be apparent", event)

                    child = node.edges[0].node
                    return cursor.replace(node=child)

                case 2:
                    if not event.apparent:
                        raise InvalidEvent("binary event should be apparent", event)

                    return cursor.replace(node=new_node)

        return Reconciliation(
            host_tree=self.host_tree,
            associate_tree=traversal.fold(
                compress_event,
                traversal.depth(self.event_tree, preorder=False),
            ),
        )

    def transfer_distance(self) -> int:
        """
        Compute the total transfer distance of this history.

        The distance of a single transfer event is the number of edges between
        the transfer origin and its destination in the host tree. The total
        transfer distance is the sum of the distances of all transfer events.

        :returns: the total transfer distance
        """

        def count_distance(
            cursor: Zipper[Event, None]
        ) -> Zipper[tuple[int, Event], None]:
            node = cursor.node
            event = node.data

            if isinstance(event, Diverge) and event.transfer:
                _, child = cursor.down(event.result).node.data
                delta = self.host_index.distance(event.host, child.host)
            else:
                delta = 0

            below = sum(edge.node.data[0] for edge in node.edges)
            return cursor.replace(node=Node((below + delta, event)))

        result = traversal.fold(count_distance, traversal.depth(self.event_tree))
        return result.data[0]

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        return History(
            host_tree=parse_tree(Host, data["host_tree"]),
            event_tree=parse_tree(Event, data["event_tree"]),
        )

    def to_mapping(self) -> dict[str, str]:
        return {
            "host_tree": write_tree(self.host_tree),
            "event_tree": write_tree(self.event_tree),
        }
