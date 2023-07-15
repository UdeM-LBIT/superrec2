"""Reconciliation and evolution model."""
from typing import Self
from abc import ABC, abstractmethod
from ast import literal_eval
from collections.abc import Mapping
from dataclasses import dataclass, field, fields, asdict, replace
from sowing.comb.binary import is_binary
from sowing.node import Node
from sowing.zipper import Zipper
from sowing.indexed import index_trees, IndexedTree
from sowing import traversal
from sowing.repr import newick
from sowing.util.dataclasses import repr_default


def parse_tree(kind: type, data: str):
    from_mapping = getattr(kind, "from_mapping")
    return traversal.map(
        lambda node, edge, *_: (from_mapping(node or {}), None),
        traversal.depth(newick.parse(data)),
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
        if not is_binary(self.host_tree):
            raise InvalidReconciliation("host tree must be binary")

        if not is_binary(self.associate_tree):
            raise InvalidReconciliation("associate tree must be binary")

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


class InvalidEvent(Exception):
    def __init__(self, message, event=None):
        if event is not None:
            message += f" (at event {event!r})"

        super().__init__(message)
        self.event = event


@dataclass(frozen=True, slots=True, repr=False)
class Event(Associate, ABC):
    """Event in a cophylogeny history."""

    def associate(self) -> Associate:
        """Get the associate corresponding to this event."""
        return Associate(name=self.name, host=self.host, contents=self.contents)

    def anon_associate(self) -> Associate:
        """Get the unnamed associate corresponding to this event."""
        return Associate(host=self.host, contents=self.contents)

    @property
    @abstractmethod
    def arity(self) -> int:
        """Number of expected children for this event."""

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
        arity = ("a leaf", "unary", "binary")[self.arity]

        if len(children) != self.arity:
            kind = self.__class__.__name__.lower()
            raise InvalidEvent(
                f"{kind} event must be {arity}, found"
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
                raise ValueError("unknown event kind {kind!r}")


@dataclass(frozen=True, slots=True, repr=False)
class Extant(Event):
    """Terminal event in an history."""

    @property
    def arity(self) -> int:
        return 0

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
        return Extant(**asdict(Associate.from_mapping(data)))


@dataclass(frozen=True, slots=True, repr=False)
class Codiverge(Event):
    """Event where an associate follows a divergence of its host."""

    @property
    def arity(self) -> int:
        return 2

    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        Event.validate(self, host_index, children)
        host_node = host_index[self.host]

        left_host = host_node.down(0).node.data.name
        right_host = host_node.down(1).node.data.name

        assoc = self.anon_associate()
        expected_children = (assoc.switch(left_host), assoc.switch(right_host))

        if set(children) != set(expected_children):
            raise InvalidEvent(
                f"codivergence event children are {children},"
                f" expected {expected_children}",
                self,
            )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        return Codiverge(**asdict(Associate.from_mapping(data)))


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

    @property
    def arity(self) -> int:
        target, remainder = self.split(self.segment)
        return 1 if self.cut and not remainder.contents else 2

    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        Event.validate(self, host_index, children)

        if self.result not in range(self.arity):
            raise InvalidEvent(
                f"divergence result index {self.result} is out of"
                f" bounds (should be in range [0..{self.arity - 1}])",
                self,
            )

        result = children[self.result]
        assoc = self.anon_associate()

        if self.transfer and host_index.is_comparable(result.host, assoc.host):
            raise InvalidEvent(
                f"transfer-divergence target host {result.host!r}"
                f" is comparable to its origin host {assoc.host!r}",
                self,
            )

        if not self.transfer and result.host != assoc.host:
            raise InvalidEvent(
                f"copy-divergence result host {result.host!r}"
                f" differs from its parent host {assoc.host!r}",
                self,
            )

        segment, remainder = assoc.split(self.segment)
        expect_result = segment.switch(result.host)

        if result != expect_result:
            raise InvalidEvent(
                f"divergence result child is {result}, expected {expect_result}"
            )

        if self.arity == 2:
            conserved = children[1 - self.result]
            expect_conserved = remainder if self.cut else assoc

            if conserved != expect_conserved:
                raise InvalidEvent(
                    f"divergence conserved child is {result},"
                    f" expected {expect_conserved}"
                )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        data = dict(data)
        attrs = {}

        if "segment" in data:
            attrs["segment"] = _contents_from_str(data.pop("segment"))

        if "cut" in data:
            attrs["cut"] = _bool_from_str(data.pop("cut"))

        if "transfer" in data:
            attrs["transfer"] = _bool_from_str(data.pop("transfer"))

        if "result" in data:
            attrs["result"] = int(data.pop("result"))

        associate = asdict(Associate.from_mapping(data))
        return Diverge(**associate, **attrs)


@dataclass(frozen=True, slots=True, repr=False)
class Gain(Event):
    """Event where an associate gains new contents."""

    # Added contents
    gained: GainedContents = ()

    @property
    def arity(self) -> int:
        return 1

    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        Event.validate(self, host_index, children)
        child = children[0]
        expect_child = self.anon_associate().gain(self.gained)

        if child != expect_child:
            raise InvalidEvent(
                f"gain event child is {child}, expected {expect_child}",
                self,
            )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        data = dict(data)
        attrs = {}

        if "gained" in data:
            attrs["gained"] = _contents_from_str(data.pop("gained"))

        associate = asdict(Associate.from_mapping(data))
        return Gain(**associate, **attrs)


@dataclass(frozen=True, slots=True, repr=False)
class Loss(Event):
    """Event where an associate loses part or all of its contents."""

    # Lost segment
    segment: Segment = ()

    @property
    def arity(self) -> int:
        lost, remainder = self.split(self.segment)
        return 1 if remainder.contents else 0

    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        Event.validate(self, host_index, children)

        if self.arity == 1:
            _, expect_child = self.anon_associate().split(self.segment)
            child = children[0]

            if child != expect_child:
                raise InvalidEvent(
                    f"loss event child is {child}, expected {expect_child}",
                    self,
                )

    @staticmethod
    def from_mapping(data: Mapping) -> Self:
        data = dict(data)
        attrs = {}

        if "segment" in data:
            attrs["segment"] = _contents_from_str(data.pop("segment"))

        associate = asdict(Associate.from_mapping(data))
        return Loss(**associate, **attrs)


@dataclass(frozen=True, slots=True)
@index_trees
class History:
    """Evolutionary history of an associate phylogeny inside an host phylogeny."""

    # Host phylogeny
    host_tree: Node[Host, None]
    host_index = field(metadata={"index_from_tree": "host_tree"})

    # History tree, with event and associate information at each node
    event_tree: Node[Event, None]

    def compress(self) -> Reconciliation:
        """
        Reduce this history to a binary associate phylogeny mapped onto
        its host phylogeny.
        """

        def compress_event(cursor: Zipper[Event, None]) -> Zipper[Associate, None]:
            node = cursor.node
            event = node.data

            sampled = self.host_index[event.host].node.data.sampled
            new_node = node.replace(data=event.associate())

            match len(node.edges):
                case 0:
                    if sampled and isinstance(event, Extant):
                        return cursor.replace(node=new_node)

                    return cursor.replace(node=None)

                case 1:
                    child = node.edges[0].node
                    return cursor.replace(node=child)

                case 2:
                    return cursor.replace(node=new_node)

        return Reconciliation(
            host_tree=self.host_tree,
            associate_tree=traversal.fold(
                compress_event,
                traversal.depth(self.event_tree, preorder=False),
            ),
        )

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
