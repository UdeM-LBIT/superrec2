"""Reconciliation and evolution model."""
from typing import Self, overload
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, replace
from sowing.comb.binary import is_binary
from sowing.node import Node
from sowing.zipper import Zipper
from sowing.indexed import index_trees
from sowphy.clade import Clade
from sowing import traversal


# Ordered or unordered contents of an associate
Contents = tuple[str, ...] | frozenset[str]


# Segment or subset of the existing contents of an associate
Segment = tuple[int, int] | frozenset[str]


# Subsegment (with position) or subset of contents gained in an associate
GainedContents = tuple[int, tuple[str, ...]] | frozenset[str]


@overload
def insert_gain(
    contents: tuple[str, ...],
    gain: tuple[int, tuple[str, ...]],
) -> tuple[str, ...]:
    ...


@overload
def insert_gain(
    contents: frozenset[str],
    gain: frozenset[str],
) -> frozenset[str]:
    ...


def insert_gain(contents: Contents, gain: GainedContents) -> Contents:
    """
    Add an ordered or unordered gain inside existing string contents.

    :param contents: original contents
    :param gain: gained segment, with insertion index if ordered
    :returns: resulting contents
    """
    if isinstance(contents, tuple) and isinstance(gain, tuple):
        position, segment = gain
        return contents[:position] + segment + contents[position:]

    if isinstance(contents, frozenset) and isinstance(gain, frozenset):
        return contents | gain

    raise TypeError


@overload
def extract_segment(
    contents: tuple[str, ...],
    segment: tuple[int, int],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    ...


@overload
def extract_segment(
    contents: frozenset[str],
    segment: frozenset[str],
) -> frozenset[str]:
    ...


def extract_segment(contents: Contents, segment: Segment) -> Contents:
    """
    Extract an ordered or unordered segment of strings from another.

    :param contents: original string segment
    :param gain: lost segment
    :returns: extracted segment and resulting segment
    """
    if isinstance(contents, tuple) and isinstance(segment, tuple):
        start, end = segment
        return contents[start:end], contents[:start] + contents[end:]

    if isinstance(contents, frozenset) and isinstance(segment, frozenset):
        return segment, contents - segment

    raise TypeError


def show_contents(contents: Contents) -> str:
    if isinstance(contents, tuple):
        return f"({', '.join(repr(item) for item in contents)})"
    else:
        return f"{{{', '.join(repr(item) for item in sorted(contents))}}}"


@dataclass(frozen=True, slots=True)
class Associate:
    """Phylogenetic entity attached to a host."""

    # Host of the associate
    host: Clade

    # Associate metadata, if any
    clade: Clade = Clade()

    # Associate contents, if applicable
    contents: Contents = ()


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
    host_tree: Node[Clade, None]
    host_index = field(metadata={"index_from_tree": "host_tree"})

    # Associate phylogeny, with host and contents information (if known)
    associate_tree: Node[Clade | Associate, None]


    def is_complete(self) -> bool:
        """Check that all ancestral and terminal nodes have associate information."""
        return all(
            isinstance(cursor.node.data, Associate)
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
            self._validate_at(cursor.node)

    def _validate_at(self, node: Node[Clade | Associate, None]) -> None:
        associate = node.data
        is_leaf = len(node.edges) == 0

        if is_leaf:
            if not isinstance(associate, Associate):
                kind = associate.__class__.__name__
                raise InvalidReconciliation(
                    "leaf associates must be labeled by Associate instances,"
                    f" found Associate of type {kind}",
                    node,
                )

        if isinstance(associate, Associate):
            if associate.host not in self.host_index:
                raise InvalidReconciliation(
                    f"associate host {associate.host!r} does not exist in host tree",
                    node,
                )

            if is_leaf:
                if len(self.host_index[associate.host].edges) != 0:
                    raise InvalidReconciliation(
                        f"leaf associate host {associate.host!r} is not terminal",
                        node,
                    )

    def erase(self) -> Self:
        """
        Reduce this reconciliation to its original input by removing all
        ancestral associate mappings.
        """

        def erase_associate(
            cursor: Zipper[Clade | Associate, None]
        ) -> Zipper[Clade | Associate, None]:
            node = cursor.node
            associate = node.data

            if len(node.edges) == 0:
                return cursor

            node = node.replace(data=associate.clade)
            return cursor.replace(node=node)

        return replace(
            self,
            associate_tree=traversal.fold(
                erase_associate,
                traversal.depth(self.associate_tree),
            ),
        )


@dataclass(frozen=True, slots=True)
class Event(ABC):
    """Event in a cophylogeny history."""

    # Event associate and host metadata
    associate: Associate

    @property
    @abstractmethod
    def arity(self) -> int:
        """Number of expected children for this event."""
        ...


@dataclass(frozen=True, slots=True)
class Extant(Event):
    """Terminal event in an history."""

    @property
    def arity(self) -> int:
        return 0


@dataclass(frozen=True, slots=True)
class Codiverge(Event):
    """Event where an associate follows a divergence of its host."""

    @property
    def arity(self) -> int:
        return 2


@dataclass(frozen=True, slots=True)
class Diverge(Event):
    """Event where an associate diverges inside its host."""

    # Index of the child node that results from the divergence
    result: int

    # Segment of the associate contents targeted by the event
    segment: Segment = ()

    # Whether the target segment is cut out (True) or duplicated (False)
    cut: bool = False

    @property
    def arity(self) -> int:
        return 2


@dataclass(frozen=True, slots=True)
class Transfer(Event):
    """Event where an associate is transferred back into the host tree."""

    @property
    def arity(self) -> int:
        return 1


@dataclass(frozen=True, slots=True)
class Gain(Event):
    """Event where an associate gains new contents."""

    # Added contents
    gain: GainedContents

    @property
    def arity(self) -> int:
        return 1


@dataclass(frozen=True, slots=True)
class Loss(Event):
    """Event where an associate loses part or all of its contents."""

    # Lost segment
    segment: Segment = ()

    @property
    def arity(self) -> int:
        lost, remainder = extract_segment(self.associate.contents, self.segment)
        return 1 if remainder else 0


@dataclass(frozen=True, slots=True)
@index_trees
class History:
    """Evolutionary history of an associate phylogeny inside an host phylogeny."""

    # Host phylogeny
    host_tree: Node[Clade, None]
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
            associate = event.associate
            new_node = node.replace(data=associate)

            match event.arity:
                case 0:
                    is_sampled = associate.host.props.get("sampled", True)
                    is_loss = isinstance(event, Loss)

                    if is_sampled and not is_loss:
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
        for cursor in traversal.depth(self.event_tree):
            self._validate_at(cursor.node)

    def _validate_at(self, node: Node[Associate, None]) -> bool:
        event = node.data
        associate = event.associate
        contents = associate.contents

        # Check associate host
        try:
            host = self.host_index[associate.host]
        except KeyError:
            raise InvalidReconciliation(
                f"undefined associate host {associate.host!r}",
                node,
            )

        # Check out-degree depending on event kind
        kind = event.__class__.__name__
        arity = ("a leaf", "unary", "binary")[event.arity]

        if len(node.edges) != event.arity:
            raise InvalidReconciliation(
                f"{kind} node must be {arity}, found"
                f" {len(node.edges)} child(ren) instead",
                node,
            )

        # Check event-specific constraints
        match event:
            case Extant():
                if len(host.edges) != 0:
                    raise InvalidReconciliation(
                        f"leaf associate host {host.data} is not terminal",
                        node,
                    )

            case Codiverge():
                left = node.edges[0].node.data.associate
                right = node.edges[1].node.data.associate

                children_hosts = (left.host, right.host)
                host_children = (host.edges[0].node.data, host.edges[1].node.data)

                if set(children_hosts) != set(host_children):
                    raise InvalidReconciliation(
                        f"codivergence children are linked to {children_hosts}"
                        f" instead of their host’s children {host_children}",
                        node,
                    )

                unequal = None
                child_contents = None

                if right.contents != contents:
                    unequal = "right"
                    child_contents = right.contents

                if left.contents != contents:
                    unequal = "left"
                    child_contents = left.contents

                if unequal is not None:
                    raise InvalidReconciliation(
                        f"codivergence {unequal} child contents"
                        f" {show_contents(child_contents)} do not equal"
                        f" its parent’s contents {show_contents(contents)}",
                        node,
                    )

            case Diverge(cut=cut, result=result_idx, segment=segment):
                if result_idx not in (0, 1):
                    raise InvalidReconciliation(
                        f"result index is {result_idx}, expected 0 or 1",
                        node,
                    )

                result = node.edges[result_idx].node.data.associate
                conserved = node.edges[1 - result_idx].node.data.associate
                target, remainder = extract_segment(contents, segment)

                unequal = None
                child_host = None

                if result.host != associate.host:
                    unequal = "result"
                    child_host = result.host

                if conserved.host != associate.host:
                    unequal = "conserved"
                    child_host = conserved.host

                if unequal is not None:
                    raise InvalidReconciliation(
                        f"divergence {unequal} child host {child_host}"
                        f" differs from its parent host {associate.host}",
                        node,
                    )

                if result.contents != target:
                    raise InvalidReconciliation(
                        f"divergence result contents {show_contents(result.contents)}"
                        f" differ from the targeted segment {show_contents(target)}",
                        node,
                    )

                if cut:
                    if conserved.contents != remainder:
                        raise InvalidReconciliation(
                            "cut-divergence conserved child contents"
                            f" {show_contents(conserved.contents)} differ from the"
                            f" remaining contents {show_contents(remainder)}",
                            node,
                        )
                else:
                    if conserved.contents != contents:
                        raise InvalidReconciliation(
                            "copy-divergence conserved child contents"
                            f" {show_contents(conserved.contents)} differ from its"
                            f" parent’s contents {show_contents(contents)}",
                            node,
                        )

            case Transfer():
                child = node.edges[0].node.data.associate
                child_host = self.host_index[child.host]

                if self.host_index.is_comparable(child.host, host):
                    raise InvalidReconciliation(
                        f"transfer child host {child.host}"
                        f" is comparable to its origin host {associate.host}",
                        node,
                    )

                if child.contents != contents:
                    raise InvalidReconciliation(
                        f"transfer child contents {show_contents(child.contents)}"
                        " differ from its parent’s contents"
                        f" {show_contents(contents)}",
                        node,
                    )

            case Gain(gain=gain):
                child = node.edges[0].node.data.associate

                if child.host != associate.host:
                    raise InvalidReconciliation(
                        f"gain child host {child.host}"
                        f" differs from its parent host {associate.host}",
                        node,
                    )

                result = insert_gain(contents, gain)

                if child.contents != result:
                    raise InvalidReconciliation(
                        f"gain child contents {show_contents(child.contents)} differ"
                        f" from the expected gain result {show_contents(result)}",
                        node,
                    )

            case Loss(segment=segment):
                if len(node.edges) == 1:
                    child = node.edges[0].node.data.associate

                    if child.host != associate.host:
                        raise InvalidReconciliation(
                            f"loss child host {child.host}"
                            f" differs from its parent host {associate.host}",
                            node,
                        )

                    _, result = extract_segment(contents, segment)

                    if child.contents != result:
                        raise InvalidReconciliation(
                            f"loss child contents {show_contents(child.contents)}"
                            " differ from the expected loss result"
                            f" {show_contents(result)}",
                            node,
                        )

            case _:
                raise InvalidReconciliation(
                    f"unexpected {kind} node in associate tree",
                    node,
                )
