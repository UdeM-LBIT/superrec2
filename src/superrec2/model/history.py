"""Reconciliation and evolution model."""
from typing import Self
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, replace
from sowing.comb.binary import is_binary
from sowing.node import Node
from sowing.zipper import Zipper
from sowing.indexed import index_trees
from sowphy.clade import Clade
from sowing import traversal


class InvalidReconciliation(Exception):
    def __init__(self, message, node=None):
        if node is not None:
            message += f" (at node {node!r})"

        super().__init__(message)
        self.node = node


@dataclass(frozen=True, slots=True)
class Associate:
    """Phylogenetic entity attached to a host."""

    # Host of the associate
    host: Clade

    # Associate metadata, if any
    clade: Clade = Clade()

    # Associate contents, if applicable
    contents: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
@index_trees
class Reconciliation:
    """Mapping of an associate phylogeny onto an host phylogeny."""

    # Host phylogeny
    host_tree: Node[Clade, None]
    host_index = field(metadata={"index_from_tree": "host_tree"})

    # Associate phylogeny, with host and contents information (if known)
    associate_tree: Node[Clade | Associate, None]

    # Whether the contents of the associates are to be interpreted as
    # strings (True) or sets (False)
    ordered: bool = True

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
    segment: tuple[int, int] = (0, 0)

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

    # Index at which the associate contents are inserted
    index: int

    # Added contents
    gain: tuple[str, ...]

    @property
    def arity(self) -> int:
        return 1


@dataclass(frozen=True, slots=True)
class Loss(Event):
    """Event where an associate loses part or all of its contents."""

    # Lost segment
    segment: tuple[int, int]

    @property
    def arity(self) -> int:
        if self.segment == (0, len(self.associate.contents)):
            return 0
        else:
            return 1


@dataclass(frozen=True, slots=True)
@index_trees
class History:
    """Evolutionary history of an associate phylogeny inside an host phylogeny."""

    # Host phylogeny
    host_tree: Node[Clade, None]
    host_index = field(metadata={"index_from_tree": "host_tree"})

    # History tree, with event and associate information at each node
    event_tree: Node[Event, None]

    # Whether the contents of the associates are to be interpreted as
    # strings (True) or sets (False)
    ordered: bool = True

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
            ordered=self.ordered,
        )

    def validate(self) -> None:
        """
        Check that this history is valid.

        A history is valid if each node is labeled with an event and
        if this labeling follows the constraints for that event.

        :raises InvalidReconciliation: if any node is invalid
        """
        # self.compress().validate()

        for cursor in traversal.depth(self.event_tree):
            self._validate_at(cursor.node)

    def _contents_equal(self, first: tuple[str, ...], second: tuple[str, ...]) -> bool:
        if self.ordered:
            return first == second
        else:
            return set(first) == set(second)

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

                if not self._contents_equal(right.contents, contents):
                    unequal = "right"
                    child_contents = right.contents

                if not self._contents_equal(left.contents, contents):
                    unequal = "left"
                    child_contents = left.contents

                if unequal is not None:
                    raise InvalidReconciliation(
                        f"codivergence {unequal} child contents {child_contents}"
                        f" does not equal its parent’s contents {contents}",
                        node,
                    )

            case Diverge(cut=cut, result=result_idx, segment=(start, end)):
                if result_idx not in (0, 1):
                    raise InvalidReconciliation(
                        f"result index is {result_idx}, expected 0 or 1",
                        node,
                    )

                result = node.edges[result_idx].node.data.associate
                conserved = node.edges[1 - result_idx].node.data.associate
                target = contents[start:end]
                remainder = contents[:start] + contents[end:]

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

                if not self._contents_equal(result.contents, target):
                    raise InvalidReconciliation(
                        f"divergence result contents {result.contents}"
                        f" differ from the targeted segment {target}",
                        node,
                    )

                if cut:
                    if not self._contents_equal(conserved.contents, remainder):
                        raise InvalidReconciliation(
                            "cut-divergence conserved child contents"
                            f" {conserved.contents} differ from the remaining"
                            f" contents {remainder}",
                            node,
                        )
                else:
                    if not self._contents_equal(conserved.contents, contents):
                        raise InvalidReconciliation(
                            "copy-divergence conserved child contents"
                            f" {conserved.contents} differ from its parent’s"
                            f" contents {contents}",
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

                if not self._contents_equal(child.contents, contents):
                    raise InvalidReconciliation(
                        f"transfer child contents {child.contents}"
                        f" differs from its parent’s contents {contents}",
                        node,
                    )

            case Gain(index=index, gain=gain):
                child = node.edges[0].node.data.associate

                if child.host != associate.host:
                    raise InvalidReconciliation(
                        f"gain child host {child.host}"
                        f" differs from its parent host {associate.host}",
                        node,
                    )

                result = contents[:index] + gain + contents[index:]

                if not self._contents_equal(child.contents, result):
                    raise InvalidReconciliation(
                        f"gain child contents {child.contents}"
                        f" differ from the expected gain result {result}",
                        node,
                    )

            case Loss(segment=(start, end)):
                if len(node.edges) == 1:
                    child = node.edges[0].node.data.associate

                    if child.host != associate.host:
                        raise InvalidReconciliation(
                            f"loss child host {child.host}"
                            f" differs from its parent host {associate.host}",
                            node,
                        )

                    result = contents[:start] + contents[end:]

                    if not self._contents_equal(child.contents, result):
                        raise InvalidReconciliation(
                            f"loss child contents {child.contents}"
                            f" differ from the expected loss result {result}",
                            node,
                        )

            case _:
                raise InvalidReconciliation(
                    f"unexpected {kind} node in associate tree",
                    node,
                )
