"""Generic utilities for implementing dynamic programming algorithms."""
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum, auto
from itertools import product
from typing import (
    Any,
    Callable,
    Generic,
    Optional,
    overload,
    Protocol,
    List,
    Iterable,
    Set,
    Tuple,
    TypeVar,
    Union,
)
from infinity import Infinity, inf, is_infinite


@dataclass
class Dimension:
    pass


@dataclass
class ListDimension(Dimension):
    """
    List dimension.

    Declares a dimension of a dynamic programming table that is a list of
    given length.
    """

    length: int


@dataclass
class DictDimension(Dimension):
    """
    Dictionary dimension.

    Declares a dimension of a dynamic programming table that is a dictionary.
    It is initially empty and accepts any hashable object as a key.
    """

    pass


def _generate_table(dimensions: Iterable[Dimension]):
    if not dimensions:
        return None

    dim = dimensions[0]
    rem = dimensions[1:]

    if isinstance(dim, ListDimension):
        return [_generate_table(rem) for i in range(dim.length)]
    elif isinstance(dim, DictDimension):
        return defaultdict(lambda: _generate_table(rem))
    else:
        raise RuntimeError(f"Unknown dimension type: {dim}")


InfoType = TypeVar("InfoType")
ValueType = TypeVar("ValueType")


@dataclass(frozen=True)
class Candidate(Generic[ValueType, InfoType]):
    """Candidate value for an entry of a dynamic programming table."""

    # Value of the candidate
    value: Union[ValueType, Infinity]

    # Arbitrary info tag attached to this value
    info: Optional[InfoType] = None


class MergePolicy(Enum):
    """Policy on choosing between values for an entry of the table."""

    # Keep the lowest values
    MIN = auto()

    # Keep the highest values
    MAX = auto()


class RetentionPolicy(Enum):
    """Policy on retaining info tags associated to values."""

    # Don’t store any info tag in the entries
    # (used when only the numeric values are of interest)
    NONE = auto()

    # When several optimal values exist for an entry, keep any info tag
    # from one of them
    ANY = auto()

    # Keep all info tags for optimal values of an entry
    ALL = auto()


class EntryProtocol(Generic[ValueType, InfoType], Protocol):
    """Protocol for entries of a dynamic programming table."""

    def is_infinite(self) -> bool:
        """Check whether the value stored in this entry is infinite."""

    def value(self) -> Union[ValueType, Infinity]:
        """Get the value currently stored in this entry."""

    def infos(self) -> Set[InfoType]:
        """Get the list of custom info tags associated to the current value."""

    def info(self) -> Optional[InfoType]:
        """Get any info tag associated to the current value."""

    def update(self, *candidates: Candidate[ValueType, InfoType]) -> None:
        """Receive a set of candidates and keep the optimal ones."""

    def combine(
        self,
        other: "EntryProtocol[ValueType, InfoType]",
        combinator: Callable[
            [Candidate[ValueType, InfoType], Candidate[ValueType, InfoType]],
            ValueType,
        ],
    ) -> "EntryProtocol[ValueType, Tuple[InfoType, InfoType]]":
        """
        Combine the candidates from two entries.

        :param other: other entry to combine
        :param combinator: called for each pair of candidates from the product
            of the two entries to produce combined candidates
        :returns: resulting combined entry
        """

    def __eq__(self, other: Any) -> bool:
        """Test equality of value and info tags to another entry."""

    def __iter__(self):
        """Get the list of candidates of this entry."""

    def __len__(self):
        """Get the number of info tags in this entry."""


class Entry(Generic[ValueType, InfoType]):
    """Type for entries stored in dynamic programming tables."""

    @overload
    def __init__(
        self,
        merge_policy: MergePolicy,
        retention_policy: RetentionPolicy,
    ):
        """Create a default-initialized entry matching given policies."""
        ...

    def __init__(
        self,
        value: Union[ValueType, Infinity],
        infos: Iterable[InfoType],
        merge_policy: MergePolicy = None,
        retention_policy: RetentionPolicy = None,
    ):
        """Create an entry with initial value and info tags."""
        if (
            merge_policy is None
            and retention_policy is None
            and type(value) is MergePolicy
            and type(infos) is RetentionPolicy
        ):
            self._value = inf if value == MergePolicy.MIN else -inf
            self._infos = set()
            self._merge_policy = value
            self._retention_policy = infos
        else:
            self._value = value
            self._infos = set(infos)
            self._merge_policy = (
                merge_policy if merge_policy is not None else MergePolicy.MIN
            )
            self._retention_policy = (
                retention_policy
                if retention_policy is not None
                else RetentionPolicy.NONE
            )

    def is_infinite(self) -> bool:
        return is_infinite(self._value)

    def value(self) -> Union[ValueType, Infinity]:
        return self._value

    def infos(self) -> Set[InfoType]:
        return self._infos

    def info(self) -> Optional[InfoType]:
        if self._infos:
            return min(self._infos)
        else:
            return None

    def update(self, *candidates: Candidate[ValueType, InfoType]) -> None:
        is_min = self._merge_policy == MergePolicy.MIN
        is_max = self._merge_policy == MergePolicy.MAX

        is_any = self._retention_policy == RetentionPolicy.ANY
        is_all = self._retention_policy == RetentionPolicy.ALL

        for candidate in candidates:
            value = candidate.value
            info = candidate.info

            if self._value == value:
                if info and (is_all or (is_any and not self._infos)):
                    self._infos.add(info)

                self._value = value

            if (is_min and self._value > value) or (
                is_max and self._value < value
            ):
                if info and (is_all or is_any):
                    self._infos = {info}

                self._value = value

    def combine(
        self,
        other: "EntryProtocol[ValueType, InfoType]",
        combinator: Callable[
            [Candidate[ValueType, InfoType], Candidate[ValueType, InfoType]],
            ValueType,
        ],
    ) -> "EntryProtocol[ValueType, Tuple[InfoType, InfoType]]":
        result = Entry(self._merge_policy, self._retention_policy)

        for ours, theirs in product(self._infos, other.infos()):
            result.update(
                combinator(
                    Candidate(self._value, ours),
                    Candidate(other.value(), theirs),
                )
            )

        return result

    def __eq__(self, other: Any):
        if not callable(getattr(other, "value", None)):
            return False

        if not callable(getattr(other, "infos", None)):
            return False

        return self.value() == other.value() and self.infos() == other.infos()

    def __iter__(self):
        for info in self._infos:
            yield Candidate(self._value, info)

    def __len__(self):
        return len(self._infos)


class Table(Generic[ValueType, InfoType]):
    """
    Generic dynamic programming table.

    This class manages a multidimensional table of numeric values, makes sure
    that each entry stores the lowest (or highest) value that it receives,
    and retains custom “info tags” attached to each optimal value.

    Info tags can be used as part of a dynamic programming algorithm to
    reconstruct an optimal solution.
    """

    def __init__(
        self,
        dimensions: Iterable[Dimension],
        merge_policy: MergePolicy = MergePolicy.MIN,
        retention_policy: RetentionPolicy = RetentionPolicy.NONE,
    ):
        """
        Initialize the table.

        :param dimensions: list of dimensions of the table, specified either as
            :class:`ListDimension`s or :class:`DictDimension`s
        :param merge_policy: select between keeping the highest or lowest
            values in each entry
        :param retention_policy: select between keeping no info tags, the info
            tag of any optimal value, or the info tags of all optimal values
        """
        self.merge_policy = merge_policy
        self.retention_policy = retention_policy
        self.dimensions = dimensions
        self._table = _generate_table(dimensions)

    def entry(
        self,
        value: Union[None, ValueType, Infinity] = None,
        infos: Union[None, Iterable[InfoType]] = None,
    ) -> Entry[ValueType, InfoType]:
        """
        Create an empty entry with policies matching those of this table.

        This is a shortcut useful for algorithms where parts of an entry are
        computed separately and then combined using :meth:`Entry.combine`.
        """
        if value is None and infos is None:
            return Entry(self.merge_policy, self.retention_policy)
        else:
            return Entry(
                value,
                infos,
                self.merge_policy,
                self.retention_policy,
            )

    def keys(self) -> Iterable[Any]:
        """Get the set of keys defined in the first dimension."""
        return TableProxy(self, ()).keys()

    def __iter__(self) -> Iterable[Any]:
        return self.keys()

    def __getitem__(
        self, key: Any
    ) -> Union[
        "EntryProtocol[ValueType, InfoType]",
        "TableProxy[ValueType, InfoType]",
    ]:
        """
        Retrieve a subset or an entry of the table.

        :param key: tuple specifying a key for the first dimension of the
            table. If the table has more than a dimension, further indexing
            operations must be chained to get to an actual entry. Incomplete
            addressing will return a proxy object that acts as a view on a
            part of the table
        :returns: table entry if the address is complete,
            or table proxy object for partial addresses
        """
        return TableProxy(self, ())[key]

    def __setitem__(
        self, key: Any, candidate: Candidate[ValueType, InfoType]
    ) -> None:
        """
        Assign a candidate to an entry of the table.

        Only a candidate that matches the current MergePolicy will be stored,
        i.e., only if it has the lowest or highest value.
        """
        TableProxy(self, ())[key] = candidate


class EntryProxy(Generic[ValueType, InfoType]):
    """
    Virtual entry returned for entries that may not yet exist.

    If the underlying entry does not exist, this proxy entry behaves like an
    infinite-value entry with no info tags. Trying to assign a non-infinity
    candidate to this entry will instantiate it in its parent dynamic
    programming table.

    If the underlying entry does exist, this instance will act as a transparent
    proxy to the real entry.
    """

    def __init__(self, parent: Table[ValueType, InfoType], key: Any):
        self._parent = parent
        self._key = key

    def _get_real(self) -> Optional[Entry[ValueType, InfoType]]:
        entry = self._parent._table

        for x in self._key:
            entry = entry[x]

        return entry

    def is_infinite(self) -> Optional[InfoType]:
        real = self._get_real()

        if real is None:
            return True
        else:
            return real.is_infinite()

    def value(self) -> Union[ValueType, Infinity]:
        real = self._get_real()

        if real is None:
            return inf if self._parent.merge_policy == MergePolicy.MIN else -inf
        else:
            return real.value()

    def infos(self) -> Set[InfoType]:
        real = self._get_real()

        if real is None:
            return set()
        else:
            return real.infos()

    def info(self) -> Optional[InfoType]:
        real = self._get_real()

        if real is None:
            return None
        else:
            return real.info()

    def update(self, *candidates: Candidate[ValueType, InfoType]) -> None:
        if any(not (is_infinite(candidate.value)) for candidate in candidates):
            entry = self._parent._table

            for x in self._key[:-1]:
                entry = entry[x]

            if entry[self._key[-1]] is None:
                entry[self._key[-1]] = self._parent.entry()

            entry[self._key[-1]].update(*candidates)

    def combine(
        self,
        other: "EntryProtocol[ValueType, InfoType]",
        combinator: Callable[
            [Candidate[ValueType, InfoType], Candidate[ValueType, InfoType]],
            ValueType,
        ],
    ) -> "EntryProtocol[ValueType, Tuple[InfoType, InfoType]]":
        real = self._get_real()

        if real is None:
            return self
        else:
            return real.combine(other, combinator)

    def __eq__(self, other: Any):
        if not callable(getattr(other, "value", None)):
            return False

        if not callable(getattr(other, "infos", None)):
            return False

        return self.value() == other.value() and self.infos() == other.infos()

    def __iter__(self):
        real = self._get_real()

        if real is not None:
            yield from real.__iter__()

    def __len__(self):
        real = self._get_real()

        if real is None:
            return 0
        else:
            return real.__len__()


class TableProxy(Generic[ValueType, InfoType]):
    """
    Proxy object for dynamic programming tables.

    An instance of this class is effectively a window on a subset of a dynamic
    programming table, where the values for dimensions 1 to i are fixed and
    dimensions i + 1 to n remain to be set. Such instances are used when
    indexing the table.
    """

    def __init__(self, parent: Table[ValueType, InfoType], prefix: Any):
        self.parent = parent
        self.prefix = prefix

    def keys(self) -> Iterable[Any]:
        """Get the set of keys defined in the next dimension."""
        entry = self.parent._table

        for x in self.prefix:
            entry = entry[x]

        if type(entry) == list:
            yield from range(len(entry))
        else:
            yield from entry.keys()

    def __iter__(self) -> Iterable[Any]:
        return self.keys()

    def __getitem__(
        self, key: Any
    ) -> Union[
        None,
        EntryProtocol[ValueType, InfoType],
        "TableProxy[ValueType, InfoType]",
    ]:
        """
        Set the value of the next dimension, or retrieve an entry.

        If more than one dimensions are unset in this proxy, this will
        return a new proxy with one more dimension set to the given key.
        Otherwise, this will return the entry corresponding to this
        proxy’s keys plus the given key.
        """
        if len(self.prefix) + 1 == len(self.parent.dimensions):
            entry = self.parent._table

            for x in self.prefix:
                entry = entry[x]

            return EntryProxy(self.parent, self.prefix + (key,))
        else:
            return TableProxy(self.parent, self.prefix + (key,))

    def __setitem__(
        self, key: Any, candidate: Candidate[ValueType, InfoType]
    ) -> None:
        """
        Update the entry at the given key.

        This only works if all dimensions but the last one are set in this
        proxy. In this case, the :meth:`Entry.update` method of the entry
        corresponding to the given key will be called with the given
        candidates. If no entry corresponds to the given key yet, one is
        created on the fly.
        """
        if len(self.prefix) + 1 == len(self.parent.dimensions):
            EntryProxy(self.parent, self.prefix + (key,)).update(candidate)
        else:
            raise TypeError(
                f"Cannot address a dynamic programming table \
with {len(self.parent.dimensions)} dimension(s) using only \
{len(self.prefix) + 1} value(s)"
            )
