from typing import Any, Self, Generic, Protocol, Callable, TypeVar
from dataclasses import dataclass
from functools import wraps
from itertools import product
from immutables import Map
import operator


T = TypeVar("T")


def check_same_type(
    fun: Callable[[T, T], T]
) -> Callable[[T, Any], T | type[NotImplemented]]:
    """
    Decorate a unary method to check that its argument has the same type
    as the object it is called on, and otherwise return NotImplemented.
    """

    @wraps(fun)
    def checked_fun(self: T, other: Any) -> T:
        if isinstance(other, type(self)):
            return fun(self, other)

        return NotImplemented

    return checked_fun


class Semiring(Generic[T], Protocol):
    """Structure with distributive addition and multiplication operations."""

    value: T

    @classmethod
    def cast(cls, other: Any) -> Self:
        """Cast an element to be part of the semiring."""
        ...

    @classmethod
    def null(cls) -> Self:
        """Return the unit element for addition."""
        ...

    @classmethod
    def unit(cls) -> Self:
        """Return the unit element for multiplication."""
        ...

    @classmethod
    def make(cls, *args, **kwargs) -> Self:
        """Make an element according to the given arguments."""
        ...

    def __add__(self, other: Self) -> Self:
        ...

    def __mul__(self, other: Self) -> Self:
        ...


def make_semiring(
    typename: str,
    cast: Callable[[Any], T],
    null: T,
    unit: T,
    add: Callable[[T, T], T],
    mul: Callable[[T, T], T],
    make: Callable[..., T],
) -> Semiring[T]:
    """Helper factory to create semiring classes."""

    @dataclass(frozen=True, slots=True)
    class Result(Semiring[T]):
        value: T

        @classmethod
        def cast(cls, other: Any) -> Self:
            return cls(cast(other.value))

        @classmethod
        def null(cls) -> Self:
            return cls(null)

        @classmethod
        def unit(cls) -> Self:
            return cls(unit)

        @classmethod
        def make(cls, *args, **kwargs) -> Self:
            return cls(make(*args, **kwargs))

        @check_same_type
        def __add__(self, other: Self) -> Self:
            return self.__class__(add(self.value, other.value))

        @check_same_type
        def __mul__(self, other: Self) -> Self:
            return self.__class__(mul(self.value, other.value))

    Result.__name__ = typename
    Result.__qualname__ = typename
    return Result


def min_plus(typename: str, make: Callable[..., T] = lambda x: x) -> Semiring[T]:
    """Semiring used for finding the minimum cost of any solution."""
    return make_semiring(
        typename,
        cast=lambda x: x,
        null=float("inf"),
        unit=0,
        add=min,
        mul=operator.add,
        make=make,
    )


def max_plus(typename: str, make: Callable[..., T] = lambda x: x) -> Semiring[T]:
    """Semiring used for finding the maximum score of any solution."""
    return make_semiring(
        typename,
        cast=lambda x: x,
        null=float("-inf"),
        unit=0,
        add=max,
        mul=operator.add,
        make=make,
    )


class UnitMagma(Generic[T], Protocol):
    """Structure with an operation and a unit element."""

    value: T

    @classmethod
    def cast(cls, other: Any) -> Self:
        """Cast an element to be part of the magma."""
        ...

    @classmethod
    def unit(cls) -> Self:
        """Return the unit element for the operation."""
        ...

    @classmethod
    def make(cls, *args, **kwargs) -> Self:
        """Make an element according to the given arguments."""
        ...

    def __mul__(self, other: Self) -> Self:
        ...


def make_unit_magma(
    typename: str,
    unit: T,
    mul: Callable[[T, T], T],
    make: Callable[..., T],
) -> UnitMagma[T]:
    @dataclass(frozen=True, slots=True)
    class Result(UnitMagma[T]):
        value: T

        @classmethod
        def cast(cls, other: Any) -> Self:
            return cls(other.value)

        @classmethod
        def unit(cls) -> Self:
            return cls(unit)

        @classmethod
        def make(cls, *args, **kwargs) -> Self:
            return cls(make(*args, **kwargs))

        @check_same_type
        def __mul__(self, other: Self) -> Self:
            return self.__class__(mul(self.value, other.value))

    Result.__name__ = typename
    Result.__qualname__ = typename
    return Result


def make_unit_generator(
    typename: str, Solution: type[UnitMagma]
) -> "Semiring[UnitMagma]":
    return make_semiring(
        typename,
        cast=Solution.cast,
        null=None,
        unit=Solution.unit(),
        add=lambda sol1, sol2: sol1 if sol1 is not None else sol2,
        mul=lambda sol1, sol2: (
            sol1 * sol2 if sol1 is not None and sol2 is not None else None
        ),
        make=Solution.make,
    )


def make_generator(
    typename: str, Solution: type[UnitMagma]
) -> "Semiring[frozenset[UnitMagma]]":
    """
    Make a semiring to generate sets of solutions.

    :param typename: name of the semiring type to create
    :param Solution: solution type, with suitable composition operator
    """
    return make_semiring(
        typename,
        cast=lambda set1: frozenset(map(Solution.cast, set1)),
        null=frozenset(),
        unit=frozenset({Solution.unit()}),
        add=operator.or_,
        mul=lambda set1, set2: frozenset(
            {value1 * value2 for value1, value2 in product(set1, set2)}
        ),
        make=lambda *args, **kwargs: frozenset({Solution.make(*args, **kwargs)}),
    )


class OrderedMagma(Generic[T], UnitMagma[T], Protocol):
    """Magma with a partial order."""

    def __le__(self, other):
        ...


def tuple_ordered_magma(cls):
    """Decorate a named tuple to make it into an ordered magma."""

    @classmethod
    def unit(cls) -> Self:
        return cls()

    def __le__(self, other):
        return all(left <= right for left, right in zip(self, other))

    def __mul__(self, other):
        return self.__class__(*(left + right for left, right in zip(self, other)))

    cls.unit = unit
    cls.__le__ = __le__
    cls.__mul__ = __mul__

    return cls


def pareto(
    typename: str, Value: type[OrderedMagma]
) -> Semiring[frozenset[OrderedMagma]]:
    def select(options):
        return frozenset(
            {
                value
                for value in options
                if not any(other <= value and other != value for other in options)
            }
        )

    return make_semiring(
        typename,
        cast=lambda x: x,
        null=frozenset(),
        unit=frozenset({Value.unit()}),
        add=lambda x, y: select(x | y),
        mul=lambda x, y: select({a * b for a, b in product(x, y)}),
        make=lambda *args, **kwargs: frozenset({Value.make(*args, **kwargs)}),
    )


def viterbi(typename: str, make: Callable[..., T] = lambda x: x) -> Semiring[T]:
    """Semiring used for finding the maximum probability of any solution."""
    return make_semiring(
        typename,
        cast=lambda x: x,
        null=0,
        unit=1,
        add=max,
        mul=operator.mul,
        make=make,
    )


def boolean(typename: str, make: Callable[..., T] = lambda x: x) -> Semiring[T]:
    """Semiring used for testing the existence of any solution."""
    return make_semiring(
        typename,
        cast=lambda x: x,
        null=False,
        unit=True,
        add=operator.or_,
        mul=operator.and_,
        make=make,
    )


def count(typename: str, make: Callable[..., T] = lambda x: x) -> Semiring[T]:
    """Semiring used for counting the number of solutions."""
    return make_semiring(
        typename,
        cast=lambda x: x,
        null=0,
        unit=1,
        add=operator.add,
        mul=operator.mul,
        make=make,
    )


def make_product(
    typename: str, Left: type[Semiring], Right: type[Semiring]
) -> Semiring:
    return make_semiring(
        typename,
        cast=lambda data: (Left.cast(data[0]), Right.cast(data[1])),
        null=(Left.null(), Right.null()),
        unit=(Left.unit(), Right.unit()),
        add=lambda x, y: (x[0] + y[0], x[1] + y[1]),
        mul=lambda x, y: (x[0] * y[0], x[1] * y[1]),
        make=lambda *args, **kwargs: (
            Left.make(*args, **kwargs),
            Right.make(*args, **kwargs),
        ),
    )


def make_single_selector(
    typename: str, Key: type[Semiring], Value: type[Semiring]
) -> Semiring:
    """
    Make a selector that keeps a single least or greatest value for a
    scalar idempotent semiring.

    :param name: name of the semiring type to create
    :param Key: scalar idempotent semiring (min-plus, max-plus, viterbi, ...)
    :param Value: semiring to augment
    """

    @dataclass(frozen=True, slots=True)
    class KeyValue:
        key: Key
        value: Value

    def select(item1, item2):
        key_sum = item1.key + item2.key

        if item1.key != item2.key:
            if item1.key == key_sum:
                return item1

            if item2.key == key_sum:
                return item2

        return KeyValue(key_sum, item1.value + item2.value)

    def combine(item1, item2):
        return KeyValue(item1.key * item2.key, item1.value * item2.value)

    return make_semiring(
        typename,
        cast=lambda item: KeyValue(Key.cast(item.key), Value.cast(item.value)),
        null=KeyValue(Key.null(), Value.null()),
        unit=KeyValue(Key.unit(), Value.unit()),
        add=select,
        mul=combine,
        make=lambda *args, **kwargs: KeyValue(
            Key.make(*args, **kwargs),
            Value.make(*args, **kwargs),
        ),
    )


def make_multiple_selector(
    typename: str,
    Key: type[Semiring],
    Value: type[Semiring],
) -> Semiring:
    """
    Make a selector that keeps all least or greatest values for an
    iterable partially-ordered semiring.

    :param name: name of the semiring type to create
    :param Key: iterable partially-ordered semiring (pareto, ...)
    :param Value: semiring to augment
    """

    def select(opts1, opts2):
        keys = Key(frozenset(opts1.keys())) + Key(frozenset(opts2.keys()))
        return Map(
            {
                key: opts1.get(key, Value.null()) + opts2.get(key, Value.null())
                for key in keys.value
            }
        )

    def combine(opts1, opts2):
        keys = Key(frozenset(opts1.keys())) * Key(frozenset(opts2.keys()))
        return Map(
            {
                key1 * key2: value1 * value2
                for (key1, value1), (key2, value2) in product(
                    opts1.items(), opts2.items()
                )
                if key1 * key2 in keys.value
            }
        )

    return make_semiring(
        typename,
        cast=lambda opts: Map({key: Value.cast(value) for key, value in opts.items()}),
        null=Map({}),
        unit=Map({key: Value.unit() for key in Key.unit().value}),
        add=select,
        mul=combine,
        make=lambda *args, **kwargs: Map(
            {
                key: Value.make(*args, **kwargs)
                for key in Key.make(*args, **kwargs).value
            }
        ),
    )
