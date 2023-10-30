from typing import Any, Self, Generic, Protocol, Callable, TypeVar
from dataclasses import dataclass
from functools import wraps
from itertools import product
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


class UnitMagma(Generic[T], Protocol):
    """Structure with an operation and a unit element."""

    value: T

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
        null=frozenset(),
        unit=frozenset({Solution.unit()}),
        add=operator.or_,
        mul=lambda set1, set2: frozenset(
            {value1 * value2 for value1, value2 in product(set1, set2)}
        ),
        make=lambda *args, **kwargs: frozenset({Solution.make(*args, **kwargs)}),
    )


def min_plus(typename: str, make: Callable[..., T] = lambda x: x) -> Semiring[T]:
    """Semiring used for finding the minimum cost of any solution."""
    return make_semiring(
        typename,
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
        null=float("-inf"),
        unit=0,
        add=max,
        mul=operator.add,
        make=make,
    )


def viterbi(typename: str, make: Callable[..., T] = lambda x: x) -> Semiring[T]:
    """Semiring used for finding the maximum probability of any solution."""
    return make_semiring(
        typename,
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
        null=(Left.null(), Right.null()),
        unit=(Left.unit(), Right.unit()),
        add=lambda x, y: (x[0] + y[0], x[1] + y[1]),
        mul=lambda x, y: (x[0] * y[0], x[1] * y[1]),
        make=lambda *args, **kwargs: (
            Left.make(*args, **kwargs),
            Right.make(*args, **kwargs),
        ),
    )


def make_selector(
    typename: str, Cost: type[Semiring], Value: type[Semiring]
) -> Semiring:
    """
    Make a selector that keeps only least- or greatest-cost values of a semiring.

    :param name: name of the semiring type to create
    :param Cost: additively-idempotent semiring (min-plus, max-plus, viterbi, ...)
    :param Value: semiring to augment
    """

    @dataclass(frozen=True, slots=True)
    class CostValue:
        cost: Cost
        selected: Value

    def select(item1, item2):
        cost_sum = item1.cost + item2.cost

        if item1.cost != item2.cost:
            if item1.cost == cost_sum:
                return item1

            if item2.cost == cost_sum:
                return item2

        return CostValue(cost_sum, item1.selected + item2.selected)

    def combine(item1, item2):
        return CostValue(item1.cost * item2.cost, item1.selected * item2.selected)

    return make_semiring(
        typename,
        null=CostValue(Cost.null(), Value.null()),
        unit=CostValue(Cost.unit(), Value.unit()),
        add=select,
        mul=combine,
        make=lambda *args, **kwargs: CostValue(
            Cost.make(*args, **kwargs),
            Value.make(*args, **kwargs),
        ),
    )
