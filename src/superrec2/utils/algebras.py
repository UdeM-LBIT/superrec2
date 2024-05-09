from math import inf
from immutables import Map
import operator
from typing import Generic, Callable, Iterable, TypeVar, ParamSpec, Self, Any


T = TypeVar("T")
U = TypeVar("U")
P = ParamSpec("P")


class SemiRingType(type):
    """Semiring metaclass defining class-level composition operations."""

    def __new__(meta, name, inherit, attrs):
        return super().__new__(meta, name, inherit, {**attrs, "__slots__": ("value",)})

    def __add__(first, second):
        """
        Create a semiring that is the direct product of multiple semirings.

        Values of the direct product are tuples of values of each combined
        semiring, and the addition and multiplication operations are
        broadcasted on each element of the tuple.

        This operation is non-commutative but associative. Combining
        more than two semi-rings yields tuples with more than two values.
        """
        args = [
            orig
            for arg in (first, second)
            for orig in (arg.__add_from__ if hasattr(arg, "__add_from__") else (arg,))
        ]

        def _add(x, y):
            return tuple(arg._add(x[i], y[i]) for i, arg in enumerate(args))

        def _mul(x, y):
            return tuple(arg._mul(x[i], y[i]) for i, arg in enumerate(args))

        result = SemiRingType(
            f'({" + ".join(arg.__name__ for arg in args)})',
            (SemiRing,),
            {
                "_zero": tuple(arg._zero for arg in args),
                "_one": tuple(arg._one for arg in args),
                "_add": _add,
                "_mul": _mul,
            },
        )
        result.__add_from__ = args
        return result

    def __mul__(first, second):
        """
        Create a semiring that is the selection product of two semirings.

        Values of the selection product are tuples (A, B) from the first and
        second semiring, such that A is minimum (as defined by the addition
        operation of the first semiring) and B is the sum of all encountered
        values associated to A.

        This operation is non-commutative and non-associative.
        """

        def _add(x, y):
            if x[0] == y[0]:
                return (x[0], second._add(x[1], y[1]))
            elif x[0] == first._add(x[0], y[0]):
                return (x[0], x[1])
            else:
                return (y[0], y[1])

        def _mul(x, y):
            return (first._mul(x[0], y[0]), second._mul(x[1], y[1]))

        return SemiRingType(
            f"({first.__name__} * {second.__name__})",
            (SemiRing,),
            {
                "_zero": (first._zero, second._zero),
                "_one": (first._one, second._one),
                "_add": _add,
                "_mul": _mul,
            },
        )

    def __matmul__(first, second):
        """
        Create a semiring that is the multiple selection product of two semirings.

        Values of the multiple selection product are dictionaries where keys
        are from the first semiring and values are from the second. The set of
        keys is always Pareto-minimum (as defined by the addition operation of
        the first semiring), and the values are sums of values of the second
        semiring for each key.

        This operation is non-commutative and non-associative.
        """

        def _add(x, y):
            keys = first._add(frozenset(x.keys()), frozenset(y.keys()))
            return Map(
                {
                    key: second._add(
                        x.get(key, second._zero),
                        y.get(key, second._zero),
                    )
                    for key in keys
                }
            )

        def _mul(x, y):
            keys = first._mul(frozenset(x.keys()), frozenset(y.keys()))
            return Map(
                {
                    key1 * key2: second._mul(value1, value2)
                    for key1, value1 in x.items()
                    for key2, value2 in y.items()
                    if key1 * key2 in keys
                }
            )

        return SemiRingType(
            f"({first.__name__} @ {second.__name__})",
            (SemiRing,),
            {
                "_zero": Map({left: second._zero for left in first._zero}),
                "_one": Map({left: second._one for left in first._one}),
                "_add": _add,
                "_mul": _mul,
            },
        )


class SemiRing(Generic[T], metaclass=SemiRingType):
    """
    Semiring type on immutable values.

    A semiring is made up of two operations, called addition and multiplication,
    both associative, with neutral elements zero and one, and such that
    multiplication distributes over addition.

    Values of the semiring are wrapped inside objects for which the addition
    and multiplication operations are overloaded to map to the corresponding
    semiring operations.

    Operations and neutral elements are to be defined on the ground type with an
    underscore prefix, a wrapped version will automatically be generated and
    exposed without the underscore prefix.

    Class-level operations are defined that automate the composition of
    multiple semirings together.
    """

    _pool: dict[T, Self] = {}

    zero: Self
    _zero: T

    one: Self
    _one: T

    value: T

    _add: Callable[[T, T], T]
    _mul: Callable[[T, T], T]

    def __init_subclass__(cls) -> None:
        cls._pool = {}
        cls.zero = cls(cls._zero)
        cls.one = cls(cls._one)

    def __new__(cls, value: T) -> "SemiRing[T]":
        # Return its argument unchanged to prevent nested wrapping
        if isinstance(value, cls):
            return value

        # Only instantiate if no instance with the same value exists
        if value not in cls._pool:
            cls._pool[value] = super().__new__(cls)

        return cls._pool[value]

    def __init__(self, value: T):
        if not hasattr(self, "value"):
            self.value = value

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.value})"

    def __eq__(self, other: Any) -> bool:
        if not hasattr(other, "value"):
            return NotImplemented

        return self.value == other.value

    def __hash__(self) -> int:
        return hash(self.value)

    def __add__(self, other: Any) -> Self:
        if isinstance(other, self.__class__):
            value = other.value
        else:
            value = other

        return self.__class__(self.__class__._add(self.value, value))

    def __mul__(self, other: Any) -> Self:
        if isinstance(other, self.__class__):
            value = other.value
        else:
            value = other

        return self.__class__(self.__class__._mul(self.value, value))


class Structure(Generic[T, P]):
    __slots__ = ("semiring", "morphism", "zero", "one", "__add_from__")

    def __init__(
        self,
        semiring: type[SemiRing[T]],
        morphism: Callable[P, T],
    ):
        self.semiring = semiring
        self.morphism = morphism
        self.zero = semiring.zero
        self.one = semiring.one
        self.__add_from__ = None

    def __call__(self, *args: P.args, **kwargs: P.kwargs) -> SemiRing[T]:
        return self.semiring(self.morphism(*args, **kwargs))

    def __add__(first, second):
        structs = [
            orig
            for struct in (first, second)
            for orig in (
                struct.__add_from__ if struct.__add_from__ is not None else (struct,)
            )
        ]

        def morphism(*args, **kwargs):
            return tuple(struct.morphism(*args, **kwargs) for struct in structs)

        result = Structure(first.semiring + second.semiring, morphism)
        result.__add_from__ = structs
        return result

    def __mul__(first, second):
        def morphism(*args, **kwargs):
            return (first.morphism(*args, **kwargs), second.morphism(*args, **kwargs))

        return Structure(first.semiring * second.semiring, morphism)

    def __matmul__(first, second):
        def morphism(*args, **kwargs):
            return Map(
                {
                    key: second.morphism(*args, **kwargs)
                    for key in first.morphism(*args, **kwargs)
                }
            )

        return Structure(first.semiring @ second.semiring, morphism)


class MinPlus(SemiRing[int | float]):
    _zero = inf
    _one = 0
    _add = min
    _mul = operator.add


class MaxPlus(SemiRing[int | float]):
    _zero = -inf
    _one = 0
    _add = max
    _mul = operator.add


class Viterbi(SemiRing[float]):
    _zero = 0
    _one = 1
    _add = max
    _mul = operator.mul


class Boolean(SemiRing[bool]):
    _zero = False
    _one = True
    _add = operator.or_
    _mul = operator.and_


class Counter(SemiRing[int]):
    _zero = 0
    _one = 1
    _add = operator.add
    _mul = operator.mul


def generator_of(unit: T) -> SemiRing[frozenset[T]]:
    class Generator(SemiRing[frozenset[T]]):
        _zero = frozenset()
        _one = frozenset({unit})
        _add = operator.or_

        def _mul(x, y):
            return frozenset({a * b for a in x for b in y})

    return Generator


def projector_of(unit: T) -> SemiRing[T | None]:
    class Projector(SemiRing):
        _zero = None
        _one = ()

        def _add(x, y):
            return x if x is not None else y

        def _mul(x, y):
            return x * y if x is not None and y is not None else None

    return Projector


def pareto_min(items: Iterable[T]) -> frozenset[T]:
    return frozenset(x for x in items if not any(y < x for y in items))


def pareto_of(unit: T) -> type[SemiRing[frozenset[T]]]:
    class Pareto(SemiRing[frozenset[T]]):
        _zero = frozenset()
        _one = frozenset({unit})

        def _add(x, y):
            return pareto_min(x | y)

        def _mul(x, y):
            return pareto_min(frozenset(a * b for a in x for b in y))

    return Pareto


def vector(cls):
    """Mixin for tuple classes to make the tuples ordered and combinable."""

    def __le__(self, other):
        return all(left <= right for left, right in zip(self, other))

    cls.__le__ = __le__

    def __lt__(self, other):
        return self <= other and self != other

    cls.__lt__ = __lt__

    def __ge__(self, other):
        return other <= self

    cls.__ge__ = __ge__

    def __gt__(self, other):
        return other <= self and self != other

    cls.__gt__ = __gt__

    def __mul__(self, other):
        return self.__class__(*(left + right for left, right in zip(self, other)))

    cls.__mul__ = __mul__
    return cls
