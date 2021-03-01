from typing import Generic, List, Tuple, TypeVar, Union
from collections.abc import Sequence
from numbers import Integral
from infinity import Infinity, inf


T = TypeVar('T')
ExtendedIntegral = Union[Integral, Infinity]


class MinSequence(Generic[T], Sequence[T]):
    """Receive ordered elements and keep only the minimal ones."""

    def __init__(self, max_keep: ExtendedIntegral = inf) -> None:
        """
        Create a minimum sequence.
        :param max_keep: Maximum count of minimal values to keep.
        """
        self.min: ExtendedIntegral = inf
        self.max_keep = max_keep
        self._items: List[T] = []

    def update(self, *elements: Tuple[ExtendedIntegral, T]) -> None:
        """Receive a set of elements."""
        for element in elements:
            if element[0] < inf:
                if self.min > element[0]:
                    self.min = element[0]
                    self._items = [element[1]]
                elif self.min == element[0]:
                    if len(self._items) < self.max_keep:
                        self._items.append(element[1])

    def __getitem__(self, idx: int) -> T:
        """Get the nth minimal element."""
        return self._items[idx]

    def __len__(self) -> int:
        """Get the number of minimal elements."""
        return len(self._items)
