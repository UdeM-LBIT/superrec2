from typing import Generic, List, overload, Tuple, TypeVar, Union
from collections.abc import Sequence
from infinity import Infinity, inf


T = TypeVar('T')
ExtendedIntegral = Union[int, Infinity]


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

    @overload
    def __getitem__(self, idx: int) -> T:
        ...

    @overload
    def __getitem__(self, idx: slice) -> "MinSequence[T]":
        ...

    def __getitem__(self, idx: Union[int, slice]) -> Union[T, "MinSequence[T]"]:
        """Get the nth minimal element."""
        if isinstance(idx, slice):
            result: MinSequence[T] = MinSequence(self.max_keep)
            result._items = self._items[idx]
            return result
        else:
            return self._items[idx]

    def __len__(self) -> int:
        """Get the number of minimal elements."""
        return len(self._items)
