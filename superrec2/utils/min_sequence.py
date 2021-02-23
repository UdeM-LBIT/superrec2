from typing import Generic, List, Tuple, TypeVar, Union
from collections.abc import Sequence


T = TypeVar('T')


class MinSequence(Generic[T], Sequence[T]):
    """Receive ordered elements and keep only the minimal ones."""

    def __init__(self) -> None:
        self.min: Union[int, float] = float('inf')
        self._items: List[T] = []

    def update(self, *elements: Tuple[int, T]) -> None:
        """Receive a set of elements."""
        for element in elements:
            if element[0] < float('inf'):
                if self.min > element[0]:
                    self.min = element[0]
                    self._items = [element[1]]
                elif self.min == element[0]:
                    self._items.append(element[1])

    def __getitem__(self, idx: int) -> T:
        """Get the nth minimal element."""
        return self._items[idx]

    def __len__(self) -> int:
        """Get the number of minimal elements."""
        return len(self._items)
