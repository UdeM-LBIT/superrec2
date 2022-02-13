"""Preprocessing to answer minimum-in-range queries in constant time."""
from typing import Any, Generic, List, Optional, Protocol, Sequence, TypeVar


def _ilog2(value: int):
    """Integral part of the base-2 logarithm of a positive integer."""
    return value.bit_length() - 1


# See <https://github.com/python/typing/issues/760>
class SupportsLessThan(Protocol):  # pylint:disable=too-few-public-methods
    """Types that support the less-than operator."""

    def __lt__(self, other: Any) -> bool:
        ...


Element = TypeVar("Element", bound=SupportsLessThan)


class RangeMinQuery(Generic[Element]):  # pylint:disable=too-few-public-methods
    """
    Fast answering of repeated range-minimum queries using a sparse table.

    Note that the input data cannot be changed after initialization.
    See <https://cp-algorithms.com/data_structures/sparse-table.html>.
    """

    def __init__(self, data: Sequence[Element]):
        """
        Pre-compute the sparse table for range-minimum queries.

        Complexity: O(N × log(N)), where N = len(data).

        :param data: input list of objects with total ordering
        """
        length = len(data)
        levels = _ilog2(length) + 1

        # sparse_table[depth][i] stores the minimum of the
        # (i, i + 2**depth) range
        self.sparse_table: List[List[Optional[Element]]] = [
            [None] * length for _ in range(levels)
        ]
        self.sparse_table[0] = list(data)

        for depth in range(1, levels):
            for i in range(length - 2**depth + 1):
                left = self.sparse_table[depth - 1][i]
                assert left is not None
                right = self.sparse_table[depth - 1][i + 2 ** (depth - 1)]
                assert right is not None
                self.sparse_table[depth][i] = min(left, right)

    def __call__(self, start: int, stop: int) -> Optional[Element]:
        """
        Find the minimum value in a range of the input data.

        Complexity: O(1).

        :param start: first index of the range
        :param stop: index following the last index of the range
        :returns: minimum value, or None if the range is empty
        """
        if start >= stop:
            return None

        depth = _ilog2(stop - start)
        return min(
            self.sparse_table[depth][start],
            self.sparse_table[depth][stop - 2**depth],
        )
