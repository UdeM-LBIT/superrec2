def _ilog2(integer):
    """Integral part of the base-2 logarithm of a positive integer."""
    return integer.bit_length() - 1


class RangeMinQuery:
    """
    Fast answering of repeated range-minimum queries using a sparse table.

    Note that the input data cannot be changed after initialization.
    See <https://cp-algorithms.com/data_structures/sparse-table.html>.
    """

    def __init__(self, data):
        """
        Pre-compute the sparse table for range-minimum queries.

        Complexity: O(N × log(N)), where N = len(data).

        :param data: input list of objects with total ordering
        """
        length = len(data)
        levels = _ilog2(length) + 1

        # sparse_table[depth][i] stores the minimum of the
        # (i, i + 2**depth) range
        self.sparse_table = [[None] * length for _ in range(levels)]
        self.sparse_table[0] = data.copy()

        for depth in range(1, levels):
            for i in range(length - 2 ** depth + 1):
                self.sparse_table[depth][i] = min(
                    self.sparse_table[depth - 1][i],
                    self.sparse_table[depth - 1][i + 2 ** (depth - 1)],
                )

    def __call__(self, start, stop):
        """
        Find the minimum value in a range of the input data.

        Complexity: O(1).

        :param start: first index of the range
        :param stop: index following the last index of the range
        :returns: minimum value
        """
        if start >= stop:
            return float('inf')

        depth = _ilog2(stop - start)
        return min(
            self.sparse_table[depth][start],
            self.sparse_table[depth][stop - 2 ** depth],
        )
