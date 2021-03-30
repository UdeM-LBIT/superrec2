"""Fast disjoint-set structure implementing the union-find algorithm."""


class UnionFind:
    """Fast disjoint-set structure implementing the union-find algorithm."""

    def __init__(self, n):
        """
        Create a disjoint-set structure.

        Elements in the set range from 0 to :param:`n` − 1. Initially, they
        all belong to a different singleton set.

        :param n: number of elements in the structure
        """
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, element: int) -> int:
        """
        Identify the set to which an element belongs.

        Complexity: O(α(n))

        :param element: element to test
        :returns: canonical element representing the set to which
            :param:`element` belongs
        """
        if self.parent[element] == element:
            return element

        self.parent[element] = self.find(self.parent[element])
        return self.parent[element]

    def unite(self, first: int, second: int) -> bool:
        """
        Unite the two sets to which two elements belong.

        Complexity: O(α(n))

        :param first: first element to unite
        :param second: second element to unite
        :returns: false if the two elements were already in the same set,
            true if two different sets actually got merged together
        """
        rep_first = self.find(first)
        rep_second = self.find(second)

        if rep_first == rep_second:
            return False

        if self.rank[rep_first] == self.rank[rep_second]:
            self.rank[rep_first] += 1
            self.parent[rep_second] = rep_first
        elif self.rank[rep_first] > self.rank[rep_second]:
            self.parent[rep_second] = rep_first
        else:
            self.parent[rep_first] = rep_second

        return True
