"""Fast disjoint-set structure implementing the union-find strategy."""
from typing import List, Optional
from copy import deepcopy


class DisjointSet:
    """Fast disjoint-set structure implementing the union-find strategy."""

    def __init__(self, count):
        """
        Create a disjoint-set structure.

        Elements in the set range from 0 to :param:`count` − 1. Initially, they
        all belong to a different singleton set.

        :param count: number of elements in the structure
        """
        self.parent = list(range(count))
        self.rank = [0] * count
        self.groups = count

    def find(self, element: int) -> int:
        """
        Identify the set to which an element belongs.

        Complexity: O(α(count))

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

        Complexity: O(α(count))

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

        self.groups -= 1
        return True

    def __len__(self) -> int:
        """Get the number of groups in this partition."""
        return self.groups

    def to_list(self) -> List[List[int]]:
        """Create a list of all groups in this partition."""
        result: List[List[int]] = [[] for i in range(len(self.parent))]

        for i in range(len(self.parent)):
            result[self.find(i)].append(i)

        return [group for group in result if group]

    def __repr__(self) -> str:
        return (
            "DisjointSet({"
            + ", ".join(f"{{{', '.join(map(str, group))}}}" for group in self.to_list())
            + "})"
        )

    def binary(self) -> List["DisjointSet"]:
        """
        Generate all possible partitions obtained from the current partition
        by merging groups so that exactly two groups remain.
        """

        def _binary(
            partition: "DisjointSet",
            groups: List[int],
            first: Optional[int],
            second: Optional[int],
        ) -> List["DisjointSet"]:
            if not groups:
                if first is None or second is None:
                    return []

                return [partition]

            part_1 = deepcopy(partition)

            if first is not None:
                part_1.unite(first, groups[0])
                results_1 = _binary(part_1, groups[1:], first, second)
            elif second is None or groups[0] < second:
                results_1 = _binary(part_1, groups[1:], groups[0], second)
            else:
                results_1 = []

            part_2 = deepcopy(partition)

            if second is not None:
                part_2.unite(second, groups[0])
                results_2 = _binary(part_2, groups[1:], first, second)
            elif first is None or groups[0] > first:
                results_2 = _binary(part_2, groups[1:], first, groups[0])
            else:
                results_2 = []

            return results_1 + results_2

        return _binary(
            partition=self,
            groups=list(set(self.find(i) for i in range(len(self.parent)))),
            first=None,
            second=None,
        )
