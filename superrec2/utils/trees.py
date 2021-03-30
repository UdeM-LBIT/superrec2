"""Preprocessing to answer lowest common ancestor queries in constant time."""
from typing import Dict, List, Tuple
from ete3 import Tree, TreeNode
from .range_min_query import RangeMinQuery


def _euler_tour(root: TreeNode, level: int = 0) -> List[Tuple[int, TreeNode]]:
    """
    Compute the Euler tour representation of a tree.

    See <https://en.wikipedia.org/wiki/Euler_tour_technique>.

    :param root: root node of the tree
    :param level: current level (used for recursive calls)
    :returns: Euler tour representation
    """
    if root.is_leaf():
        return [(level, root)]

    tour = [(level, root)]

    for child in root.children:
        tour.extend(_euler_tour(child, level + 1))
        tour.append((level, root))

    return tour


class LowestCommonAncestor:
    """
    Fast answering of repeated lowest common ancestor queries using a
    sparse table range-minimum query implementation.

    Note that the input tree cannot be changed after initialization.
    See <https://cp-algorithms.com/graph/lca.html>.
    """

    def __init__(self, tree: Tree):
        """
        Pre-compute the sparse table for lowest common ancestor queries.

        Complexity: O(V × log(V)), with V the number of nodes in `tree`.

        :param tree: input tree to make requests on
        """
        self.tree = tree
        self.traversal = _euler_tour(tree)
        self.range_min_query = RangeMinQuery(self.traversal)
        self.traversal_index: Dict[TreeNode, int] = {}

        for i, (_, node) in enumerate(self.traversal):
            if node not in self.traversal_index:
                self.traversal_index[node] = i

    def __call__(self, first: TreeNode, second: TreeNode) -> TreeNode:
        """
        Find the lowest common ancestor of two nodes.

        Complexity: O(1).

        :param first: first node
        :param second: second node
        :returns: the ancestor of both `first` and `second` that is most
            distant from the tree root
        """
        a_index = self.traversal_index[first]
        b_index = self.traversal_index[second]

        if a_index <= b_index:
            result = self.range_min_query(a_index, b_index + 1)
        else:
            result = self.range_min_query(b_index, a_index + 1)

        assert result is not None
        return result[1]

    def is_ancestor_of(self, first: TreeNode, second: TreeNode) -> bool:
        """
        Check whether a node is an ancestor of another.

        Complexity: O(1).

        :param first: ancestor node
        :param second: descendant node
        :returns: True if and only if `second` is on the path from the tree
            root to `first` (i.e., `first` is an ancestor of `second`)
        """
        return self(first, second) == first

    def is_strict_ancestor_of(self, first: TreeNode, second: TreeNode) -> bool:
        """
        Check whether a node is a strict an ancestor of another
        (i.e. is an ancestor distinct from the other node).

        Complexity: O(1).

        :param first: ancestor node
        :param second: descendant node
        :returns: True if and only if `second` is on the path from the tree
            root to `first` (i.e., `first` is an ancestor of `second`)
        """
        return self(first, second) == first and first != second

    def is_comparable(self, first: TreeNode, second: TreeNode) -> bool:
        """
        Check whether two nodes are in the same subtree.

        Complexity: O(1).

        :param first: first node
        :param second: second node
        :returns: True if and only if either `first` is an ancestor of `second`
            or descendant of `second` (i.e., `first` and `second` are
            comparable)
        """
        return self.is_ancestor_of(
            first, second
        ) or self.is_ancestor_of(  # pylint:disable=arguments-out-of-order
            second, first
        )

    def level(self, node: TreeNode) -> int:
        """
        Find the level of a node in the tree.

        Complexity: O(1).
        """
        return self.traversal[self.traversal_index[node]][0]

    def distance(self, first: TreeNode, second: TreeNode) -> int:
        """
        Find the distance between two nodes in the tree.

        Complexity: O(1).
        """
        return (
            self.level(first)
            + self.level(second)
            - 2 * self.level(self(first, second))
        )
