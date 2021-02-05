from .range_min_query import RangeMinQuery


def _euler_tour(root, level=0):
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
    def __init__(self, tree):
        """
        Pre-compute the sparse table for lowest common ancestor queries.

        Complexity: O(V × log(V)), with V the number of nodes in `tree`.

        :param tree: input tree to make requests on
        """
        self.tree = tree
        self.traversal = _euler_tour(tree)
        self.range_min_query = RangeMinQuery(self.traversal)
        self.traversal_index = {}

        for i, (_, node) in enumerate(self.traversal):
            if node not in self.traversal_index:
                self.traversal_index[node] = i

    def __call__(self, a, b):
        """
        Find the lowest common ancestor of two nodes.

        Complexity: O(1).

        :param a: first node
        :param b: second node
        :returns: the ancestor of both `a` and `b` that is most distant
            from the tree root
        """
        a_index = self.traversal_index[a]
        b_index = self.traversal_index[b]

        if a_index <= b_index:
            return self.range_min_query(a_index, b_index + 1)[1]
        else:
            return self.range_min_query(b_index, a_index + 1)[1]

    def is_ancestor_of(self, a, b):
        """
        Check whether a node is an ancestor of another.

        Complexity: O(1).

        :param a: ancestor node
        :param b: descendant node
        :returns: True if and only if `b` is on the path from the tree root
            to `a` (i.e., `a` is an ancestor of `b`)
        """
        return self(a, b) == a

    def is_strict_ancestor_of(self, a, b):
        """
        Check whether a node is a strict an ancestor of another
        (i.e. is an ancestor distinct from the other node).

        Complexity: O(1).

        :param a: ancestor node
        :param b: descendant node
        :returns: True if and only if `b` is on the path from the tree root
            to `a` (i.e., `a` is an ancestor of `b`)
        """
        return self(a, b) == a and a != b

    def is_comparable(self, a, b):
        """
        Check whether two nodes are in the same subtree.

        Complexity: O(1).

        :param a: first node
        :param b: second node
        :returns True if and only if either `a` is an ancestor of `b` or
            or descendant of `b` (i.e., `a` and `b` are comparable)
        """
        return self.is_ancestor_of(a, b) or self.is_ancestor_of(b, a)

    def level(self, node):
        """
        Find the level of a node in the tree.

        Complexity: O(1).
        """
        return self.traversal[self.traversal_index[node]][0]

    def distance(self, a, b):
        """
        Find the distance between two nodes in the tree.

        Complexity: O(1).
        """
        return self.level(a) + self.level(b) - 2 * self.level(self(a, b))
