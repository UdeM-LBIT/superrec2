"""Common operations on trees."""
from itertools import product
from typing import Dict, List, Optional, Tuple
from ete3 import Tree, TreeNode
from .range_min_query import RangeMinQuery
from .disjoint_set import DisjointSet


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


Triple = Tuple[str, str, str]


def tree_to_triples(tree: Tree) -> Tuple[List[str], List[Triple]]:
    """
    Compute a set of triples that encode the topology of a binary
    phylogenetic tree.

    A triple is a binary tree containing exactly three leaves and two internal
    nodes. Each triple is represented as a canonical tuple containing its three
    leaves in the following order:

    * the deepest leaf that comes first in lexicographic order
    * the deepest leaf that comes second
    * the shallowest leaf

    This implements the BreakUp algorithm from [Ng and Wormald, 1996],
    restricted to binary trees and triples only. Feeding the output of this
    function directly to :func:`tree_from_triples` will reconstruct the
    original tree.

    :param tree: input tree to break up in triples
    :returns: a tuple containing the labels of leaves in the original tree and
        the list of triples that characterize the tree’s topology
    """
    result = []
    leaves = [leaf.name for leaf in tree.get_leaves()]
    tree = tree.copy()
    minimal_int_nodes = set()

    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue

        if not all(child.is_leaf() for child in node.children):
            continue

        minimal_int_nodes.add(node)

    while minimal_int_nodes:
        other = minimal_int_nodes.pop()
        parent = other.up

        if parent is None:
            break

        leaf = other.get_sisters()[0].get_leaves()[0].name
        left_node, right_node = other.children
        left, right = left_node.name, right_node.name

        parent.add_child(right_node.detach())
        parent.remove_child(other)

        if left <= right:
            result.append((left, right, leaf))
        else:
            result.append((right, left, leaf))

        if all(child.is_leaf() for child in parent.children):
            minimal_int_nodes.add(parent)

    return leaves, result


def tree_from_triples(
    leaves: List[str], triples: List[Triple]
) -> Optional[Tree]:
    """
    Reconstruct a phylogenetic tree that respects the constraints given
    by a set of triples.

    Note that the reconstructed tree may not be a binary tree, if the set of
    triples is not specific enough.

    This implements the OneTree algorithm from [Ng and Wormald, 1996],
    restricted to triples only.

    :param leaves: set of leaf labels
    :param triples: set of triples that constrain the tree topology
    :returns: either a phylogenetic tree, or None if the set of triples
        is not consistent
    """
    if not leaves:
        return None

    if len(leaves) == 1:
        return Tree(name=leaves[0])

    if len(leaves) == 2:
        root = Tree()
        root.add_child(Tree(name=leaves[0]))
        root.add_child(Tree(name=leaves[1]))
        return root

    partition = DisjointSet(len(leaves))
    leaf_index = {leaf: i for i, leaf in enumerate(leaves)}

    for left, right, _ in triples:
        partition.unite(leaf_index[left], leaf_index[right])

    if len(partition) <= 1:
        return None

    root = Tree()

    for group in partition.to_list():
        group_leaves = [leaves[item] for item in group]
        group_triples = [
            triple
            for triple in triples
            if all(leaf in group_leaves for leaf in triple)
        ]

        subtree = tree_from_triples(group_leaves, group_triples)

        if not subtree:
            return None

        root.add_child(subtree)

    return root


def all_trees_from_triples(
    leaves: List[str], triples: List[Triple]
) -> List[Tree]:
    """
    Find all the phylogenetic binary trees that respect the constraints given
    by a set of triples.

    This implements the AllTrees algorithm from [Ng and Wormald, 1996],
    restricted to triples and binary trees only.

    :param leaves: set of leaf labels
    :param triples: set of triples that constrain the tree topology
    :returns: list of compatible either a phylogenetic trees, potentially
        empty if the set of triples is not consistent
    """

    def _all_trees_from_triples(
        leaves: List[str], triples: List[Triple]
    ) -> List[Tree]:
        if not leaves:
            return []

        if len(leaves) == 1:
            return [Tree(name=leaves[0])]

        if len(leaves) == 2:
            root = Tree()
            root.add_child(Tree(name=leaves[0]))
            root.add_child(Tree(name=leaves[1]))
            return [root]

        partition = DisjointSet(len(leaves))
        leaf_index = {leaf: index for index, leaf in enumerate(leaves)}

        for left, right, _ in triples:
            partition.unite(leaf_index[left], leaf_index[right])

        results = []

        for bin_partition in partition.binary():
            groups = bin_partition.to_list()
            groups_leaves = [
                [leaves[item] for item in group] for group in groups
            ]
            groups_triples = [
                [
                    triple
                    for triple in triples
                    if all(leaf in group_leaves for leaf in triple)
                ]
                for group_leaves in groups_leaves
            ]

            for left_tree, right_tree in product(
                _all_trees_from_triples(groups_leaves[0], groups_triples[0]),
                _all_trees_from_triples(groups_leaves[1], groups_triples[1]),
            ):
                root = Tree()
                root.add_child(left_tree.copy())
                root.add_child(right_tree.copy())
                results.append(root)

        return results

    # Detect whether the input set of triples is inconsistent right away
    if tree_from_triples(leaves, triples) is None:
        return []

    return _all_trees_from_triples(leaves, triples)
