import unittest
from ete3 import Tree
from .lowest_common_ancestor import LowestCommonAncestor


def _is_ancestor_of_naive(a, b):
    while b.up != None and a != b:
        b = b.up
    return a == b


def _get_common_ancestor_naive(a, b):
    while a.up != None:
        b_ancestor = b

        while b_ancestor.up != None and a != b_ancestor:
            b_ancestor = b_ancestor.up

        if a == b_ancestor:
            return a

        a = a.up

    return a

def _level_naive(node):
    level = 0

    while node.up != None:
        level += 1
        node = node.up

    return level


def _distance_naive(start, stop, visited=None):
    if visited is None:
        visited = set((start,))
    else:
        visited.add(start)

    if start == stop:
        return 0

    distance = float('inf')

    if start.up is not None and start.up not in visited:
        distance = min(distance, 1 + _distance_naive(start.up, stop, visited))

    for child in start.children:
        if child not in visited:
            distance = min(distance, 1 + _distance_naive(child, stop, visited))

    return distance


class TestLowestCommonAncestor(unittest.TestCase):
    def test_lca(self):
        tree = Tree("((2,(4,5)3)1,(7,8,(10)9)6)0;", format=8)
        lca = LowestCommonAncestor(tree)

        for a in tree.traverse():
            for b in tree.traverse():
                ab_lca = lca(a, b)
                self.assertEqual(ab_lca, tree.get_common_ancestor(a, b))
                self.assertEqual(ab_lca, _get_common_ancestor_naive(a, b))

    def test_ancestor_relation(self):
        tree = Tree("((2,(4,5)3)1,(7,8,(10)9)6)0;", format=8)
        lca = LowestCommonAncestor(tree)

        for a in tree.traverse():
            for b in tree.traverse():
                ab_lca = lca(a, b)

                self.assertTrue(a in ab_lca or a == ab_lca)
                self.assertTrue(b in ab_lca or b == ab_lca)
                self.assertEqual(lca.is_ancestor_of(a, b), b in a or a == b)
                self.assertEqual(lca.is_ancestor_of(b, a), a in b or a == b)
                self.assertEqual(lca.is_strict_ancestor_of(a, b), b in a)
                self.assertEqual(lca.is_strict_ancestor_of(b, a), a in b)
                self.assertEqual(
                    lca.is_comparable(a, b),
                    a in b or b in a or a == b
                )

                self.assertTrue(_is_ancestor_of_naive(ab_lca, a))
                self.assertTrue(_is_ancestor_of_naive(ab_lca, b))
                self.assertEqual(
                    lca.is_ancestor_of(a, b),
                    _is_ancestor_of_naive(a, b)
                )
                self.assertEqual(
                    lca.is_ancestor_of(b, a),
                    _is_ancestor_of_naive(b, a)
                )
                self.assertEqual(
                    lca.is_strict_ancestor_of(a, b),
                    _is_ancestor_of_naive(a, b) and a != b
                )
                self.assertEqual(
                    lca.is_strict_ancestor_of(b, a),
                    _is_ancestor_of_naive(b, a) and b != a
                )
                self.assertEqual(
                    lca.is_comparable(a, b),
                    _is_ancestor_of_naive(a, b) or _is_ancestor_of_naive(b, a)
                )

    def test_level(self):
        tree = Tree("((2,(4,5)3)1,(7,8,(10)9)6)0;", format=8)
        lca = LowestCommonAncestor(tree)

        for node in tree.traverse():
            self.assertEqual(lca.level(node), _level_naive(node))

    def test_distance(self):
        tree = Tree("((2,(4,5)3)1,(7,8,(10)9)6)0;", format=8)
        lca = LowestCommonAncestor(tree)

        for a in tree.traverse():
            for b in tree.traverse():
                self.assertEqual(
                    _distance_naive(a, b),
                    lca.distance(a, b),
                )
