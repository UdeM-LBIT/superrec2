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


class TestLowestCommonAncestor(unittest.TestCase):
    def test_valid_lca(self):
        tree = Tree("((2,(4,5)3)1,(7,8,(10)9)6)0;", format=8)
        lca = LowestCommonAncestor(tree)

        for a in tree.traverse():
            for b in tree.traverse():
                ab_lca = lca(a, b)

                # Check with native ete3 operators
                self.assertTrue(a in ab_lca or a == ab_lca)
                self.assertTrue(b in ab_lca or b == ab_lca)
                self.assertEqual(ab_lca, tree.get_common_ancestor(a, b))
                self.assertEqual(lca.is_ancestor_of(a, b), b in a or a == b)
                self.assertEqual(lca.is_ancestor_of(b, a), a in b or a == b)
                self.assertEqual(
                    lca.is_comparable(a, b),
                    a in b or b in a or a == b
                )

                # Check with naive implementation
                self.assertTrue(_is_ancestor_of_naive(ab_lca, a))
                self.assertTrue(_is_ancestor_of_naive(ab_lca, b))
                self.assertEqual(ab_lca, _get_common_ancestor_naive(a, b))
                self.assertEqual(
                    lca.is_ancestor_of(a, b),
                    _is_ancestor_of_naive(a, b)
                )
                self.assertEqual(
                    lca.is_ancestor_of(b, a),
                    _is_ancestor_of_naive(b, a)
                )
                self.assertEqual(
                    lca.is_comparable(a, b),
                    _is_ancestor_of_naive(a, b) or _is_ancestor_of_naive(b, a)
                )

