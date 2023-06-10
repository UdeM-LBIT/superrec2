import unittest
from itertools import chain, combinations
from ete3 import Tree
from infinity import inf
from superrec2.utils.trees import (
    LowestCommonAncestor,
    tree_to_triples,
    tree_from_triples,
    all_trees_from_triples,
    supertree,
    all_supertrees,
    is_binary,
    graft,
    arrange_leaves,
    binarize,
)


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

    distance = inf

    if start.up is not None and start.up not in visited:
        distance = min(distance, 1 + _distance_naive(start.up, stop, visited))

    for child in start.children:
        if child not in visited:
            distance = min(distance, 1 + _distance_naive(child, stop, visited))

    return distance


def _powerset(iterable):
    seq = list(iterable)
    return chain.from_iterable(combinations(seq, r) for r in range(len(seq) + 1))


class TestLowestCommonAncestor(unittest.TestCase):
    def test_lca(self):
        tree = Tree("((2,(4,5)3)1,(7,8,(10)9)6)0;", format=8)
        lca = LowestCommonAncestor(tree)

        for a in tree.traverse():
            for b in tree.traverse():
                ab_lca = lca(a, b)
                self.assertEqual(ab_lca, tree.get_common_ancestor(a, b))
                self.assertEqual(ab_lca, _get_common_ancestor_naive(a, b))

        for subset in _powerset(tree.traverse()):
            if not subset:
                with self.assertRaises(TypeError):
                    subset_lca = lca(*subset)
            elif len(subset) == 1:
                self.assertEqual(lca(*subset), subset[0])
            else:
                subset_lca = lca(*subset)
                naive_lca = subset[0]

                for node in subset[1:]:
                    naive_lca = _get_common_ancestor_naive(naive_lca, node)

                self.assertEqual(subset_lca, tree.get_common_ancestor(*subset))
                self.assertEqual(subset_lca, naive_lca)

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
                self.assertEqual(lca.is_comparable(a, b), a in b or b in a or a == b)

                self.assertTrue(_is_ancestor_of_naive(ab_lca, a))
                self.assertTrue(_is_ancestor_of_naive(ab_lca, b))
                self.assertEqual(lca.is_ancestor_of(a, b), _is_ancestor_of_naive(a, b))
                self.assertEqual(lca.is_ancestor_of(b, a), _is_ancestor_of_naive(b, a))
                self.assertEqual(
                    lca.is_strict_ancestor_of(a, b),
                    _is_ancestor_of_naive(a, b) and a != b,
                )
                self.assertEqual(
                    lca.is_strict_ancestor_of(b, a),
                    _is_ancestor_of_naive(b, a) and b != a,
                )
                self.assertEqual(
                    lca.is_comparable(a, b),
                    _is_ancestor_of_naive(a, b) or _is_ancestor_of_naive(b, a),
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


class TestSupertree(unittest.TestCase):
    def test_tree_to_triples(self):
        tree1 = Tree("(a,(b,(c,(d,(e,f)))));")
        self.assertEqual(
            tree_to_triples(tree1),
            (
                list("abcdef"),
                [
                    ("e", "f", "d"),
                    ("d", "f", "c"),
                    ("c", "f", "b"),
                    ("b", "f", "a"),
                ],
            ),
        )

        tree2 = Tree("(a,b);")
        self.assertEqual(
            tree_to_triples(tree2),
            (
                list("ab"),
                [],
            ),
        )

    def test_tree_from_triples(self):
        self.assertEqual(
            tree_from_triples(
                list("abcd"),
                [
                    ("a", "b", "c"),
                    ("c", "d", "a"),
                ],
            ).write(format=9),
            "((a,b),(c,d));",
        )

        self.assertEqual(
            tree_from_triples(
                list("abc"),
                [
                    ("a", "b", "c"),
                    ("b", "c", "a"),
                ],
            ),
            None,
        )

    def test_all_trees_from_triples(self):
        self.assertCountEqual(
            (
                result.write(format=9)
                for result in all_trees_from_triples(
                    list("abcd"),
                    [
                        ("a", "b", "c"),
                        ("a", "b", "d"),
                    ],
                )
            ),
            (
                "(((a,b),c),d);",
                "(((a,b),d),c);",
                "((a,b),(c,d));",
            ),
        )
        self.assertEqual(
            all_trees_from_triples(
                list("abc"),
                [
                    ("a", "b", "c"),
                    ("b", "c", "a"),
                ],
            ),
            [],
        )

    def test_tree_feedback(self):
        tree = Tree("((a,(b,c)),((d,e),(f,g)));")
        self.assertEqual(
            tree.get_topology_id(),
            tree_from_triples(*tree_to_triples(tree)).get_topology_id(),
        )
        self.assertCountEqual(
            (
                result.get_topology_id()
                for result in all_trees_from_triples(*tree_to_triples(tree))
            ),
            (tree.get_topology_id(),),
        )

    def test_supertree(self):
        trees = (
            Tree("((Y1,Z1),(X2,X3));"),
            Tree("((Z1,X1),X2);"),
            Tree("(X2,X3);"),
        )

        self.assertEqual(
            supertree(trees).get_topology_id(),
            Tree("((Z1,Y1,X1),(X3,X2));").get_topology_id(),
        )
        self.assertCountEqual(
            (result.get_topology_id() for result in all_supertrees(trees)),
            (
                Tree("(((Z1,X1),Y1),(X3,X2));").get_topology_id(),
                Tree("(((Z1,Y1),X1),(X3,X2));").get_topology_id(),
                Tree("((Z1,(X1,Y1)),(X3,X2));").get_topology_id(),
            ),
        )


class TestBinary(unittest.TestCase):
    def assertTreeResults(self, expect, actual):
        self.assertCountEqual(
            list(map(lambda tree: tree.get_topology_id(), expect)),
            list(map(lambda tree: tree.get_topology_id(), actual)),
        )

    def test_is_binary(self):
        self.assertTrue(is_binary(Tree("A;")))
        self.assertTrue(is_binary(Tree("(A,B);")))
        self.assertTrue(is_binary(Tree("(A,(B,C));")))
        self.assertTrue(is_binary(Tree("((A,B),C);")))
        self.assertTrue(is_binary(Tree("((A,B),(C,D));")))
        self.assertTrue(is_binary(Tree("((((A,B),C),D),E);")))
        self.assertTrue(is_binary(Tree("(A,(B,(C,(D,E))));")))
        self.assertTrue(is_binary(Tree("(A,(B,(C,(D,E))));")))
        self.assertFalse(is_binary(Tree("(A);")))
        self.assertFalse(is_binary(Tree("(A,B,C);")))
        self.assertFalse(is_binary(Tree("(((((A,B,C),D),E),F),G);")))
        self.assertFalse(is_binary(Tree("(A,(B,(C,(D,(E,F,G)))));")))
        self.assertFalse(is_binary(Tree("(A,(B,(C,(D),E),F),G);")))

    def test_graft(self):
        self.assertTreeResults(
            graft(Tree("A;"), Tree("X;")),
            (Tree("(A,X);"),),
        )
        self.assertTreeResults(
            graft(Tree("(A,B);"), Tree("X;")),
            (Tree("((A,X),B);"), Tree("(A,(B,X));"), Tree("((A,B),X);")),
        )
        self.assertTreeResults(
            graft(Tree("(A,(B,C));"), Tree("X;")),
            (
                Tree("(X,(A,(B,C)));"),
                Tree("((A,X),(B,C));"),
                Tree("(A,(X,(B,C)));"),
                Tree("(A,((B,X),C));"),
                Tree("(A,(B,(C,X)));"),
            ),
        )

    def test_arrange_leaves(self):
        self.assertTreeResults(arrange_leaves(()), ())
        self.assertTreeResults(arrange_leaves((Tree("A;"),)), (Tree("A;"),))
        self.assertTreeResults(
            arrange_leaves((Tree("A;"), Tree("B;"))), (Tree("(A,B);"),)
        )
        self.assertTreeResults(
            arrange_leaves((Tree("A;"), Tree("B;"), Tree("C;"))),
            (
                Tree("(A,(B,C));"),
                Tree("(B,(A,C));"),
                Tree("(C,(A,B));"),
            ),
        )

    def test_binarize(self):
        self.assertTreeResults(
            binarize(Tree("(((A,B),(C,D,E,F)P,G),H);", format=1)),
            (
                Tree("((((A,B),(((C,D),E),F)P),G),H);", format=1),
                Tree("((((A,B),(((C,D),F),E)P),G),H);", format=1),
                Tree("((((A,B),((C,D),(E,F))P),G),H);", format=1),
                Tree("((((A,B),(((C,E),D),F)P),G),H);", format=1),
                Tree("((((A,B),(((C,E),F),D)P),G),H);", format=1),
                Tree("((((A,B),((C,E),(D,F))P),G),H);", format=1),
                Tree("((((A,B),(((C,F),D),E)P),G),H);", format=1),
                Tree("((((A,B),(((C,F),E),D)P),G),H);", format=1),
                Tree("((((A,B),((C,F),(E,D))P),G),H);", format=1),
                Tree("((((A,B),(((D,E),C),F)P),G),H);", format=1),
                Tree("((((A,B),(((D,E),F),C)P),G),H);", format=1),
                Tree("((((A,B),(((D,F),C),E)P),G),H);", format=1),
                Tree("((((A,B),(((D,F),E),C)P),G),H);", format=1),
                Tree("((((A,B),(((E,F),C),D)P),G),H);", format=1),
                Tree("((((A,B),(((E,F),D),C)P),G),H);", format=1),
                Tree("((((A,B),G),(((C,D),E),F)P),H);", format=1),
                Tree("((((A,B),G),(((C,D),F),E)P),H);", format=1),
                Tree("((((A,B),G),((C,D),(E,F))P),H);", format=1),
                Tree("((((A,B),G),(((C,E),D),F)P),H);", format=1),
                Tree("((((A,B),G),(((C,E),F),D)P),H);", format=1),
                Tree("((((A,B),G),((C,E),(D,F))P),H);", format=1),
                Tree("((((A,B),G),(((C,F),D),E)P),H);", format=1),
                Tree("((((A,B),G),(((C,F),E),D)P),H);", format=1),
                Tree("((((A,B),G),((C,F),(E,D))P),H);", format=1),
                Tree("((((A,B),G),(((D,E),C),F)P),H);", format=1),
                Tree("((((A,B),G),(((D,E),F),C)P),H);", format=1),
                Tree("((((A,B),G),(((D,F),C),E)P),H);", format=1),
                Tree("((((A,B),G),(((D,F),E),C)P),H);", format=1),
                Tree("((((A,B),G),(((E,F),C),D)P),H);", format=1),
                Tree("((((A,B),G),(((E,F),D),C)P),H);", format=1),
                Tree("(((A,B),((((C,D),E),F)P,G)),H);", format=1),
                Tree("(((A,B),((((C,D),F),E)P,G)),H);", format=1),
                Tree("(((A,B),(((C,D),(E,F))P,G)),H);", format=1),
                Tree("(((A,B),((((C,E),D),F)P,G)),H);", format=1),
                Tree("(((A,B),((((C,E),F),D)P,G)),H);", format=1),
                Tree("(((A,B),(((C,E),(D,F))P,G)),H);", format=1),
                Tree("(((A,B),((((C,F),D),E)P,G)),H);", format=1),
                Tree("(((A,B),((((C,F),E),D)P,G)),H);", format=1),
                Tree("(((A,B),(((C,F),(E,D))P,G)),H);", format=1),
                Tree("(((A,B),((((D,E),C),F)P,G)),H);", format=1),
                Tree("(((A,B),((((D,E),F),C)P,G)),H);", format=1),
                Tree("(((A,B),((((D,F),C),E)P,G)),H);", format=1),
                Tree("(((A,B),((((D,F),E),C)P,G)),H);", format=1),
                Tree("(((A,B),((((E,F),C),D)P,G)),H);", format=1),
                Tree("(((A,B),((((E,F),D),C)P,G)),H);", format=1),
            ),
        )
        self.assertTreeResults(
            binarize(Tree("(((A,B),(C,D,E,F)P,G),H);", format=1)),
            (
                Tree(newick, quoted_node_names=True)
                for newick in Tree(
                    "(((A,B),(C,D,E,F)P,G),H);", format=1
                ).expand_polytomies()
            ),
        )
