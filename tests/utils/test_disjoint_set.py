import unittest
from itertools import product
from superrec2.utils.disjoint_set import DisjointSet


class TestDisjointSet(unittest.TestCase):
    def test_union_find(self):
        uf = DisjointSet(13)
        self.assertEqual(len(uf), 13)

        for i, j in product(range(13), repeat=2):
            if i != j:
                self.assertNotEqual(uf.find(i), uf.find(j))

        self.assertTrue(uf.unite(4, 6))
        self.assertEqual(len(uf), 12)
        self.assertEqual(uf.find(4), uf.find(6))

        self.assertTrue(uf.unite(0, 2))
        self.assertEqual(len(uf), 11)
        self.assertTrue(uf.unite(2, 4))
        self.assertEqual(len(uf), 10)
        self.assertTrue(uf.unite(6, 8))
        self.assertEqual(len(uf), 9)
        self.assertTrue(uf.unite(8, 10))
        self.assertEqual(len(uf), 8)
        self.assertTrue(uf.unite(10, 12))
        self.assertEqual(len(uf), 7)

        self.assertTrue(uf.unite(1, 3))
        self.assertEqual(len(uf), 6)
        self.assertTrue(uf.unite(3, 5))
        self.assertEqual(len(uf), 5)
        self.assertTrue(uf.unite(5, 7))
        self.assertEqual(len(uf), 4)
        self.assertTrue(uf.unite(7, 9))
        self.assertEqual(len(uf), 3)
        self.assertTrue(uf.unite(9, 11))
        self.assertEqual(len(uf), 2)

        self.assertTrue(all(uf.find(i) == uf.find(0) for i in range(0, 13, 2)))
        self.assertTrue(all(uf.find(i) == uf.find(1) for i in range(1, 13, 2)))

        self.assertFalse(uf.unite(1, 5))
        self.assertFalse(uf.unite(1, 1))
        self.assertTrue(uf.unite(0, 1))
        self.assertEqual(len(uf), 1)

    def test_repr(self):
        uf = DisjointSet(4)
        self.assertEqual(repr(uf), "DisjointSet({{0}, {1}, {2}, {3}})")
        uf.unite(0, 3)
        self.assertEqual(repr(uf), "DisjointSet({{0, 3}, {1}, {2}})")
        uf.unite(1, 0)
        self.assertEqual(repr(uf), "DisjointSet({{0, 1, 3}, {2}})")
        uf.unite(2, 0)
        self.assertEqual(repr(uf), "DisjointSet({{0, 1, 2, 3}})")

    def test_binary(self):
        uf = DisjointSet(4)
        uf_bin = uf.binary()

        self.assertCountEqual(
            [repr(element) for element in uf_bin],
            (
                "DisjointSet({{0}, {1, 2, 3}})",
                "DisjointSet({{0, 2, 3}, {1}})",
                "DisjointSet({{0, 1, 3}, {2}})",
                "DisjointSet({{0, 1, 2}, {3}})",
                "DisjointSet({{0, 1}, {2, 3}})",
                "DisjointSet({{0, 2}, {1, 3}})",
                "DisjointSet({{0, 3}, {1, 2}})",
            ),
        )
