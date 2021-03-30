import unittest
from itertools import product
from .disjoint_set import UnionFind


class TestDisjointSet(unittest.TestCase):
    def test_valid(self):
        uf = UnionFind(13)

        for i, j in product(range(13), repeat=2):
            if i != j:
                self.assertNotEqual(uf.find(i), uf.find(j))

        uf.unite(4, 6)
        self.assertEqual(uf.find(4), uf.find(6))

        uf.unite(0, 2)
        uf.unite(2, 4)
        uf.unite(6, 8)
        uf.unite(8, 10)
        uf.unite(10, 12)

        uf.unite(1, 3)
        uf.unite(3, 5)
        uf.unite(5, 7)
        uf.unite(7, 9)
        uf.unite(9, 11)

        self.assertTrue(all(uf.find(i) == uf.find(0) for i in range(0, 13, 2)))

        self.assertTrue(all(uf.find(i) == uf.find(1) for i in range(1, 13, 2)))
