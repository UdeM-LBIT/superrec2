from itertools import permutations
import unittest
from .toposort import toposort, toposort_all


class TestUtilsToposort(unittest.TestCase):
    def test_toposort(self):
        graph = {
            0: [],
            1: [],
            2: [3],
            3: [1],
            4: [0, 1],
            5: [0, 2],
        }

        sort = toposort(graph)
        self.assertIn(sort, toposort_all(graph))

        for node, node_succs in graph.items():
            for node_succ in node_succs:
                self.assertTrue(
                    sort.index(node) < sort.index(node_succ),
                    f"{node} should come before {node_succ}",
                )

        self.assertEqual(
            toposort(
                {
                    0: [1],
                    1: [2],
                    2: [0],
                }
            ),
            None,
        )

        self.assertEqual(
            toposort(
                {
                    0: [0],
                    1: [],
                    2: [3],
                    3: [1],
                    4: [0, 1],
                    5: [0, 2],
                }
            ),
            None,
        )

    def test_toposort_all(self):
        self.assertEqual(
            toposort_all(
                {
                    0: [],
                    1: [],
                    2: [3],
                    3: [1],
                    4: [0, 1],
                    5: [0, 2],
                }
            ),
            [
                [4, 5, 0, 2, 3, 1],
                [4, 5, 2, 0, 3, 1],
                [4, 5, 2, 3, 0, 1],
                [4, 5, 2, 3, 1, 0],
                [5, 2, 3, 4, 0, 1],
                [5, 2, 3, 4, 1, 0],
                [5, 2, 4, 0, 3, 1],
                [5, 2, 4, 3, 0, 1],
                [5, 2, 4, 3, 1, 0],
                [5, 4, 0, 2, 3, 1],
                [5, 4, 2, 0, 3, 1],
                [5, 4, 2, 3, 0, 1],
                [5, 4, 2, 3, 1, 0],
            ],
        )

        self.assertEqual(
            toposort_all({i: [] for i in range(6)}),
            [list(perm) for perm in permutations(range(6))],
        )

        self.assertEqual(
            toposort_all(
                {
                    0: [1],
                    1: [2],
                    2: [0],
                }
            ),
            [],
        )

        self.assertEqual(
            toposort_all(
                {
                    0: [0],
                    1: [],
                    2: [3],
                    3: [1],
                    4: [0, 1],
                    5: [0, 2],
                }
            ),
            [],
        )
