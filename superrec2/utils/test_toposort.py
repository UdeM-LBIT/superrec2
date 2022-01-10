from itertools import permutations
from typing import List, TypeVar
import unittest
from .toposort import toposort, toposort_all, find_cycle


T = TypeVar("T")


def assertRotationOf(l1: List[T], l2: List[T]):
    if len(l1) != len(l2):
        raise AssertionError(f"'{l1}' is not a rotation of '{l2}': lengths \
differ")

    for i in range(len(l1)):
        if l1[i:] + l1[:i] == l2:
            return

    for i in range(len(l1)):
        if l1[i::-1] + l1[:i:-1] == l2:
            return

    raise AssertionError(f"'{l1}' is not a rotation of '{l2}'")


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

    def test_find_cycle(self):
        assertRotationOf(
            find_cycle(
                {
                    0: [1],
                    1: [2],
                    2: [0],
                }
            ),
            [0, 1, 2],
        )

        self.assertEqual(
            find_cycle(
                {
                    0: [0],
                    1: [],
                    2: [3],
                    3: [1],
                    4: [0, 1],
                    5: [0, 2],
                }
            ),
            [0],
        )

        assertRotationOf(
            find_cycle(
                {
                    0: [1],
                    1: [2],
                    2: [3, 4],
                    3: [4, 5],
                    4: [6],
                    5: [4],
                    6: [1],
                }
            ),
            [1, 2, 4, 6],
        )

        self.assertIsNone(
            find_cycle(
                {
                    0: [],
                    1: [],
                    2: [3],
                    3: [1],
                    4: [0, 1],
                    5: [0, 2],
                }
            )
        )
