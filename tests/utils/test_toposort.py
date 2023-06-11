from itertools import permutations
from typing import List, TypeVar
from superrec2.utils.toposort import toposort, toposort_all, find_cycle


T = TypeVar("T")


def _assert_rotation_of(l1: List[T], l2: List[T]):
    if len(l1) != len(l2):
        raise AssertionError(
            f"'{l1}' is not a rotation of '{l2}': lengths \
differ"
        )

    for i in range(len(l1)):
        if l1[i:] + l1[:i] == l2:
            return

    for i in range(len(l1)):
        if l1[i::-1] + l1[:i:-1] == l2:
            return

    raise AssertionError(f"'{l1}' is not a rotation of '{l2}'")


def test_toposort():
    graph = {
        0: [],
        1: [],
        2: [3],
        3: [1],
        4: [0, 1],
        5: [0, 2],
    }

    sort = toposort(graph)
    assert sort in toposort_all(graph)

    for node, node_succs in graph.items():
        for node_succ in node_succs:
            assert sort.index(node) < sort.index(node_succ)

    assert (
        toposort(
            {
                0: [1],
                1: [2],
                2: [0],
            }
        )
        is None
    )

    assert (
        toposort(
            {
                0: [0],
                1: [],
                2: [3],
                3: [1],
                4: [0, 1],
                5: [0, 2],
            }
        )
        is None
    )


def test_toposort_all():
    assert toposort_all(
        {
            0: [],
            1: [],
            2: [3],
            3: [1],
            4: [0, 1],
            5: [0, 2],
        }
    ) == [
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
    ]

    assert toposort_all({i: [] for i in range(6)}) == [
        list(perm) for perm in permutations(range(6))
    ]

    assert (
        toposort_all(
            {
                0: [1],
                1: [2],
                2: [0],
            }
        )
        == []
    )

    assert (
        toposort_all(
            {
                0: [0],
                1: [],
                2: [3],
                3: [1],
                4: [0, 1],
                5: [0, 2],
            }
        )
        == []
    )


def test_find_cycle():
    _assert_rotation_of(
        find_cycle(
            {
                0: [1],
                1: [2],
                2: [0],
            }
        ),
        [0, 1, 2],
    )

    assert find_cycle(
        {
            0: [0],
            1: [],
            2: [3],
            3: [1],
            4: [0, 1],
            5: [0, 2],
        }
    ) == [0]

    _assert_rotation_of(
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

    assert (
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
        is None
    )
