import random
import time
from superrec2.utils.range_min_query import RangeMinQuery


def _compute_all_naive(data):
    size = len(data)
    results = [[None] * size for _ in range(size)]

    for i in range(size):
        for j in range(i + 1, size):
            results[i][j] = min(data[i:j])

    return results


def _compute_all_rmq(data):
    size = len(data)
    results = [[None] * size for _ in range(size)]
    rmq = RangeMinQuery(data)

    for i in range(len(data)):
        for j in range(i + 1, len(data)):
            results[i][j] = rmq(i, j)

    return results


def test_valid():
    data = [3, 1, 5, 3, 4, 7, 6, 1]
    rmq = RangeMinQuery(data)

    for i in range(len(data)):
        for j in range(len(data)):
            if i < j:
                assert rmq(i, j) == min(data[i:j])
            else:
                assert rmq(i, j) is None


def test_change_source():
    data = [8, 8, 8]
    rmq = RangeMinQuery(data)
    data[0] = 0

    for i in range(len(data)):
        for j in range(i + 1, len(data)):
            assert rmq(i, j) == 8


def test_faster():
    size = 1_000
    data = random.sample(range(100_000), size)

    naive_start = time.time()
    naive_results = _compute_all_naive(data)
    naive_dur = time.time() - naive_start

    rmq_start = time.time()
    rmq_results = _compute_all_rmq(data)
    rmq_dur = time.time() - rmq_start

    for i in range(len(data)):
        for j in range(i + 1, len(data)):
            assert naive_results[i][j] == rmq_results[i][j]

    print("Naive time:", naive_dur)
    print("RMQ time:", rmq_dur)
    assert rmq_dur < naive_dur
