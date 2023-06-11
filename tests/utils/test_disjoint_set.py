from itertools import product
from superrec2.utils.disjoint_set import DisjointSet


def test_union_find():
    uf = DisjointSet(13)
    assert len(uf) == 13

    for i, j in product(range(13), repeat=2):
        if i != j:
            assert uf.find(i) != uf.find(j)

    assert uf.unite(4, 6)
    assert len(uf) == 12
    assert uf.find(4) == uf.find(6)

    assert uf.unite(0, 2)
    assert len(uf) == 11
    assert uf.unite(2, 4)
    assert len(uf) == 10
    assert uf.unite(6, 8)
    assert len(uf) == 9
    assert uf.unite(8, 10)
    assert len(uf) == 8
    assert uf.unite(10, 12)
    assert len(uf) == 7

    assert uf.unite(1, 3)
    assert len(uf) == 6
    assert uf.unite(3, 5)
    assert len(uf) == 5
    assert uf.unite(5, 7)
    assert len(uf) == 4
    assert uf.unite(7, 9)
    assert len(uf) == 3
    assert uf.unite(9, 11)
    assert len(uf) == 2

    assert all(uf.find(i) == uf.find(0) for i in range(0, 13, 2))
    assert all(uf.find(i) == uf.find(1) for i in range(1, 13, 2))

    assert not uf.unite(1, 5)
    assert not uf.unite(1, 1)
    assert uf.unite(0, 1)
    assert len(uf) == 1


def test_repr():
    uf = DisjointSet(4)
    assert repr(uf) == "DisjointSet({{0}, {1}, {2}, {3}})"
    uf.unite(0, 3)
    assert repr(uf) == "DisjointSet({{0, 3}, {1}, {2}})"
    uf.unite(1, 0)
    assert repr(uf) == "DisjointSet({{0, 1, 3}, {2}})"
    uf.unite(2, 0)
    assert repr(uf) == "DisjointSet({{0, 1, 2, 3}})"


def test_binary():
    uf = DisjointSet(4)
    uf_bin = uf.binary()

    assert {repr(element) for element in uf_bin} == {
        "DisjointSet({{0}, {1, 2, 3}})",
        "DisjointSet({{0, 2, 3}, {1}})",
        "DisjointSet({{0, 1, 3}, {2}})",
        "DisjointSet({{0, 1, 2}, {3}})",
        "DisjointSet({{0, 1}, {2, 3}})",
        "DisjointSet({{0, 2}, {1, 3}})",
        "DisjointSet({{0, 3}, {1, 2}})",
    }
