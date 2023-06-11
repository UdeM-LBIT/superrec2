import pytest
from infinity import inf
from superrec2.utils.dynamic_programming import (
    Table,
    TableProxy,
    Candidate,
    Entry,
    EntryProxy,
    DictDimension,
    ListDimension,
    MergePolicy,
    RetentionPolicy,
)


def test_bound_checking():
    table = Table([ListDimension(10), ListDimension(20)])

    assert type(table[2]) == TableProxy
    assert type(table[2][8]) == EntryProxy
    assert table[2][8] == table.entry(inf, [])

    with pytest.raises(TypeError):
        table[2][8][0]

    table[2][8].update(Candidate(5))
    assert table[2][8] == table.entry(5, [])

    table[2][9] = Candidate(6)
    assert table[2][9] == table.entry(6, [])

    with pytest.raises(TypeError):
        table[2] = Candidate(5)

    with pytest.raises(TypeError):
        table[2][8][0] = Candidate(5)


def test_dimensions():
    table_l = Table([ListDimension(10)])
    table_l[8] = Candidate(0)
    assert list(table_l.keys()) == list(range(10))
    assert type(table_l[7]) == EntryProxy
    assert table_l[8] == table_l.entry(0, [])

    i = 0

    for j in table_l:
        assert i == j
        i += 1

    assert i == 10

    table_lll = Table([ListDimension(10), ListDimension(20), ListDimension(30)])
    table_lll[3][2][1] = Candidate(0)
    assert list(table_lll.keys()) == list(range(10))
    assert list(table_lll[1].keys()) == list(range(20))
    assert list(table_lll[1][2].keys()) == list(range(30))
    assert type(table_lll[1][2][3]) == EntryProxy
    assert table_lll[3][2][1] == table_lll.entry(0, [])

    i = 0

    for j in table_lll[3][2]:
        assert i == j
        i += 1

    assert i == 30

    table_d = Table([DictDimension()])
    key = object()
    assert list(table_d.keys()) == []
    table_d[key] = Candidate(0)
    assert list(table_d.keys()) == [key]
    assert type(table_d[object()]) == EntryProxy
    assert table_d[key] == table_d.entry(0, [])

    table_ddd = Table([DictDimension(), DictDimension(), DictDimension()])
    x = object()
    y1 = object()
    y2 = object()
    z = object()
    assert list(table_ddd.keys()) == []
    table_ddd[x][y1][z] = Candidate(0)
    assert list(table_ddd.keys()) == [x]
    assert list(table_ddd[x].keys()) == [y1]
    assert list(table_ddd[x][y1].keys()) == [z]
    assert type(table_ddd[x][y2][z]) == EntryProxy
    assert table_ddd[x][y1][z] == table_ddd.entry(0, [])


def test_min_all():
    table = Table(
        [ListDimension(10), ListDimension(20)],
        MergePolicy.MIN,
        RetentionPolicy.ALL,
    )

    assert type(table[4][18]) == EntryProxy

    table[4][18] = Candidate(4)
    assert table[4][18] == table.entry(4, [])

    table[4][18] = Candidate(3)
    assert table[4][18] == table.entry(3, [])

    table[4][18] = Candidate(3)
    assert table[4][18] == table.entry(3, [])

    table[4][18] = Candidate(3, "tag1")
    assert table[4][18] == table.entry(3, ["tag1"])
    assert table[4][18].info() == "tag1"

    table[4][18] = Candidate(3, "tag2")
    assert table[4][18] == table.entry(3, ["tag1", "tag2"])
    assert table[4][18].info() in ["tag1", "tag2"]

    table[4][18] = Candidate(3, "tagA")
    table[4][18] = Candidate(3, "tagB")
    table[4][18] = Candidate(3, "tag0")
    assert table[4][18] == table.entry(3, ["tag1", "tag2", "tagA", "tagB", "tag0"])

    table[4][18] = Candidate(3, "tag4")
    table[4][18] = Candidate(2, "tagα")
    table[4][18] = Candidate(3, "tag5")
    assert table[4][18] == table.entry(2, ["tagα"])


def test_max_all():
    table = Table(
        [ListDimension(10), ListDimension(20)],
        MergePolicy.MAX,
        RetentionPolicy.ALL,
    )

    assert type(table[4][18]) == EntryProxy

    table[4][18] = Candidate(4)
    assert table[4][18] == table.entry(4, [])

    table[4][18] = Candidate(5)
    assert table[4][18] == table.entry(5, [])

    table[4][18] = Candidate(5)
    assert table[4][18] == table.entry(5, [])

    table[4][18] = Candidate(5, "tag1")
    assert table[4][18] == table.entry(5, ["tag1"])
    assert table[4][18].info() == "tag1"

    table[4][18] = Candidate(5, "tag2")
    assert table[4][18] == table.entry(5, ["tag1", "tag2"])
    assert table[4][18].info() in ["tag1", "tag2"]

    table[4][18] = Candidate(5, "tagA")
    table[4][18] = Candidate(5, "tagB")
    table[4][18] = Candidate(5, "tag0")

    assert table[4][18] == table.entry(5, ["tag1", "tag2", "tagA", "tagB", "tag0"])

    table[4][18] = Candidate(5, "tag4")
    table[4][18] = Candidate(6, "tagα")
    table[4][18] = Candidate(5, "tag5")
    assert table[4][18] == table.entry(6, ["tagα"])


def test_max_any():
    table = Table(
        [ListDimension(10), ListDimension(20)],
        MergePolicy.MAX,
        RetentionPolicy.ANY,
    )

    assert type(table[4][18]) == EntryProxy

    table[4][18] = Candidate(4)
    assert table[4][18] == table.entry(4, [])

    table[4][18] = Candidate(5)
    assert table[4][18] == table.entry(5, [])

    table[4][18] = Candidate(5)
    assert table[4][18] == table.entry(5, [])

    table[4][18] = Candidate(5, "tag1")
    assert table[4][18] == table.entry(5, ["tag1"])
    assert table[4][18].info() == "tag1"

    table[4][18] = Candidate(5, "tag2")
    assert table[4][18].value() == 5
    assert len(table[4][18].infos()) == 1
    assert table[4][18].info() in ["tag1", "tag2"]

    table[4][18] = Candidate(5, "tagA")
    table[4][18] = Candidate(5, "tagB")
    table[4][18] = Candidate(5, "tag0")
    assert table[4][18].value() == 5
    assert len(table[4][18].infos()) == 1
    assert table[4][18].info() in ["tag1", "tag2", "tagA", "tagB", "tag0"]

    table[4][18] = Candidate(5, "tag4")
    table[4][18] = Candidate(6, "tagα")
    table[4][18] = Candidate(5, "tag5")
    assert table[4][18] == table.entry(infos=["tagα"], value=6)


def test_standalone_entries():
    entry_1 = Entry(MergePolicy.MIN, RetentionPolicy.ALL)
    entry_2 = Entry(MergePolicy.MIN, RetentionPolicy.ALL)

    entry_1.update(
        Candidate(5, "tag1"),
        Candidate(5, "tag2"),
        Candidate(6, "tag3"),
        Candidate(5, "tag4"),
    )

    entry_2.update(
        Candidate(2, "tag1"),
        Candidate(3, "tag2"),
        Candidate(3, "tag3"),
        Candidate(2, "tag4"),
    )

    entry = entry_1.combine(
        entry_2,
        lambda x, y: Candidate(x.value + y.value, (x.info, y.info)),
    )
    assert entry == Entry(
        7,
        [
            ("tag1", "tag1"),
            ("tag1", "tag4"),
            ("tag2", "tag1"),
            ("tag2", "tag4"),
            ("tag4", "tag1"),
            ("tag4", "tag4"),
        ],
        MergePolicy.MIN,
        RetentionPolicy.ALL,
    )

    assert set(entry_2) == {
        Candidate(2, "tag1"),
        Candidate(2, "tag4"),
    }

    entry_1.update(*entry_2)
    assert entry_1.value() == 2
    assert entry_1.infos() == {"tag1", "tag4"}

    entry_1 = Entry(MergePolicy.MIN, RetentionPolicy.ANY)
    entry_2 = Entry(MergePolicy.MIN, RetentionPolicy.ANY)

    entry_1.update(
        Candidate(5, "tag1"),
        Candidate(5, "tag2"),
        Candidate(6, "tag3"),
        Candidate(5, "tag4"),
    )

    entry_2.update(
        Candidate(2, "tag1"),
        Candidate(3, "tag2"),
        Candidate(3, "tag3"),
        Candidate(2, "tag4"),
    )

    entry = entry_1.combine(
        entry_2,
        lambda x, y: Candidate(x.value + y.value, (x.info, y.info)),
    )
    assert entry.value() == 7
    assert len(entry.infos()) == 1
    assert entry.info() in [
        ("tag1", "tag1"),
        ("tag1", "tag4"),
        ("tag2", "tag1"),
        ("tag2", "tag4"),
        ("tag4", "tag1"),
        ("tag4", "tag4"),
    ]


def _edit_distance(w1, w2, retention):
    table = Table(
        [
            ListDimension(len(w1) + 1),
            ListDimension(len(w2) + 1),
        ],
        MergePolicy.MIN,
        retention,
    )

    table[0][0] = Candidate(0)

    for i in range(1, len(w1) + 1):
        table[i][0] = Candidate(i, (i - 1, 0))

    for j in range(1, len(w2) + 1):
        table[0][j] = Candidate(j, (0, j - 1))

    for i in range(1, len(w1) + 1):
        for j in range(1, len(w2) + 1):
            table[i][j] = Candidate(
                table[i - 1][j - 1].value() + (1 if w1[i - 1] != w2[j - 1] else 0),
                (i - 1, j - 1),
            )
            table[i][j] = Candidate(table[i][j - 1].value() + 1, (i, j - 1))
            table[i][j] = Candidate(table[i - 1][j].value() + 1, (i - 1, j))

    return table


@staticmethod
def _get_alignments(w1, w2, table, i=None, j=None):
    if i is None:
        i = len(w1)

    if j is None:
        j = len(w2)

    if i == 0 and j == 0:
        yield ("", "")

    if not table[i][j].infos():
        return

    for prec in table[i][j].infos():
        for align in _get_alignments(w1, w2, table, *prec):
            if prec == (i - 1, j - 1):
                yield (
                    align[0] + w1[i - 1],
                    align[1] + w2[j - 1],
                )
            elif prec == (i, j - 1):
                yield (
                    align[0] + "-",
                    align[1] + w2[j - 1],
                )
            else:
                yield (
                    align[0] + w1[i - 1],
                    align[1] + "-",
                )


def test_edit_distance():
    w1 = "elephant"
    w2 = "relevant"

    table_all = _edit_distance(w1, w2, RetentionPolicy.ALL)
    table_any = _edit_distance(w1, w2, RetentionPolicy.ANY)
    table_none = _edit_distance(w1, w2, RetentionPolicy.NONE)

    assert table_all[len(w1)][len(w2)].value() == 3
    assert table_any[len(w1)][len(w2)].value() == 3
    assert table_none[len(w1)][len(w2)].value() == 3

    expected = [
        [0, 1, 2, 3, 4, 5, 6, 7, 8],
        [1, 1, 1, 2, 3, 4, 5, 6, 7],
        [2, 2, 2, 1, 2, 3, 4, 5, 6],
        [3, 3, 2, 2, 1, 2, 3, 4, 5],
        [4, 4, 3, 3, 2, 2, 3, 4, 5],
        [5, 5, 4, 4, 3, 3, 3, 4, 5],
        [6, 6, 5, 5, 4, 4, 3, 4, 5],
        [7, 7, 6, 6, 5, 5, 4, 3, 4],
        [8, 8, 7, 7, 6, 6, 5, 4, 3],
    ]

    for i in range(9):
        for j in range(9):
            assert expected[i][j] == table_all[i][j].value()
            assert expected[i][j] == table_any[i][j].value()
            assert expected[i][j] == table_none[i][j].value()

    all_solutions = set(_get_alignments(w1, w2, table_all))
    any_solution = list(_get_alignments(w1, w2, table_any))
    no_solution = list(_get_alignments(w1, w2, table_none))

    assert all_solutions == {("-elephant", "rele-vant"), ("-elephant", "relev-ant")}
    assert len(any_solution) == 1
    assert any_solution[0] in all_solutions
    assert no_solution == []
