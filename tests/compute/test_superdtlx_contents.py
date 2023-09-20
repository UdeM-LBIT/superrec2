from sowing.indexed import IndexedTree
from superrec2.compute.superdtlx.contents import compute_min_contents
from superrec2.model.history import Associate, parse_tree


def test_compute_min_contents():
    associate_tree = parse_tree(
        Associate,
        """
        (
            (
                (4[&contents='{"a","b","c"}'],5[&contents='{"a","b"}'])3,
                6[&contents='{"c"}']
            )2,
            (
                8[&contents='{"e","c"}'],
                (10[&contents='{"e","b","c"}'],11[&contents='{"d","c"}'])9
            )7
        )1;
        """,
    )
    nodes = IndexedTree(associate_tree)
    assert compute_min_contents(associate_tree) == {
        nodes["1"]: frozenset("bc"),
        nodes["2"]: frozenset("bc"),
        nodes["3"]: frozenset("abc"),
        nodes["4"]: frozenset("abc"),
        nodes["5"]: frozenset("ab"),
        nodes["6"]: frozenset("c"),
        nodes["7"]: frozenset("bce"),
        nodes["8"]: frozenset("ec"),
        nodes["9"]: frozenset("bce"),
        nodes["10"]: frozenset("ebc"),
        nodes["11"]: frozenset("dc"),
    }

    associate_tree = parse_tree(
        Associate,
        """
        (
            (
                (
                    (a[&contents='{"a","b"}'],b[&contents='{"a","c"}'])1,
                    c[&contents='{"a","d"}']
                )2,
                d[&contents='{"b","d"}']
            )3,
            (e[&contents='{"c","d","e"}'],f[&contents='{"e","f"}'])4
        )5;
        """,
    )
    nodes = IndexedTree(associate_tree)
    assert compute_min_contents(associate_tree) == {
        nodes["a"]: frozenset("ab"),
        nodes["b"]: frozenset("ac"),
        nodes["c"]: frozenset("ad"),
        nodes["d"]: frozenset("bd"),
        nodes["e"]: frozenset("cde"),
        nodes["f"]: frozenset("ef"),
        nodes["1"]: frozenset("abc"),
        nodes["2"]: frozenset("abcd"),
        nodes["3"]: frozenset("bcd"),
        nodes["4"]: frozenset("cde"),
        nodes["5"]: frozenset("cd"),
    }
