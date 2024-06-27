from sowing.indexed import IndexedTree
from superrec2.compute.synesth.contents import compute_min_contents, propagate_contents
from superrec2.model.history import Host, Associate, parse_tree, Event, History


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


def test_propagate_contents():
    host_tree = parse_tree(Host, "(D,(B,C)A)E;")

    hist_flat = History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              2[&host=B,contents='{"a","__extra__"}'],
              3[&host=C,contents='{"a","b"}']
            )
            1[&kind=codiverge,host=A,contents='{"a","b"}'];
            """,
        ),
    )
    prop_flat = propagate_contents(hist_flat)
    prop_flat.validate()

    assert prop_flat == History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              2[&host=B,contents='{"a","b"}'],
              3[&host=C,contents='{"a","b"}']
            )
            1[&kind=codiverge,host=A,contents='{"a","b"}'];
            """,
        ),
    )

    hist_nested = History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              (
                4[&host=B,contents='{"__extra__"}'],
                5[&host=C,contents='{"b","__extra__"}']
              )2[&kind=codiverge,host=A,contents='{"a","__extra__"}'],
              3[&host=D,contents='{"a","b","c"}']
            )
            1[&kind=codiverge,host=E,contents='{"a","b","c"}'];
            """,
        ),
    )
    prop_nested = propagate_contents(hist_nested)
    prop_nested.validate()

    assert prop_nested == History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              (
                4[&host=B,contents='{"a","b","c"}'],
                5[&host=C,contents='{"a","b","c"}']
              )2[&kind=codiverge,host=A,contents='{"a","b","c"}'],
              3[&host=D,contents='{"a","b","c"}']
            )
            1[&kind=codiverge,host=E,contents='{"a","b","c"}'];
            """,
        ),
    )

    hist_loss_extra = History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              (
                2[&host=B,contents='{"a"}']
              )[&kind=loss,host=B,contents='{"a","__extra__"}',segment='{"__extra__"}'],
              3[&host=C,contents='{"__extra__"}']
            )
            1[&kind=codiverge,host=A,contents='{"a","b","c"}'];
            """,
        ),
    )
    prop_loss_extra = propagate_contents(hist_loss_extra)
    prop_loss_extra.validate()

    assert prop_loss_extra == History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              (
                2[&host=B,contents='{"a"}']
              )[&kind=loss,host=B,contents='{"a","b","c"}',segment='{"b","c"}'],
              3[&host=C,contents='{"a","b","c"}']
            )
            1[&kind=codiverge,host=A,contents='{"a","b","c"}'];
            """,
        ),
    )

    hist_loss_other = History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              2[&host=B,contents='{"a","b","c"}'],
              (
                3[&host=C,contents='{"__extra__"}']
              )[&kind=loss,host=C,contents='{"b","__extra__"}',segment='{"b"}'],
            )
            1[&kind=codiverge,host=A,contents='{"a","b","c"}'];
            """,
        ),
    )
    prop_loss_other = propagate_contents(hist_loss_other)
    prop_loss_other.validate()

    assert prop_loss_other == History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              2[&host=B,contents='{"a","b","c"}'],
              (
                3[&host=C,contents='{"a","c"}']
              )[&kind=loss,host=C,contents='{"a","b","c"}',segment='{"b"}']
            )
            1[&kind=codiverge,host=A,contents='{"a","b","c"}'];
            """,
        ),
    )

    hist_dup = History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              2[&host=C,contents='{"a","__extra__"}'],
              3[&host=C,contents='{"a"}']
            )
            1[&kind=diverge,host=C,contents='{"a","b"}',segment={"a"},result=1];
            """,
        ),
    )
    prop_dup = propagate_contents(hist_dup)
    prop_dup.validate()

    assert prop_dup == History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              2[&host=C,contents='{"a","b"}'],
              3[&host=C,contents='{"a"}']
            )
            1[&kind=diverge,host=C,contents='{"a","b"}',segment={"a"},result=1];
            """,
        ),
    )

    hist_cut = History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              2[&host=C,contents='{"a"}'],
              3[&host=C,contents='{"b","__extra__"}']
            )
            1[&kind=diverge,host=C,contents='{"a","b","c"}',segment={"a"},cut=True];
            """,
        ),
    )
    prop_cut = propagate_contents(hist_cut)
    prop_cut.validate()

    assert prop_cut == History(
        host_tree=host_tree,
        event_tree=parse_tree(
            Event,
            """
            (
              2[&host=C,contents='{"a"}'],
              3[&host=C,contents='{"b","c"}']
            )
            1[&kind=diverge,host=C,contents='{"a","b","c"}',segment={"a"},cut=True];
            """,
        ),
    )
