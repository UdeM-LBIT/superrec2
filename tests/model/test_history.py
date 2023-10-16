from sowing.indexed import IndexedTree
from superrec2.model.history import (
    parse_tree,
    graft_unsampled_hosts,
    Host,
    Associate as Assoc,
    InvalidReconciliation,
    Reconciliation,
    InvalidEvent,
    Event,
    Extant,
    Codiverge,
    Diverge,
    Gain,
    Loss,
    History,
)
from superrec2.model.graph import CycleError
import pytest


def test_associate():
    empty_assoc = Assoc()
    assert not empty_assoc.is_complete()
    assert repr(empty_assoc) == "Associate()"
    assert empty_assoc == Assoc.from_mapping({})

    ord_assoc = Assoc(name="a", host="1", contents=tuple("abc"))
    assert ord_assoc.is_complete()
    assert repr(ord_assoc) == "Associate(name='a', host='1', contents=('a', 'b', 'c'))"
    assert ord_assoc == Assoc.from_mapping(
        {
            "name": "a",
            "host": "1",
            "contents": "('a', 'b', 'c')",
        }
    )

    assert ord_assoc.gain((0, tuple())) == ord_assoc
    assert ord_assoc.gain((2, tuple("de"))) == Assoc(
        name="a",
        host="1",
        contents=tuple("abdec"),
    )

    assert ord_assoc.split((1, 2)) == (
        Assoc(name="a", host="1", contents=tuple("b")),
        Assoc(name="a", host="1", contents=tuple("ac")),
    )
    assert ord_assoc.split((0, 3)) == (
        Assoc(name="a", host="1", contents=tuple("abc")),
        Assoc(name="a", host="1", contents=tuple()),
    )

    assert ord_assoc.switch("2") == Assoc(name="a", host="2", contents=tuple("abc"))

    und_assoc = Assoc(name="b", host="2", contents=frozenset("abc"))
    assert und_assoc.is_complete()
    assert repr(und_assoc) == "Associate(name='b', host='2', contents={'a', 'b', 'c'})"
    assert und_assoc == Assoc.from_mapping(
        {
            "name": "b",
            "host": "2",
            "contents": "{'a', 'b', 'c'}",
        }
    )

    assert und_assoc.gain(frozenset()) == und_assoc
    assert und_assoc.gain(frozenset("de")) == Assoc(
        name="b",
        host="2",
        contents=frozenset("abcde"),
    )

    assert und_assoc.split(frozenset()) == (
        Assoc(name="b", host="2", contents=frozenset()),
        Assoc(name="b", host="2", contents=frozenset("abc")),
    )
    assert und_assoc.split(frozenset("bc")) == (
        Assoc(name="b", host="2", contents=frozenset("bc")),
        Assoc(name="b", host="2", contents=frozenset("a")),
    )

    with pytest.raises(ValueError) as err:
        und_assoc.split(frozenset("d"))

    assert (
        "split argument {'d'} is not a subset of existing contents {'a', 'b', 'c'}"
        in str(err.value)
    )


def test_host():
    host_tree = parse_tree(Host, "((4,5)3,(6,(8,9)7)2)1;")
    grafted_tree = parse_tree(
        Host,
        """
        (
          '1[U]'[&sampled=false],
          (
            (
              '3[U]'[&sampled=false],
              (
                ('4[U]'[&sampled=false],4)'4[P]',
                ('5[U]'[&sampled=false],5)'5[P]'
              )3
            )'3[P]',
            (
              '2[U]'[&sampled=false],
              (
                ('6[U]'[&sampled=false],6)'6[P]',
                (
                  '7[U]'[&sampled=false],
                  (
                    ('8[U]'[&sampled=false],8)'8[P]',
                    ('9[U]'[&sampled=false],9)'9[P]'
                  )7
                )'7[P]'
              )2
            )'2[P]'
          )1
        )'1[P]';
        """,
    )
    assert graft_unsampled_hosts(host_tree) == grafted_tree


def test_reconciliation_valid():
    rec_empty = Reconciliation(host_tree=None, associate_tree=None)
    rec_empty.validate()
    assert rec_empty.is_complete()
    assert rec_empty.erase() == rec_empty

    rec_full = Reconciliation(
        host_tree=parse_tree(Host, "((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=parse_tree(
            Assoc,
            "((1[&host=8],2[&host=9])3[&host=3],4[&host=9])5[&host=4];",
        ),
    )
    rec_full.validate()
    assert rec_full.is_complete()

    rec_part = Reconciliation(
        host_tree=parse_tree(Host, "((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=parse_tree(Assoc, "((1[&host=8],2[&host=9])3,4[&host=9])5;"),
    )
    rec_part.validate()
    assert not rec_part.is_complete()

    assert rec_full.erase() == rec_part


def test_reconciliation_invalid():
    with pytest.raises(InvalidReconciliation) as err:
        rec = Reconciliation(
            host_tree=parse_tree(Host, "((4,5,10)3,(6,(8,9)7)2)1;"),
            associate_tree=parse_tree(Assoc, "((1[&host=8],2[&host=9])3,4[&host=9])5;"),
        )
        rec.validate()

    assert "host tree must be binary" in str(err.value)
    assert err.value.node is None

    with pytest.raises(InvalidReconciliation) as err:
        rec = Reconciliation(
            host_tree=parse_tree(Host, "((4,5)3,(6,(8,9)7)2)1;"),
            associate_tree=parse_tree(
                Assoc,
                "((1[&host=8],2[&host=9])3,4[&host=9],5[&host=4])6;",
            ),
        )
        rec.validate()

    assert "associate tree must be binary" in str(err.value)
    assert err.value.node is None

    with pytest.raises(InvalidReconciliation) as err:
        rec = Reconciliation(
            host_tree=parse_tree(Host, "((4,5)3,(6,(8,9)7)2)1;"),
            associate_tree=parse_tree(Assoc, "((1[&host=8],2)3,4[&host=9])6;"),
        )
        rec.validate()

    assert "leaf associates must have complete information" in str(err.value)
    assert err.value.node.data == Assoc("2")

    with pytest.raises(InvalidReconciliation) as err:
        rec = Reconciliation(
            host_tree=parse_tree(Host, "((4,5)3,(6,(8,9)7)2)1;"),
            associate_tree=parse_tree(Assoc, "((1[&host=8],2[&host=7])3,4[&host=9])6;"),
        )
        rec.validate()

    assert "leaf associate host '7' is not terminal" in str(err.value)
    assert err.value.node.data == Assoc(host="7", name="2")


def test_event():
    with pytest.raises(TypeError) as err:
        Event()

    assert "abstract class Event" in str(err.value)


def test_event_extant():
    host_index1 = IndexedTree(parse_tree(Host, "(2,3)1;"))
    host_index2 = IndexedTree(parse_tree(Host, "(1,2)3;"))

    ordered = Extant(name="x", host="3", contents=tuple("abc"))
    assert ordered == Extant.from_mapping(
        {
            "name": "x",
            "host": "3",
            "contents": "('a', 'b', 'c')",
        }
    )
    assert ordered.arity == 0
    assert ordered.associate() == Assoc(name="x", host="3", contents=tuple("abc"))
    assert ordered.anon_associate() == Assoc(host="3", contents=tuple("abc"))
    ordered.validate(host_index1, ())

    unordered = Extant(name="x", host="3", contents=frozenset("abc"))
    assert unordered == Extant.from_mapping(
        {"name": "x", "host": "3", "contents": "{'a', 'b', 'c'}"}
    )
    assert unordered.arity == 0
    assert unordered.associate() == Assoc(name="x", host="3", contents=frozenset("abc"))
    assert unordered.anon_associate() == Assoc(host="3", contents=frozenset("abc"))
    unordered.validate(host_index1, ())

    with pytest.raises(InvalidEvent) as err:
        ordered.validate(host_index2, ())

    assert "extant host '3' is not terminal" in str(err.value)
    assert err.value.event == ordered

    with pytest.raises(InvalidEvent) as err:
        unordered.validate(host_index2, ())

    assert "extant host '3' is not terminal" in str(err.value)
    assert err.value.event == unordered


def test_event_codiverge():
    host_index = IndexedTree(parse_tree(Host, "(2,3)1;"))

    spe = Codiverge(name="x", host="1", contents=tuple("abc"))
    assert spe == Codiverge.from_mapping(
        {
            "name": "x",
            "host": "1",
            "contents": "('a', 'b', 'c')",
        }
    )
    assert spe.arity == 2
    spe.validate(
        host_index,
        (
            Assoc(host="2", contents=tuple("abc")),
            Assoc(host="3", contents=tuple("abc")),
        ),
    )
    spe.validate(
        host_index,
        (
            Assoc(host="3", contents=tuple("abc")),
            Assoc(host="2", contents=tuple("abc")),
        ),
    )

    with pytest.raises(InvalidEvent) as err:
        spe.validate(
            host_index,
            (
                Assoc(host="2", contents=tuple("abc")),
                Assoc(host="2", contents=tuple("abc")),
            ),
        )

    assert (
        "codivergence event children are "
        "(Associate(host='2', contents=('a', 'b', 'c')), "
        "Associate(host='2', contents=('a', 'b', 'c'))), expected "
        "(Associate(host='2', contents=('a', 'b', 'c')), "
        "Associate(host='3', contents=('a', 'b', 'c')))"
    ) in str(err.value)

    with pytest.raises(InvalidEvent) as err:
        spe.validate(
            host_index,
            (
                Assoc(host="2", contents=tuple("abc")),
                Assoc(host="3", contents=tuple("ab")),
            ),
        )

    assert (
        "codivergence event children are "
        "(Associate(host='2', contents=('a', 'b', 'c')), "
        "Associate(host='3', contents=('a', 'b'))), expected "
        "(Associate(host='2', contents=('a', 'b', 'c')), "
        "Associate(host='3', contents=('a', 'b', 'c')))"
    ) in str(err.value)


def test_event_duplicate():
    host_index = IndexedTree(parse_tree(Host, "(2,3)1;"))

    dup = Diverge(
        name="x",
        host="3",
        contents=tuple("abc"),
        result=1,
        segment=(0, 2),
        cut=False,
    )
    assert dup == Diverge.from_mapping(
        {
            "name": "x",
            "host": "3",
            "contents": "('a', 'b', 'c')",
            "result": "1",
            "segment": "(0, 2)",
            "cut": "false",
        }
    )
    assert dup.arity == 2

    dup.validate(
        host_index,
        (
            Assoc(host="3", contents=tuple("abc")),
            Assoc(host="3", contents=tuple("ab")),
        ),
    )

    with pytest.raises(InvalidEvent) as err:
        dup.validate(
            host_index,
            (
                Assoc(host="3", contents=tuple("abc")),
                Assoc(host="3", contents=tuple("abc")),
            ),
        )

    assert (
        "divergence result child is "
        "Associate(host='3', contents=('a', 'b', 'c')), "
        "expected Associate(host='3', contents=('a', 'b'))"
    ) in str(err.value)

    with pytest.raises(InvalidEvent) as err:
        dup.validate(
            host_index,
            (
                Assoc(host="3", contents=tuple("ab")),
                Assoc(host="3", contents=tuple("ab")),
            ),
        )

    assert (
        "divergence conserved child is "
        "Associate(host='3', contents=('a', 'b')), "
        "expected Associate(host='3', contents=('a', 'b', 'c'))"
    ) in str(err.value)

    with pytest.raises(InvalidEvent) as err:
        dup.validate(
            host_index,
            (
                Assoc(host="3", contents=tuple("abc")),
                Assoc(host="2", contents=tuple("ab")),
            ),
        )

    assert "copy-divergence result host '2' differs from its parent host '3'" in str(
        err.value
    )

    cut = Diverge(
        name="x",
        host="3",
        contents=tuple("abc"),
        result=1,
        segment=(0, 2),
        cut=True,
    )
    assert cut == Diverge.from_mapping(
        {
            "name": "x",
            "host": "3",
            "contents": "('a', 'b', 'c')",
            "result": "1",
            "segment": "(0, 2)",
            "cut": "true",
        }
    )
    assert cut.arity == 2

    cut.validate(
        host_index,
        (
            Assoc(host="3", contents=tuple("c")),
            Assoc(host="3", contents=tuple("ab")),
        ),
    )

    ucut = Diverge(
        host="3",
        contents=frozenset("abc"),
        result=1,
        segment=frozenset("ab"),
        cut=True,
    )
    assert ucut == Diverge.from_mapping(
        {
            "host": "3",
            "contents": "{'a', 'b', 'c'}",
            "result": "1",
            "segment": "{'a', 'b'}",
            "cut": "true",
        }
    )
    assert ucut.arity == 2

    ucut.validate(
        host_index,
        (
            Assoc(host="3", contents=frozenset("c")),
            Assoc(host="3", contents=frozenset("ab")),
        ),
    )

    fcut = Diverge(
        name="x",
        host="3",
        contents=tuple("abc"),
        result=0,
        segment=(0, 3),
        cut=True,
    )
    assert fcut == Diverge.from_mapping(
        {
            "name": "x",
            "host": "3",
            "contents": "('a', 'b', 'c')",
            "result": "0",
            "segment": "(0, 3)",
            "cut": "true",
        }
    )
    assert fcut.arity == 1
    fcut.validate(host_index, (Assoc(host="3", contents=tuple("abc")),))

    inv = Diverge(
        name="x",
        host="3",
        contents=tuple("abc"),
        result=2,
        segment=(0, 2),
        cut=False,
    )

    with pytest.raises(InvalidEvent) as err:
        inv.validate(
            host_index,
            (
                Assoc(host="3", contents=tuple("abc")),
                Assoc(host="3", contents=tuple("abc")),
            ),
        )

    assert (
        "divergence result index 2 is out of bounds (should be in range [0..1])"
        in str(err.value)
    )

    inv = Diverge(
        name="x",
        host="3",
        contents=tuple("abc"),
        result=1,
        segment=(0, 3),
        cut=True,
    )

    with pytest.raises(InvalidEvent) as err:
        inv.validate(host_index, (Assoc(host="3", contents=tuple("abc")),))

    assert (
        "divergence result index 1 is out of bounds (should be in range [0..0])"
        in str(err.value)
    )


def test_event_transfer():
    host_index = IndexedTree(parse_tree(Host, "(2,3)1;"))

    tra = Diverge(
        name="x",
        host="3",
        contents=tuple("abc"),
        result=1,
        segment=(0, 2),
        cut=False,
        transfer=True,
    )
    assert tra == Diverge.from_mapping(
        {
            "name": "x",
            "host": "3",
            "contents": "('a', 'b', 'c')",
            "result": "1",
            "segment": "(0, 2)",
            "cut": "false",
            "transfer": "true",
        }
    )
    assert tra.arity == 2

    tra.validate(
        host_index,
        (
            Assoc(host="3", contents=tuple("abc")),
            Assoc(host="2", contents=tuple("ab")),
        ),
    )

    with pytest.raises(InvalidEvent) as err:
        tra.validate(
            host_index,
            (
                Assoc(host="3", contents=tuple("abc")),
                Assoc(host="1", contents=tuple("ab")),
            ),
        )

    assert (
        "transfer-divergence target host '1' is comparable " "to its origin host '3'"
    ) in str(err.value)

    ftra = Diverge(
        name="x",
        host="3",
        contents=tuple("abc"),
        result=0,
        segment=(0, 3),
        cut=True,
        transfer=True,
    )
    assert ftra == Diverge.from_mapping(
        {
            "name": "x",
            "host": "3",
            "contents": "('a', 'b', 'c')",
            "result": "0",
            "segment": "(0, 3)",
            "cut": "true",
            "transfer": "true",
        }
    )
    assert ftra.arity == 1

    ftra.validate(host_index, (Assoc(host="2", contents=tuple("abc")),))


def test_event_gain():
    host_index = IndexedTree(parse_tree(Host, "(2,3)1;"))

    gain = Gain(name="x", host="1", contents=tuple("ab"), gained=(1, tuple("bc")))
    assert gain == Gain.from_mapping(
        {
            "name": "x",
            "host": "1",
            "contents": "('a', 'b')",
            "gained": "(1, ('b', 'c'))",
        }
    )
    assert gain.arity == 1
    gain.validate(host_index, (Assoc(host="1", contents=tuple("abcb")),))

    with pytest.raises(InvalidEvent) as err:
        gain.validate(host_index, (Assoc(host="2", contents=tuple("abcb")),))

    assert (
        "gain event child is Associate(host='2', contents=('a', 'b', 'c', 'b')), "
        "expected Associate(host='1', contents=('a', 'b', 'c', 'b'))"
    ) in str(err.value)

    with pytest.raises(InvalidEvent) as err:
        gain.validate(host_index, (Assoc(host="1", contents=tuple("acbb")),))

    assert (
        "gain event child is Associate(host='1', contents=('a', 'c', 'b', 'b')), "
        "expected Associate(host='1', contents=('a', 'b', 'c', 'b'))"
    ) in str(err.value)

    ugain = Gain(name="x", host="1", contents=frozenset("ab"), gained=frozenset("c"))
    assert ugain == Gain.from_mapping(
        {
            "name": "x",
            "host": "1",
            "contents": "{'a', 'b'}",
            "gained": "{'c'}",
        }
    )
    assert ugain.arity == 1
    ugain.validate(host_index, (Assoc(host="1", contents=frozenset("abc")),))

    with pytest.raises(InvalidEvent) as err:
        ugain.validate(host_index, (Assoc(host="2", contents=frozenset("abc")),))

    assert (
        "gain event child is Associate(host='2', contents={'a', 'b', 'c'}), "
        "expected Associate(host='1', contents={'a', 'b', 'c'})"
    ) in str(err.value)

    with pytest.raises(InvalidEvent) as err:
        ugain.validate(host_index, (Assoc(host="1", contents=frozenset("abcd")),))

    assert (
        "gain event child is Associate(host='1', contents={'a', 'b', 'c', 'd'}), "
        "expected Associate(host='1', contents={'a', 'b', 'c'})"
    ) in str(err.value)


def test_event_loss():
    host_index = IndexedTree(parse_tree(Host, "(2,3)1;"))

    loss = Loss(name="x", host="1", contents=tuple("abc"), segment=(1, 2))
    assert loss == Loss.from_mapping(
        {
            "name": "x",
            "host": "1",
            "contents": "('a', 'b', 'c')",
            "segment": "(1, 2)",
        }
    )
    assert loss.arity == 1
    loss.validate(host_index, (Assoc(host="1", contents=tuple("ac")),))

    with pytest.raises(InvalidEvent) as err:
        loss.validate(host_index, (Assoc(host="2", contents=tuple("ac")),))

    assert (
        "loss event child is Associate(host='2', contents=('a', 'c')), "
        "expected Associate(host='1', contents=('a', 'c'))"
    ) in str(err.value)

    with pytest.raises(InvalidEvent) as err:
        loss.validate(host_index, (Assoc(host="1", contents=tuple("a")),))

    assert (
        "loss event child is Associate(host='1', contents=('a',)), "
        "expected Associate(host='1', contents=('a', 'c'))"
    ) in str(err.value)

    floss = Loss(name="x", host="1", contents=tuple("abc"), segment=(0, 3))
    assert floss == Loss.from_mapping(
        {
            "name": "x",
            "host": "1",
            "contents": "('a', 'b', 'c')",
            "segment": "(0, 3)",
        }
    )
    assert floss.arity == 0
    floss.validate(host_index, ())

    uloss = Loss(name="x", host="1", contents=frozenset("abc"), segment=frozenset("ab"))
    assert uloss == Loss.from_mapping(
        {
            "name": "x",
            "host": "1",
            "contents": "{'a', 'b', 'c'}",
            "segment": "{'a', 'b'}",
        }
    )
    assert uloss.arity == 1
    uloss.validate(host_index, (Assoc(host="1", contents=frozenset("c")),))


def test_history_validate():
    valid = History(
        host_tree=parse_tree(Host, "(1,2)3;"),
        event_tree=parse_tree(
            Event,
            '(x[&host=1,contents=\'{"a","b"}\'],y[&host=2,contents=\'{"a","b"}\'])'
            '[&kind=codiverge,host=3,contents=\'{"a","b"}\'];',
        ),
    )
    valid.validate()

    unknown_host = History(
        host_tree=parse_tree(Host, "(a,2)3;"),
        event_tree=parse_tree(
            Event,
            '(x[&host=1,contents=\'{"a","b"}\'],y[&host=2,contents=\'{"a","b"}\'])'
            '[&kind=codiverge,host=3,contents=\'{"a","b"}\'];',
        ),
    )

    with pytest.raises(InvalidEvent) as err:
        unknown_host.validate()

    assert "undefined event host '1'" in str(err.value)

    invalid_arity = History(
        host_tree=parse_tree(Host, "(1,2)3;"),
        event_tree=parse_tree(
            Event,
            '(x[&host=1,contents=\'{"a","b"}\'])'
            '[&kind=codiverge,host=3,contents=\'{"a","b"}\'];',
        ),
    )

    with pytest.raises(InvalidEvent) as err:
        invalid_arity.validate()

    assert "codiverge event must be binary, found 1 child(ren) instead" in str(
        err.value
    )

    event_error = History(
        host_tree=parse_tree(Host, "(1,2)3;"),
        event_tree=parse_tree(
            Event,
            '(x[&host=1,contents=\'{"a","b"}\'],y[&host=1,contents=\'{"a","b"}\'])'
            '[&kind=codiverge,host=3,contents=\'{"a","b"}\'];',
        ),
    )

    with pytest.raises(InvalidEvent) as err:
        event_error.validate()

    assert (
        "codivergence event children are "
        "(Associate(host='1', contents={'a', 'b'}), "
        "Associate(host='1', contents={'a', 'b'})), expected "
        "(Associate(host='1', contents={'a', 'b'}), "
        "Associate(host='2', contents={'a', 'b'}))"
    ) in str(err.value)


def test_history_compress():
    compress_unary = History(
        host_tree=parse_tree(Host, "((1,2)3,4)5;"),
        event_tree=parse_tree(
            Event,
            """(
                ((x[&host=1,contents=\'{"a","c"}\'])
                    [&kind=loss,host=1,contents=\'{"a","b","c"}\',segment=\'{"b"}\'])
                    [&kind=gain,host=1,contents=\'{"a","b"}\',gained=\'{"c"}\'],
                ((y[&host=4,contents=\'{"b"}\'])
                    [&kind=loss,host=4,contents=\'{"a","b"}\',segment=\'{"a"}\'])
                    [&kind=diverge,host=2,contents=\'{"a","b"}\',segment=\'{"a","b"}\',
                        transfer=true,cut=true]
            )[&kind=codiverge,host=3,contents=\'{"a","b"}\'];""",
        ),
    )

    compress_unary.validate()
    assert compress_unary.compress() == Reconciliation(
        host_tree=parse_tree(Host, "((1,2)3,4)5;"),
        associate_tree=parse_tree(
            Assoc,
            '(x[&host=1,contents=\'{"a","c"}\'],y[&host=4,contents=\'{"b"}\'])'
            '[&host=3,contents=\'{"a","b"}\'];',
        ),
    )

    compress_loss = History(
        host_tree=parse_tree(Host, "(1,2)3;"),
        event_tree=parse_tree(
            Event,
            """(
                x[&host=1,contents=\'{"a","b"}\'],
                [&kind=loss,host=2,contents=\'{"a","b"}\',segment=\'{"a","b"}\']
            )[&kind=codiverge,host=3,contents=\'{"a","b"}\'];""",
        ),
    )

    compress_loss.validate()
    assert compress_loss.compress() == Reconciliation(
        host_tree=parse_tree(Host, "(1,2)3;"),
        associate_tree=parse_tree(
            Assoc,
            'x[&host=1,contents=\'{"a","b"}\'];',
        ),
    )

    compress_unsampled = History(
        host_tree=parse_tree(Host, "((1,(2,2U[&sampled=false])2P)3,4)5;"),
        event_tree=parse_tree(
            Event,
            """(
                a[&host=1,contents=\'{"a","b"}\'],
                (
                    b[&host=2,contents=\'{"a","b"}\'],
                    (
                        c[&host=4,contents=\'{"b"}\'],
                        [&host=2U,contents=\'{"a","b"}\']
                    )
                    [&
                        kind=diverge,
                        result=0,
                        host=2U,
                        contents=\'{"a","b"}\',
                        transfer=true,
                        segment=\'{"b"}\'
                    ]
                )[&kind=codiverge,host=2P,contents=\'{"a","b"}\']
            )[&kind=codiverge,host=3,contents=\'{"a","b"}\'];""",
        ),
    )

    compress_unsampled.validate()
    assert compress_unsampled.compress() == Reconciliation(
        host_tree=parse_tree(Host, "((1,(2,2U[&sampled=false])2P)3,4)5;"),
        associate_tree=parse_tree(
            Assoc,
            """(
                a[&host=1,contents=\'{"a","b"}\'],
                (
                    b[&host=2,contents=\'{"a","b"}\'],
                    c[&host=4,contents=\'{"b"}\']
                )[&host=2P,contents=\'{"a","b"}\']
            )[&host=3,contents=\'{"a","b"}\'];""",
        ),
    )


def test_history_epochs():
    simple = History(
        host_tree=parse_tree(Host, "(1,(2,(3,4)t)s)r;"),
        event_tree=parse_tree(
            Event,
            """
            (
              x[&host=1,contents='{"a","b"}'],
              (
                [&kind=loss,host=t,contents='{"a","b"}',segment='{"a","b"}'],
                y[&host=2,contents='{"a","b"}']
              )[&kind=codiverge,host=s,contents='{"a","b"}']
            )[&kind=codiverge,host=r,contents='{"a","b"}'];
            """,
        ),
    )
    simple.validate()
    assert simple.epochs() == {
        "1": 3,
        "2": 3,
        "3": 3,
        "4": 3,
        "t": 2,
        "s": 1,
        "r": 0,
    }

    with_transfer = History(
        host_tree=parse_tree(Host, "((1,2)s,(3,4)t)r;"),
        event_tree=parse_tree(
            Event,
            """
            (
              (
                [&kind=loss,host=s,contents='{"a","b"}',segment='{"a","b"}'],
                (
                  [&kind=loss,host=3,contents='{"a","b"}',segment='{"a","b"}'],
                  [&host=2,contents='{"a","b"}']
                )[&kind=diverge,host=3,contents='{"a","b"}',segment='{"a","b"}',result=1,transfer=true]
              )[&kind=diverge,host=s,contents='{"a","b"}',segment='{"a","b"}',result=1,transfer=true],
              (
                [&kind=loss,host=3,contents='{"a","b"}',segment='{"a","b"}'],
                [&host=4,contents='{"a","b"}']
              )[&kind=codiverge,host=t,contents='{"a","b"}']
            )[&kind=codiverge,host=r,contents='{"a","b"}'];
            """,
        ),
    )
    with_transfer.validate()
    assert with_transfer.epochs() == {
        "1": 3,
        "2": 3,
        "3": 3,
        "4": 3,
        "s": 2,
        "t": 1,
        "r": 0,
    }

    with_cycle = History(
        host_tree=parse_tree(Host, "((1,2)s,(3,4)t)r;"),
        event_tree=parse_tree(
            Event,
            """
            (
              [&kind=loss,host=t,contents='{"a","b"}',segment='{"a","b"}'],
              (
                [&kind=loss,host=1,contents='{"a","b"}',segment='{"a","b"}'],
                (
                  [&kind=loss,host=2,contents='{"a","b"}',segment='{"a","b"}'],
                  (
                    [&kind=loss,host=4,contents='{"a","b"}',segment='{"a","b"}'],
                    (
                      [&kind=loss,host=3,contents='{"a","b"}',segment='{"a","b"}'],
                      (
                        [&kind=loss,host=2,contents='{"a","b"}',segment='{"a","b"}'],
                        [&host=1,contents='{"a","b"}']
                      )[&kind=codiverge,host=s,contents='{"a","b"}']
                    )[&
                      kind=diverge,
                      host=3,
                      contents='{"a","b"}',
                      segment='{"a","b"}',
                      result=1,
                      transfer=true
                    ]
                  )[&kind=codiverge,host=t,contents='{"a","b"}']
                )[&
                  kind=diverge,
                  host=2,
                  contents='{"a","b"}',
                  segment='{"a","b"}',
                  result=1,
                  transfer=true
                ]
              )[&kind=codiverge,host=s,contents='{"a","b"}']
            )[&kind=codiverge,host=r,contents='{"a","b"}'];
            """,
        ),
    )
    with_cycle.validate()

    with pytest.raises(CycleError, match="negative-weight cycle exists") as err:
        with_cycle.epochs()

    assert err.value.args[1] == ["s", "t"]
