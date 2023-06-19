from sowing.node import Node as N
from sowing import traversal
from sowphy.clade import Clade as C
from sowphy import newick
from superrec2.model.history import (
    InvalidReconciliation,
    Associate as Assoc,
    Reconciliation,
    Event,
    Extant,
    Codiverge,
    Diverge,
    Transfer,
    Gain,
    Loss,
    History,
)
from immutables import Map
import pytest


def test_reconciliation_valid():
    rec_empty = Reconciliation(host_tree=None, associate_tree=None)
    rec_empty.validate()
    assert rec_empty.is_complete()
    assert rec_empty.erase() == rec_empty

    rec_full = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("4"), clade=C("5")))
            .add(N(Assoc(host=C("3"), clade=C("3")))
                 .add(N(Assoc(host=C("8"), clade=C("1"))))
                 .add(N(Assoc(host=C("9"), clade=C("2"))))
            )
            .add(N(Assoc(host=C("9"), clade=C("4"))))
        ),
    )
    rec_full.validate()
    assert rec_full.is_complete()

    rec_part = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(C("5"))
            .add(N(C("3"))
                 .add(N(Assoc(host=C("8"), clade=C("1"))))
                 .add(N(Assoc(host=C("9"), clade=C("2"))))
            )
            .add(N(Assoc(host=C("9"), clade=C("4"))))
        ),
    )
    rec_part.validate()
    assert not rec_part.is_complete()

    assert rec_full.erase() == rec_part


def test_reconciliation_invalid():
    with pytest.raises(InvalidReconciliation) as err:
        rec = Reconciliation(
            host_tree=newick.parse("((4,5,10)3,(6,(8,9)7)2)1;"),
            associate_tree=(
                N(C("5"))
                .add(N(C("3"))
                    .add(N(Assoc(host=C("8"), clade=C("1"))))
                    .add(N(Assoc(host=C("9"), clade=C("2"))))
                )
                .add(N(Assoc(host=C("9"), clade=C("4"))))
            ),
        )
        rec.validate()

    assert "host tree must be binary" in str(err.value)
    assert err.value.node is None

    with pytest.raises(InvalidReconciliation) as err:
        rec = Reconciliation(
            host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
            associate_tree=(
                N(C("6"))
                .add(N(C("3"))
                    .add(N(Assoc(host=C("8"), clade=C("1"))))
                    .add(N(Assoc(host=C("9"), clade=C("2"))))
                )
                .add(N(Assoc(host=C("9"), clade=C("4"))))
                .add(N(Assoc(host=C("4"), clade=C("5"))))
            ),
        )
        rec.validate()

    assert "associate tree must be binary" in str(err.value)
    assert err.value.node is None

    with pytest.raises(InvalidReconciliation) as err:
        rec = Reconciliation(
            host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
            associate_tree=(
                N(C("6"))
                .add(N(C("3"))
                    .add(N(Assoc(host=C("8"), clade=C("1"))))
                    .add(N(C("2")))
                )
                .add(N(Assoc(host=C("9"), clade=C("4"))))
            ),
        )
        rec.validate()

    assert "leaf associates must be labeled by Associate instances" in str(err.value)
    assert err.value.node.data == C("2")

    with pytest.raises(InvalidReconciliation) as err:
        rec = Reconciliation(
            host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
            associate_tree=(
                N(C("6"))
                .add(N(C("3"))
                    .add(N(Assoc(host=C("8"), clade=C("1"))))
                    .add(N(Assoc(host=C("7"), clade=C("2"))))
                )
                .add(N(Assoc(host=C("9"), clade=C("4"))))
            ),
        )
        rec.validate()

    assert "leaf associate host Clade(name='7') is not terminal" in str(err.value)
    assert err.value.node.data == Assoc(host=C("7"), clade=C("2"))


def test_events():
    with pytest.raises(TypeError) as err:
        Event(Assoc(host=C()))

    assert "abstract class Event" in str(err.value)
    assert Codiverge(Assoc(host=C())).arity == 2
    assert Diverge(Assoc(host=C()), result=2, segment=(0, 2), cut=False).arity == 2
    assert Transfer(Assoc(host=C())).arity == 1
    assert Gain(Assoc(host=C()), index=1, gain=tuple("abc")).arity == 1
    assert Loss(Assoc(host=C()), segment=(0, 0)).arity == 0
    assert Loss(Assoc(host=C(), contents=tuple("abc")), segment=(1, 2)).arity == 1
    assert Loss(Assoc(host=C(), contents=tuple("abc")), segment=(0, 3)).arity == 0


def test_history_extant():
    # Valid extant associate in terminal host
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=N(Assoc(host=C("8"), clade=C("1"), contents=tuple("ab"))),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=N(Extant(Assoc(host=C("8"), clade=C("1"), contents=tuple("ab")))),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Valid extant associate in unsampled host
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8[&&NHX:sampled=],9)7)2)1;"),
        associate_tree=None,
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8[&&NHX:sampled=],9)7)2)1;"),
        event_tree=N(Extant(Assoc(
            host=C("8", props=Map({"sampled": ""})),
            clade=C("1"),
            contents=tuple("ab"),
        ))),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Invalid extant: Non-terminal host
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=N(Extant(Assoc(host=C("7"), clade=C("1"), contents=tuple("ab")))),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert "leaf associate host Clade(name='7') is not terminal" in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid extant: Non-leaf node
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Extant(Assoc(host=C("8"), clade=C("1"), contents=tuple("ab"))))
            .add(N(Extant(Assoc(host=C("8"), clade=C("2"), contents=tuple("ab")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert "Extant node must be a leaf, found 1 child(ren) instead" in str(err.value)
    assert err.value.node == hist.event_tree


def test_history_codiverge():
    # Valid codivergence matching host order
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc")))
            .add(N(Assoc(host=C("4"), clade=C("2"), contents=tuple("abc"))))
            .add(N(Assoc(host=C("5"), clade=C("3"), contents=tuple("abc"))))
        ),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Codiverge(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc"))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("abc")))))
            .add(N(Extant(Assoc(host=C("5"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Valid codivergence in reverse host order
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc")))
            .add(N(Assoc(host=C("5"), clade=C("2"), contents=tuple("abc"))))
            .add(N(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc"))))
        ),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Codiverge(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc"))))
            .add(N(Extant(Assoc(host=C("5"), clade=C("2"), contents=tuple("abc")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Invalid codivergence: Mismatched child associate hosts
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Codiverge(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc"))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("abc")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "codivergence children are linked to (Clade(name='4'), Clade(name='4'))"
        " instead of their host’s children (Clade(name='4'), Clade(name='5'))"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid codivergence: Mismatched left child contents
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Codiverge(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc"))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("ab")))))
            .add(N(Extant(Assoc(host=C("5"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "codivergence left child contents ('a', 'b') does not equal"
        " its parent’s contents ('a', 'b', 'c')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid codivergence: Mismatched right child contents
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Codiverge(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc"))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("abc")))))
            .add(N(Extant(Assoc(host=C("5"), clade=C("3"), contents=tuple("ab")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "codivergence right child contents ('a', 'b') does not equal"
        " its parent’s contents ('a', 'b', 'c')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid codivergence: Invalid arity
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Codiverge(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc"))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("abc")))))
            .add(N(Extant(Assoc(host=C("5"), clade=C("3"), contents=tuple("abc")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("4"), contents=tuple("abc")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert "Codiverge node must be binary, found 3 child(ren) instead" in str(err.value)
    assert err.value.node == hist.event_tree


def test_history_diverge():
    # Valid full duplication
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")))
            .add(N(Assoc(host=C("4"), clade=C("2"), contents=tuple("abc"))))
            .add(N(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc"))))
        ),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=0, segment=(0, 3), cut=False,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("abc")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Valid partial duplication with result in left child
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")))
            .add(N(Assoc(host=C("4"), clade=C("2"), contents=tuple("a"))))
            .add(N(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc"))))
        ),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=0, segment=(0, 1), cut=False,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("a")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Valid partial duplication with result in right child
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")))
            .add(N(Assoc(host=C("4"), clade=C("2"), contents=tuple("abc"))))
            .add(N(Assoc(host=C("4"), clade=C("3"), contents=tuple("a"))))
        ),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=1, segment=(0, 1), cut=False,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("abc")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("a")))))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Valid cut with result in left child
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")))
            .add(N(Assoc(host=C("4"), clade=C("2"), contents=tuple("a"))))
            .add(N(Assoc(host=C("4"), clade=C("3"), contents=tuple("bc"))))
        ),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=0, segment=(0, 1), cut=True,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("a")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("bc")))))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Valid cut with result in right child
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")))
            .add(N(Assoc(host=C("4"), clade=C("2"), contents=tuple("bc"))))
            .add(N(Assoc(host=C("4"), clade=C("3"), contents=tuple("a"))))
        ),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=1, segment=(0, 1), cut=True,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("bc")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("a")))))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Invalid divergence: Out of range result index
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=2, segment=(0, 1), cut=True,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("a")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert "result index is 2, expected 0 or 1" in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid divergence: Result child in different host
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=0, segment=(0, 1), cut=True,
            ))
            .add(N(Extant(Assoc(host=C("5"), clade=C("2"), contents=tuple("a")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "divergence result child host Clade(name='5')"
        " differs from its parent host Clade(name='4')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid divergence: Conserved child in different host
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=0, segment=(0, 1), cut=True,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("a")))))
            .add(N(Extant(Assoc(host=C("5"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "divergence conserved child host Clade(name='5')"
        " differs from its parent host Clade(name='4')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid divergence: Unequal copied contents
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=0, segment=(0, 1), cut=False,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("ab")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "divergence result contents ('a', 'b') differ from"
        " the targeted segment ('a',)"
    ) in str(err.value)
    assert err.value.node == hist.event_tree


    # Invalid duplication: Unequal conserved contents
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=0, segment=(0, 1), cut=False,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("a")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("ab")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "copy-divergence conserved child contents ('a', 'b') differ from"
        " its parent’s contents ('a', 'b', 'c')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid cut: Mismatched conserved contents
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=0, segment=(0, 1), cut=True,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("a")))))
            .add(N(Extant(Assoc(host=C("4"), clade=C("3"), contents=tuple("abc")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "cut-divergence conserved child contents ('a', 'b', 'c') differ from"
        " the remaining contents ('b', 'c')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid divergence: Invalid arity
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Diverge(
                Assoc(host=C("4"), clade=C("1"), contents=tuple("abc")),
                result=0, segment=(0, 1), cut=True,
            ))
            .add(N(Extant(Assoc(host=C("4"), clade=C("2"), contents=tuple("a")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert "Diverge node must be binary, found 1 child(ren) instead" in str(err.value)
    assert err.value.node == hist.event_tree


def test_history_transfer():
    # Valid transfer
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("6"), clade=C("2"), contents=tuple("abc")))
        ),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Transfer(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc"))))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("abc")))))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Invalid transfer: Targets descendant
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Transfer(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc"))))
            .add(N(Extant(Assoc(host=C("5"), clade=C("2"), contents=tuple("abc")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "transfer child host Clade(name='5') is comparable to its origin host"
        " Clade(name='3')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid transfer: Targets ancestor
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Transfer(Assoc(host=C("5"), clade=C("1"), contents=tuple("abc")))).add(
                N(Transfer(Assoc(host=C("3"), clade=C("1"), contents=tuple("abc"))))
                .add(N(Extant(Assoc(host=C("8"), clade=C("2"), contents=tuple("abc")))))
            )
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "transfer child host Clade(name='3') is comparable to its origin host"
        " Clade(name='5')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid transfer: Differing contents
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Transfer(Assoc(host=C("5"), clade=C("1"), contents=tuple("abc"))))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("ab")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "transfer child contents ('a', 'b') differs from its"
        " parent’s contents ('a', 'b', 'c')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid transfer: Invalid arity
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Transfer(Assoc(host=C("5"), clade=C("1"), contents=tuple("abc"))))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("ab")))))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("ab")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert "Transfer node must be unary, found 2 child(ren) instead" in str(err.value)
    assert err.value.node == hist.event_tree


def test_history_gain():
    # Valid gain
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("6"), clade=C("2"), contents=tuple("adefb")))
        ),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Gain(
                Assoc(host=C("6"), clade=C("1"), contents=tuple("ab")),
                index=1, gain=tuple("def"),
            ))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("adefb")))))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Invalid gain: Different host
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Gain(
                Assoc(host=C("6"), clade=C("1"), contents=tuple("ab")),
                index=1, gain=tuple("def"),
            ))
            .add(N(Extant(Assoc(host=C("8"), clade=C("2"), contents=tuple("adefb")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "gain child host Clade(name='8') differs from its"
        " parent host Clade(name='6')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid gain: Mismatched contents
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Gain(
                Assoc(host=C("6"), clade=C("1"), contents=tuple("ab")),
                index=1, gain=tuple("def"),
            ))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("adeb")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "gain child contents ('a', 'd', 'e', 'b') differ from the expected"
        " gain result ('a', 'd', 'e', 'f', 'b')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid gain: Invalid arity
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Gain(
                Assoc(host=C("6"), clade=C("1"), contents=tuple("ab")),
                index=1, gain=tuple("def"),
            ))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("adefb")))))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("adefb")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert "Gain node must be unary, found 2 child(ren) instead" in str(err.value)


def test_history_loss():
    # Valid partial loss
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=(
            N(Assoc(host=C("6"), clade=C("2"), contents=tuple("ad")))
        ),
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Loss(
                Assoc(host=C("6"), clade=C("1"), contents=tuple("abcd")),
                segment=(1, 3),
            ))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("ad")))))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Valid complete loss
    rec = Reconciliation(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        associate_tree=None,
    )

    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Loss(
                Assoc(host=C("6"), clade=C("1"), contents=tuple("abcd")),
                segment=(0, 4),
            ))
        ),
    )

    rec.validate()
    hist.validate()
    assert hist.compress() == rec

    # Invalid loss: Different hosts
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Loss(
                Assoc(host=C("6"), clade=C("1"), contents=tuple("abcd")),
                segment=(1, 3),
            ))
            .add(N(Extant(Assoc(host=C("8"), clade=C("2"), contents=tuple("ad")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "loss child host Clade(name='8') differs from its"
        " parent host Clade(name='6')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid loss: Mismatched contents
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Loss(
                Assoc(host=C("6"), clade=C("1"), contents=tuple("abcd")),
                segment=(1, 3),
            ))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("abd")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert (
        "loss child contents ('a', 'b', 'd') differ from the expected"
        " loss result ('a', 'd')"
    ) in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid loss: Invalid partial arity
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Loss(
                Assoc(host=C("6"), clade=C("1"), contents=tuple("abcd")),
                segment=(1, 3),
            ))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert "Loss node must be unary, found 0 child(ren) instead" in str(err.value)
    assert err.value.node == hist.event_tree

    # Invalid loss: Invalid complete arity
    hist = History(
        host_tree=newick.parse("((4,5)3,(6,(8,9)7)2)1;"),
        event_tree=(
            N(Loss(
                Assoc(host=C("6"), clade=C("1"), contents=tuple("abcd")),
                segment=(0, 4),
            ))
            .add(N(Extant(Assoc(host=C("6"), clade=C("2"), contents=tuple("abd")))))
        ),
    )

    with pytest.raises(InvalidReconciliation) as err:
        hist.validate()

    assert "Loss node must be a leaf, found 1 child(ren) instead" in str(err.value)
    assert err.value.node == hist.event_tree
