from superrec2.model.history import (
    parse_tree,
    graft_unsampled_hosts,
    Associate,
    Host,
    Reconciliation,
    Event,
    History,
)
from superrec2.compute.superdtlx.recurrence import reconcile
from superrec2.utils.algebras import make_selector
from superrec2.compute.util import (
    make_cost_algebra,
    history_builder,
    history_generator,
)


unit_cost = make_cost_algebra("unit_cost", costs={})
best_unit_cost = make_selector("best_unit_cost", unit_cost, history_generator)

scaled_cost = make_cost_algebra(
    "scaled_cost",
    costs={
        "speciation": 0,
        "loss": 1,
        "dup": 2,
        "cut": 2.5,
        "transfer-dup": 4,
        "transfer-cut": 4.5,
    },
)
best_scaled_cost = make_selector("best_scaled_cost", scaled_cost, history_generator)


def _build_event_tree(source):
    return history_builder(parse_tree(Event, source))


def _tree_set(*trees):
    return frozenset(map(_build_event_tree, trees))


def _history_match_input(setting, event_tree):
    return (
        setting
        == History(
            host_tree=setting.host_tree,
            event_tree=event_tree.value,
        )
        .compress()
        .erase()
    )


def test_reconcile_simple():
    host_tree = parse_tree(Host, "(a,b)c;")

    # Two leaves with same contents in divergent hosts: Single speciation
    associate_tree = parse_tree(
        Associate, "(1[&host=a,contents='{\"x\"}'],2[&host=b,contents='{\"x\"}']);"
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, unit_cost) == 0

    results = reconcile(setting, best_unit_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            1[&host=a,contents='{"x"}'],
            2[&host=b,contents='{"x"}']
          )[&kind=codiverge,host=c,contents='{"x"}']
        )[&kind=gain,host=c,contents='set()',gained='{"x"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)

    # Two leaves with subsumed contents in divergent hosts: Speciation with gain
    associate_tree = parse_tree(
        Associate, '(1[&host=a,contents=\'{"x","y"}\'],2[&host=b,contents=\'{"x"}\']);'
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, unit_cost) == 0

    results = reconcile(setting, best_unit_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            (1[&host=a,contents='{"x","y"}'])
            [&kind=gain,host=a,contents='{"x"}',gained='{"y"}'],
            2[&host=b,contents='{"x"}']
          )[&kind=codiverge,host=c,contents='{"x"}']
        )[&kind=gain,host=c,contents='set()',gained='{"x"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)

    # Two leaves with disjoint contents in divergent hosts: Empty speciation with gains
    associate_tree = parse_tree(
        Associate, "(1[&host=a,contents='{\"x\"}'],2[&host=b,contents='{\"y\"}']);"
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, unit_cost) == 0

    results = reconcile(setting, best_unit_cost).selected.value
    assert results == _tree_set(
        """
        (
          (1[&host=a,contents='{"x"}'])
          [&kind=gain,host=a,contents='set()',gained='{"x"}'],
          (2[&host=b,contents='{"y"}'])
          [&kind=gain,host=b,contents='set()',gained='{"y"}']
        )[&kind=codiverge,host=c,contents='set()'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)


def test_reconcile_extra_contents():
    host_tree = parse_tree(Host, "(((a,b)c,d)e,f)g;")

    # Three leaves with disjoint contents and intermediate loss:
    # Factorized loss in leaves
    associate_tree = parse_tree(
        Associate,
        """
        (
          (
            (
              1[&host=a,contents='{"x","y"}'],
              2[&host=b,contents='{"y","z"}']
            ),
            3[&host=d,contents='{"w","x","y","z"}']
          ),
          4[&host=f,contents='{"w","x","y","z"}']
        );
        """,
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, scaled_cost) == 2

    results = reconcile(setting, best_unit_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            (
              (
                (
                  1[&host=a,contents='{"x","y"}']
                )[&kind=loss,host=a,contents='{"x","y","__extra__"}',segment='{"__extra__"}'],
                (
                  2[&host=b,contents='{"y","z"}']
                )[&kind=loss,host=b,contents='{"y","z","__extra__"}',segment='{"__extra__"}']
              )[&kind=codiverge,host=c,contents='{"x","y","z","__extra__"}'],
              3[&host=d,contents='{"w","x","y","z"}']
            )[&kind=codiverge,host=e,contents='{"w","x","y","z"}'],
            4[&host=f,contents='{"w","x","y","z"}']
          )[&kind=codiverge,host=g,contents='{"w","x","y","z"}']
        )[&kind=gain,host=g,contents='set()',gained='{"w","x","y","z"}'];
        """
    )


def test_reconcile_dup_cut():
    host_tree = parse_tree(Host, "(a,b)c;")

    # Three leaves in same host with subsumed contents: Partial duplications
    associate_tree = parse_tree(
        Associate,
        """
        (
          (1[&host=a,contents='{"x","y"}'],2[&host=a,contents='{"x","y","z"}']),
          3[&host=a,contents='{"x","y","z"}']
        );
        """,
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, scaled_cost) == 4

    results = reconcile(setting, best_scaled_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            (
              1[&host=a,contents='{"x","y"}'],
              2[&host=a,contents='{"x","y","z"}']
            )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y"}'],
            3[&host=a,contents='{"x","y","z"}']
          )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y","z"}']
        )[&kind=gain,host=a,contents='set()',gained='{"x","y","z"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)

    # Three leaves in same host with distinct and non-subsumed contents:
    # Partial duplication and partial loss
    associate_tree = parse_tree(
        Associate,
        """
        (
          (1[&host=a,contents='{"x","y"}'],2[&host=a,contents='{"x","z"}']),
          3[&host=a,contents='{"x","y","z"}']
        );
        """,
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, scaled_cost) == 5

    results = reconcile(setting, best_scaled_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            (
              (
                1[&host=a,contents='{"x","y"}']
              )[&kind=loss,host=a,contents='{"x","y","z"}',segment='{"z"}'],
              2[&host=a,contents='{"x","z"}']
            )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","z"}',result=1],
            3[&host=a,contents='{"x","y","z"}']
          )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y","z"}']
        )[&kind=gain,host=a,contents='set()',gained='{"x","y","z"}'];
        """,
        """
        (
          (
            (
              (
                1[&host=a,contents='{"x","y"}']
              )[&kind=loss,host=a,contents='{"x","y","__extra__"}',segment='{"__extra__"}'],
              2[&host=a,contents='{"x","z"}']
            )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","z"}',result=1],
            3[&host=a,contents='{"x","y","z"}']
          )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y","z"}']
        )[&kind=gain,host=a,contents='set()',gained='{"x","y","z"}'];
        """,
    )
    assert all(_history_match_input(setting, history) for history in results)

    # Three leaves in same host with disjoint contents: Cuts
    associate_tree = parse_tree(
        Associate,
        """
        (
          (1[&host=a,contents='{"x","y"}'],2[&host=a,contents='{"z"}']),
          3[&host=a,contents='{"x","y","z"}']
        );
        """,
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, scaled_cost) == 4.5

    results = reconcile(setting, best_scaled_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            (
              1[&host=a,contents='{"x","y"}'],
              2[&host=a,contents='{"z"}']
            )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y"}',cut=True],
            3[&host=a,contents='{"x","y","z"}']
          )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y","z"}']
        )[&kind=gain,host=a,contents='set()',gained='{"x","y","z"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)

    # Three leaves with disjoint contents and isolated genes: Cuts and gains
    associate_tree = parse_tree(
        Associate,
        """
        (
            (1[&host=a,contents='{"x","y","l"}'],2[&host=a,contents='{"z","r"}']),
            3[&host=a,contents='{"x","y","z","t"}']
        );
        """,
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, scaled_cost) == 4.5

    results = reconcile(setting, best_scaled_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            (
              (
                1[&host=a,contents='{"x","y","l"}']
              )[&kind=gain,host=a,contents='{"x","y"}',gained='{"l"}'],
              (
                2[&host=a,contents='{"z","r"}']
              )[&kind=gain,host=a,contents='{"z"}',gained='{"r"}']
            )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y"}',cut=True],
            (
              3[&host=a,contents='{"x","y","z","t"}']
            )[&kind=gain,host=a,contents='{"x","y","z"}',gained='{"t"}']
          )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y","z"}']
        )[&kind=gain,host=a,contents='set()',gained='{"x","y","z"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)


def test_reconcile_transfer():
    host_tree = parse_tree(Host, "((a,b)c,(d,e)f)g;")

    # Three leaves with a distant leaf: Transfer
    associate_tree = parse_tree(
        Associate,
        """
        (
          (1[&host=a,contents='{"x","y","z"}'],2[&host=d,contents='{"x"}']),
          3[&host=b,contents='{"x","y","z"}']
        );
        """,
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, scaled_cost) == 4

    results = reconcile(setting, best_scaled_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            (
              1[&host=a,contents='{"x","y","z"}'],
              2[&host=d,contents='{"x"}']
            )[&
              kind=diverge,
              host=a,
              contents='{"x","y","z"}',
              segment='{"x"}',
              transfer=True,
              result=1
            ],
            3[&host=b,contents='{"x","y","z"}']
          )[&kind=codiverge,host=c,contents='{"x","y","z"}']
        )[&kind=gain,host=c,contents='set()',gained='{"x","y","z"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)

    # Three leaves with a distant leaf and disjoint contents: Cut and transfer
    associate_tree = parse_tree(
        Associate,
        """
        (
          (1[&host=a,contents='{"y","z"}'],2[&host=d,contents='{"x"}']),
          3[&host=b,contents='{"x","y","z"}']
        );
        """,
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, scaled_cost) == 4.5

    results = reconcile(setting, best_scaled_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            (
              1[&host=a,contents='{"y","z"}'],
              2[&host=d,contents='{"x"}']
            )[&
              kind=diverge,
              host=a,
              contents='{"x","y","z"}',
              segment='{"x"}',
              transfer=True,
              cut=True,
              result=1
            ],
            3[&host=b,contents='{"x","y","z"}']
          )[&kind=codiverge,host=c,contents='{"x","y","z"}']
        )[&kind=gain,host=c,contents='set()',gained='{"x","y","z"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)


def test_reconcile_unsampled():
    host_tree = parse_tree(Host, "(((A,B)ab,C)abc,(X,(Y,Z)yz)xyz)root;")
    unsampled_host_tree = graft_unsampled_hosts(host_tree)

    associate_tree = parse_tree(
        Associate,
        """
        (
          (
            1[&contents='{"a"}',host=A],
            (
              2[&contents='{"a"}',host=X],
              (
                3[&contents='{"a"}',host=Y],
                4[&contents='{"a"}',host=Z]
              )
            )
          ),
          5[&contents='{"a"}',host=Z]
        );
        """,
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert reconcile(setting, scaled_cost) == 7

    results = reconcile(setting, best_scaled_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            (
              (
                [&kind=loss,host=C,contents='{"a"}',segment='{"a"}'],
                (
                  [&kind=loss,host=B,contents='{"a"}',segment='{"a"}'],
                  1[&contents='{"a"}',host=A]
                )[&kind=codiverge,host=ab,contents='{"a"}']
              )[&kind=codiverge,host=abc,contents='{"a"}'],
              (
                2[&contents='{"a"}',host=X],
                (
                  3[&contents='{"a"}',host=Y],
                  4[&contents='{"a"}',host=Z]
                )[&kind=codiverge,host=yz,contents='{"a"}']
              )[&kind=codiverge,host=xyz,contents='{"a"}']
            )[&kind=codiverge,host=root,contents='{"a"}'],
            (
              [&kind=loss,host=abc,contents='{"a"}',segment='{"a"}'],
              (
                [&kind=loss,host=X,contents='{"a"}',segment='{"a"}'],
                (
                  [&kind=loss,host=Y,contents='{"a"}',segment='{"a"}'],
                  5[&contents='{"a"}',host=Z]
                )[&kind=codiverge,host=yz,contents='{"a"}']
              )[&kind=codiverge,host=xyz,contents='{"a"}']
            )[&kind=codiverge,host=root,contents='{"a"}']
          )[&kind=diverge,host=root,contents='{"a"}',segment='{"a"}']
        )[&kind=gain,host=root,contents='set()',gained='{"a"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)

    setting = Reconciliation(unsampled_host_tree, associate_tree)
    assert reconcile(setting, scaled_cost) == 6

    results = reconcile(setting, best_scaled_cost).selected.value
    assert results == _tree_set(
        """
        (
          (
            (
              (
                [&contents='{"a"}',host='abc[U]'],
                (
                  [&kind=loss,contents='{"a"}',segment='{"a"}',host='C[P]'],
                  (
                    [&contents='{"a"}',host='ab[U]'],
                    (
                      [&kind=loss,contents='{"a"}',segment='{"a"}',host='B[P]'],
                      (
                        [&contents='{"a"}',host='A[U]'],
                        1[&contents='{"a"}',host=A]
                      )[&kind=codiverge,host='A[P]',contents='{"a"}']
                    )[&kind=codiverge,host=ab,contents='{"a"}']
                  )[&kind=codiverge,host='ab[P]',contents='{"a"}']
                )[&kind=codiverge,host=abc,contents='{"a"}']
              )[&kind=codiverge,host='abc[P]',contents='{"a"}'],
              (
                [&contents='{"a"}',host='xyz[U]'],
                (
                  (
                    [&contents='{"a"}',host='X[U]'],
                    2[&contents='{"a"}',host=X]
                  )[&kind=codiverge,host='X[P]',contents='{"a"}'],
                  (
                    [&contents='{"a"}',host='yz[U]'],
                    (
                      (
                        [&contents='{"a"}',host='Y[U]'],
                        3[&contents='{"a"}',host=Y]
                      )[&kind=codiverge,host='Y[P]',contents='{"a"}'],
                      (
                        [&contents='{"a"}',host='Z[U]'],
                        4[&contents='{"a"}',host=Z]
                      )[&kind=codiverge,host='Z[P]',contents='{"a"}']
                    )[&kind=codiverge,host=yz,contents='{"a"}']
                  )[&kind=codiverge,host='yz[P]',contents='{"a"}']
                )[&kind=codiverge,host=xyz,contents='{"a"}']
              )[&kind=codiverge,host='xyz[P]',contents='{"a"}']
            )[&kind=codiverge,host=root,contents='{"a"}'],
            (
              [&contents='{"a"}',host='root[U]'],
              5[&contents='{"a"}',host=Z]
            )[&kind=diverge,host='root[U]',contents='{"a"}',segment='{"a"}',transfer=True,result=1]
          )[&kind=codiverge,host='root[P]',contents='{"a"}']
        )[&kind=gain,host='root[P]',contents='set()',gained='{"a"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)
