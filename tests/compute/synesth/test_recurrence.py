from superrec2.model.history import (
    parse_tree,
    graft_unsampled_hosts,
    Associate,
    Host,
    Reconciliation,
    Event,
    History,
)
from superrec2.utils.algebras import Structure, MinPlus
from superrec2.compute.synesth import (
    solve_binary,
    EventCosts,
    HistoryBuilder,
    history_generator,
)


_unit_cost = EventCosts()
_scaled_cost = EventCosts(
    loss=1,
    duplication=2,
    cut=2.5,
    transfer_duplication=4,
    transfer_cut=4.5,
)


min_unit_cost = Structure(MinPlus, _unit_cost.morphism)
min_scaled_cost = Structure(MinPlus, _scaled_cost.morphism)

best_unit_cost = min_unit_cost * history_generator
best_scaled_cost = min_scaled_cost * history_generator


def _build_event_tree(source):
    return HistoryBuilder(parse_tree(Event, source))


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
    assert solve_binary(setting, min_unit_cost).value == 0

    cost, results = solve_binary(setting, best_unit_cost).value
    assert cost == 0
    assert results == _tree_set(
        """
        (
          (
            1[&host=a,contents='{"x"}',apparent=True],
            2[&host=b,contents='{"x"}',apparent=True]
          )[&kind=codiverge,host=c,contents='{"x"}',apparent=True]
        )[&kind=gain,host=c,contents='set()',gained='{"x"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)

    # Two leaves with subsumed contents in divergent hosts: Speciation with gain
    associate_tree = parse_tree(
        Associate, '(1[&host=a,contents=\'{"x","y"}\'],2[&host=b,contents=\'{"x"}\']);'
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert solve_binary(setting, min_unit_cost).value == 0

    cost, results = solve_binary(setting, best_unit_cost).value
    assert cost == 0
    assert results == _tree_set(
        """
        (
          (
            (1[&host=a,contents='{"x","y"}',apparent=True])
            [&kind=gain,host=a,contents='{"x"}',gained='{"y"}'],
            2[&host=b,contents='{"x"}',apparent=True]
          )[&kind=codiverge,host=c,contents='{"x"}',apparent=True]
        )[&kind=gain,host=c,contents='set()',gained='{"x"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)

    # Two leaves with disjoint contents in divergent hosts: Empty speciation with gains
    associate_tree = parse_tree(
        Associate, "(1[&host=a,contents='{\"x\"}'],2[&host=b,contents='{\"y\"}']);"
    )
    setting = Reconciliation(host_tree, associate_tree)
    assert solve_binary(setting, min_unit_cost).value == 0

    cost, results = solve_binary(setting, best_unit_cost).value
    assert cost == 0
    assert results == _tree_set(
        """
        (
          (1[&host=a,contents='{"x"}',apparent=True])
          [&kind=gain,host=a,contents='set()',gained='{"x"}'],
          (2[&host=b,contents='{"y"}',apparent=True])
          [&kind=gain,host=b,contents='set()',gained='{"y"}']
        )[&kind=codiverge,host=c,contents='set()',apparent=True];
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
    assert solve_binary(setting, min_scaled_cost).value == 2

    cost, results = solve_binary(setting, best_unit_cost).value
    assert cost == 2
    assert results == _tree_set(
        """
        (
          (
            (
              (
                (
                  1[&host=a,contents='{"x","y"}',apparent=True]
                )[&kind=loss,host=a,contents='{"x","y","__extra__"}',segment='{"__extra__"}'],
                (
                  2[&host=b,contents='{"y","z"}',apparent=True]
                )[&kind=loss,host=b,contents='{"y","z","__extra__"}',segment='{"__extra__"}']
              )[&kind=codiverge,host=c,contents='{"x","y","z","__extra__"}',apparent=True],
              3[&host=d,contents='{"w","x","y","z"}',apparent=True]
            )[&kind=codiverge,host=e,contents='{"w","x","y","z"}',apparent=True],
            4[&host=f,contents='{"w","x","y","z"}',apparent=True]
          )[&kind=codiverge,host=g,contents='{"w","x","y","z"}',apparent=True]
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
    assert solve_binary(setting, min_scaled_cost).value == 4

    cost, results = solve_binary(setting, best_scaled_cost).value
    assert cost == 4
    assert results == _tree_set(
        """
        (
          (
            (
              1[&host=a,contents='{"x","y"}',apparent=True],
              2[&host=a,contents='{"x","y","z"}',apparent=True]
            )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y"}',apparent=True],
            3[&host=a,contents='{"x","y","z"}',apparent=True]
          )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y","z"}',apparent=True]
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
    assert solve_binary(setting, min_scaled_cost).value == 5

    cost, results = solve_binary(setting, best_scaled_cost).value
    assert cost == 5
    assert results == _tree_set(
        """
        (
          (
            (
              (
                1[&host=a,contents='{"x","y"}',apparent=True]
              )[&kind=loss,host=a,contents='{"x","y","z"}',segment='{"z"}'],
              2[&host=a,contents='{"x","z"}',apparent=True]
            )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","z"}',result=1,apparent=True],
            3[&host=a,contents='{"x","y","z"}',apparent=True]
          )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y","z"}',apparent=True]
        )[&kind=gain,host=a,contents='set()',gained='{"x","y","z"}'];
        """,
        """
        (
          (
            (
              (
                1[&host=a,contents='{"x","y"}',apparent=True]
              )[&kind=loss,host=a,contents='{"x","y","__extra__"}',segment='{"__extra__"}'],
              2[&host=a,contents='{"x","z"}',apparent=True]
            )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","z"}',result=1,apparent=True],
            3[&host=a,contents='{"x","y","z"}',apparent=True]
          )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y","z"}',apparent=True]
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
    assert solve_binary(setting, min_scaled_cost).value == 4.5

    cost, results = solve_binary(setting, best_scaled_cost).value
    assert cost == 4.5
    assert results == _tree_set(
        """
        (
          (
            (
              1[&host=a,contents='{"x","y"}',apparent=True],
              2[&host=a,contents='{"z"}',apparent=True]
            )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y"}',cut=True,apparent=True],
            3[&host=a,contents='{"x","y","z"}',apparent=True]
          )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y","z"}',apparent=True]
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
    assert solve_binary(setting, min_scaled_cost).value == 4.5

    cost, results = solve_binary(setting, best_scaled_cost).value
    assert cost == 4.5
    assert results == _tree_set(
        """
        (
          (
            (
              (
                1[&host=a,contents='{"x","y","l"}',apparent=True]
              )[&kind=gain,host=a,contents='{"x","y"}',gained='{"l"}'],
              (
                2[&host=a,contents='{"z","r"}',apparent=True]
              )[&kind=gain,host=a,contents='{"z"}',gained='{"r"}']
            )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y"}',cut=True,apparent=True],
            (
              3[&host=a,contents='{"x","y","z","t"}',apparent=True]
            )[&kind=gain,host=a,contents='{"x","y","z"}',gained='{"t"}']
          )[&kind=diverge,host=a,contents='{"x","y","z"}',segment='{"x","y","z"}',apparent=True]
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
    assert solve_binary(setting, min_scaled_cost).value == 4

    cost, results = solve_binary(setting, best_scaled_cost).value
    assert cost == 4
    assert results == _tree_set(
        """
        (
          (
            (
              1[&host=a,contents='{"x","y","z"}',apparent=True],
              2[&host=d,contents='{"x"}',apparent=True]
            )[&
              kind=diverge,
              host=a,
              contents='{"x","y","z"}',
              segment='{"x"}',
              transfer=True,
              result=1,
              apparent=True
            ],
            3[&host=b,contents='{"x","y","z"}',apparent=True]
          )[&kind=codiverge,host=c,contents='{"x","y","z"}',apparent=True]
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
    assert solve_binary(setting, min_scaled_cost).value == 4.5

    cost, results = solve_binary(setting, best_scaled_cost).value
    assert cost == 4.5
    assert results == _tree_set(
        """
        (
          (
            (
              1[&host=a,contents='{"y","z"}',apparent=True],
              2[&host=d,contents='{"x"}',apparent=True]
            )[&
              kind=diverge,
              host=a,
              contents='{"x","y","z"}',
              segment='{"x"}',
              transfer=True,
              cut=True,
              result=1,
              apparent=True
            ],
            3[&host=b,contents='{"x","y","z"}',apparent=True]
          )[&kind=codiverge,host=c,contents='{"x","y","z"}',apparent=True]
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
    assert solve_binary(setting, min_scaled_cost).value == 7

    cost, results = solve_binary(setting, best_scaled_cost).value
    assert cost == 7
    assert results == _tree_set(
        """
        (
          (
            (
              (
                [&kind=loss,host=C,contents='{"a"}',segment='{"a"}'],
                (
                  [&kind=loss,host=B,contents='{"a"}',segment='{"a"}'],
                  1[&contents='{"a"}',host=A,apparent=True]
                )[&kind=codiverge,host=ab,contents='{"a"}']
              )[&kind=codiverge,host=abc,contents='{"a"}'],
              (
                2[&contents='{"a"}',host=X,apparent=True],
                (
                  3[&contents='{"a"}',host=Y,apparent=True],
                  4[&contents='{"a"}',host=Z,apparent=True]
                )[&kind=codiverge,host=yz,contents='{"a"}',apparent=True]
              )[&kind=codiverge,host=xyz,contents='{"a"}',apparent=True]
            )[&kind=codiverge,host=root,contents='{"a"}',apparent=True],
            (
              [&kind=loss,host=abc,contents='{"a"}',segment='{"a"}'],
              (
                [&kind=loss,host=X,contents='{"a"}',segment='{"a"}'],
                (
                  [&kind=loss,host=Y,contents='{"a"}',segment='{"a"}'],
                  5[&contents='{"a"}',host=Z,apparent=True]
                )[&kind=codiverge,host=yz,contents='{"a"}']
              )[&kind=codiverge,host=xyz,contents='{"a"}']
            )[&kind=codiverge,host=root,contents='{"a"}']
          )[&kind=diverge,host=root,contents='{"a"}',segment='{"a"}',apparent=True]
        )[&kind=gain,host=root,contents='set()',gained='{"a"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)

    setting = Reconciliation(unsampled_host_tree, associate_tree)
    assert solve_binary(setting, min_scaled_cost).value == 6

    cost, results = solve_binary(setting, best_scaled_cost).value
    assert cost == 6
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
                        1[&contents='{"a"}',host=A,apparent=True]
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
                    2[&contents='{"a"}',host=X,apparent=True]
                  )[&kind=codiverge,host='X[P]',contents='{"a"}'],
                  (
                    [&contents='{"a"}',host='yz[U]'],
                    (
                      (
                        [&contents='{"a"}',host='Y[U]'],
                        3[&contents='{"a"}',host=Y,apparent=True]
                      )[&kind=codiverge,host='Y[P]',contents='{"a"}'],
                      (
                        [&contents='{"a"}',host='Z[U]'],
                        4[&contents='{"a"}',host=Z,apparent=True]
                      )[&kind=codiverge,host='Z[P]',contents='{"a"}']
                    )[&kind=codiverge,host=yz,contents='{"a"}',apparent=True]
                  )[&kind=codiverge,host='yz[P]',contents='{"a"}']
                )[&kind=codiverge,host=xyz,contents='{"a"}',apparent=True]
              )[&kind=codiverge,host='xyz[P]',contents='{"a"}']
            )[&kind=codiverge,host=root,contents='{"a"}',apparent=True],
            (
              [&contents='{"a"}',host='root[U]'],
              5[&contents='{"a"}',host=Z,apparent=True]
            )[&kind=diverge,host='root[U]',contents='{"a"}',segment='{"a"}',transfer=True,result=1]
          )[&kind=codiverge,host='root[P]',contents='{"a"}',apparent=True]
        )[&kind=gain,host='root[P]',contents='set()',gained='{"a"}'];
        """
    )
    assert all(_history_match_input(setting, history) for history in results)
