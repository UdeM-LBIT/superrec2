from math import inf
from sowing.indexed import IndexedTree
from superrec2.model.history import parse_tree, Host, Event, Extant
from superrec2.compute.superdtlx.contents import EXTRA_CONTENTS
from superrec2.compute.superdtlx.paths import (
    make_codiv_path,
    make_transfer_path,
    make_path,
)
from superrec2.compute.util import (
    EventCosts,
    make_cost_algebra,
    history_builder,
    history_generator,
)


unit_cost = make_cost_algebra("unit_cost", costs=EventCosts())
host_tree = parse_tree(
    Host,
    """
        (
            (
                (
                    (4,2LU[&sampled=False])2L,
                    (
                        ((6,5LU[&sampled=False])5L,(7,5RU[&sampled=False])5R)5,
                        2RU[&sampled=False]
                    )2R
                )2,
                1LU[&sampled=False]
            )1L,
            (3,1RU[&sampled=False])1R
        )1;
    """,
)
host_index = IndexedTree(host_tree)


def _build_event_tree(source):
    return history_builder(parse_tree(Event, source))


def _tree_set(*trees):
    return frozenset(map(_build_event_tree, trees))


def test_make_codiv_path():
    # Towards ancestor: No solution
    leaf = Extant(name="x", host="1", contents=frozenset("abc"))
    args = {
        "start_host": "6",
        "end_host": "1",
        "contents": frozenset("abc"),
        "host_index": host_index,
    }

    assert (
        make_codiv_path(**args, structure=unit_cost, path=unit_cost.make(leaf)).value
        == inf
    )

    assert (
        make_codiv_path(
            **args,
            structure=history_generator,
            path=history_generator.make(leaf),
        ).value
        == frozenset()
    )

    # Separate hosts: No solution
    leaf = Extant(name="x", host="7", contents=frozenset("abc"))
    args = {
        "start_host": "6",
        "end_host": "7",
        "contents": frozenset("abc"),
        "host_index": host_index,
    }

    assert (
        make_codiv_path(
            **args,
            structure=unit_cost,
            path=unit_cost.make(leaf),
        ).value
        == inf
    )

    assert (
        make_codiv_path(
            **args,
            structure=history_generator,
            path=history_generator.make(leaf),
        ).value
        == frozenset()
    )

    # Equal hosts: No events needed
    leaf = Extant(name="x", host="6", contents=frozenset("abc"))
    args = {
        "start_host": "6",
        "end_host": "6",
        "contents": frozenset("abc"),
        "host_index": host_index,
    }

    assert (
        make_codiv_path(
            **args,
            structure=unit_cost,
            path=unit_cost.make(leaf),
        ).value
        == 0
    )

    assert make_codiv_path(
        **args,
        structure=history_generator,
        path=history_generator.make(leaf),
    ).value == frozenset({history_builder.make(leaf)})

    # Different hosts: Codivergence path
    leaf = Extant(name="x", host="6", contents=frozenset("abc"))
    args = {
        "start_host": "1",
        "end_host": "6",
        "contents": frozenset("abc"),
        "host_index": host_index,
    }

    assert (
        make_codiv_path(
            **args,
            structure=unit_cost,
            path=unit_cost.make(leaf),
        ).value
        == 3
    )

    assert make_codiv_path(
        **args,
        structure=history_generator,
        path=history_generator.make(leaf),
    ).value == _tree_set(
        """
        (
            [&kind=loss,host=1R,contents='{"a","b","c"}',segment='{"a","b","c"}'],
            (
                [&kind=extant,host=1LU,contents='{"a","b","c"}'],
                (
                    [&
                        kind=loss,
                        host=2L,
                        contents='{"a","b","c"}',
                        segment='{"a","b","c"}'
                    ],
                    (
                        [&kind=extant,host=2RU,contents='{"a","b","c"}'],
                        (
                            [&
                                kind=loss,
                                host=5R,
                                contents='{"a","b","c"}',
                                segment='{"a","b","c"}'
                            ],
                            (
                                [&kind=extant,host=5LU,contents='{"a","b","c"}'],
                                x[&kind=extant,host=6,contents='{"a","b","c"}']
                            )[&kind=codiverge,host=5L,contents='{"a","b","c"}']
                        )[&kind=codiverge,host=5,contents='{"a","b","c"}']
                    )[&kind=codiverge,host=2R,contents='{"a","b","c"}']
                )[&kind=codiverge,host=2,contents='{"a","b","c"}']
            )[&kind=codiverge,host=1L,contents='{"a","b","c"}']
        )[&kind=codiverge,host=1,contents='{"a","b","c"}'];
        """
    )


def test_make_transfer_path():
    # Towards ancestor: No solution
    leaf = Extant(name="x", host="1", contents=frozenset("abc"))
    args = {
        "start_host": "6",
        "end_host": "1",
        "start_contents": frozenset("abc"),
        "end_contents": frozenset("abc"),
        "host_index": host_index,
    }

    assert (
        make_transfer_path(
            **args,
            structure=unit_cost,
            path=unit_cost.make(leaf),
        ).value
        == inf
    )

    assert (
        make_transfer_path(
            **args,
            structure=history_generator,
            path=history_generator.make(leaf),
        ).value
        == frozenset()
    )

    # Separate sampled hosts, same contents: Direct transfer-cut
    leaf = Extant(name="x", host="7", contents=frozenset("abc"))
    args = {
        "start_host": "6",
        "end_host": "7",
        "start_contents": frozenset("abc"),
        "end_contents": frozenset("abc"),
        "host_index": host_index,
    }

    assert (
        make_transfer_path(
            **args,
            structure=unit_cost,
            path=unit_cost.make(leaf),
        ).value
        == 1
    )

    assert make_transfer_path(
        **args,
        structure=history_generator,
        path=history_generator.make(leaf),
    ).value == _tree_set(
        """
        (x[&kind=extant,host=7,contents='{"a","b","c"}'])
        [&
            kind=diverge,
            host=6,
            contents='{"a","b","c"}',
            segment='{"a","b","c"}',
            transfer=True,
            cut=True
        ];
        """,
        """
        (
            [&kind=loss,host=6,contents='{"a","b","c"}',segment='{"a","b","c"}'],
            x[&kind=extant,host=7,contents='{"a","b","c"}']
        )
        [&
            kind=diverge,
            host=6,
            contents='{"a","b","c"}',
            segment='{"a","b","c"}',
            transfer=True,
            result=1
        ];
        """,
    )

    # Separate sampled hosts, subset of contents: Direct transfer-cut or copy
    leaf = Extant(name="x", host="7", contents=frozenset("ab"))
    args = {
        "start_host": "6",
        "end_host": "7",
        "start_contents": frozenset("abc"),
        "end_contents": frozenset("ab"),
        "host_index": host_index,
    }

    cost = make_transfer_path(**args, structure=unit_cost, path=unit_cost.make(leaf))
    assert cost.value == 2

    assert make_transfer_path(
        **args,
        structure=history_generator,
        path=history_generator.make(leaf),
    ).value == _tree_set(
        """
        (
            [&kind=loss,host=6,contents='{"c"}',segment='{"c"}'],
            x[&kind=extant,host=7,contents='{"a","b"}']
        )
        [&
            kind=diverge,
            host=6,
            contents='{"a","b","c"}',
            segment='{"a","b"}',
            transfer=True,
            cut=True,
            result=1
        ];
        """,
        """
        (
            [&kind=loss,host=6,contents='{"a","b","c"}',segment='{"a","b","c"}'],
            x[&kind=extant,host=7,contents='{"a","b"}']
        )
        [&
            kind=diverge,
            host=6,
            contents='{"a","b","c"}',
            segment='{"a","b"}',
            transfer=True,
            result=1
        ];
        """,
    )

    # Separate sampled hosts, different contents: No solution
    leaf = Extant(name="x", host="7", contents=frozenset("bcd"))
    args = {
        "start_host": "6",
        "end_host": "7",
        "start_contents": frozenset("abc"),
        "end_contents": frozenset("bcd"),
        "host_index": host_index,
    }

    assert (
        make_transfer_path(
            **args,
            structure=unit_cost,
            path=unit_cost.make(leaf),
        ).value
        == inf
    )

    assert (
        make_transfer_path(
            **args,
            structure=history_generator,
            path=history_generator.make(leaf),
        ).value
        == frozenset()
    )

    # Equal hosts: No solution
    leaf = Extant(name="x", host="6", contents=frozenset("abc"))
    args = {
        "start_host": "6",
        "end_host": "6",
        "start_contents": frozenset("abc"),
        "end_contents": frozenset("abc"),
        "host_index": host_index,
    }

    assert (
        make_transfer_path(
            **args,
            structure=unit_cost,
            path=unit_cost.make(leaf),
        ).value
        == inf
    )

    assert (
        make_transfer_path(
            **args,
            structure=history_generator,
            path=history_generator.make(leaf),
        ).value
        == frozenset()
    )

    # Sampled host towards descendant: Start with codivergence path
    leaf = Extant(name="x", host="7", contents=frozenset("ab"))
    args = {
        "start_host": "1",
        "end_host": "7",
        "start_contents": frozenset("abc"),
        "end_contents": frozenset("ab"),
        "host_index": host_index,
    }

    assert (
        make_transfer_path(
            **args,
            structure=unit_cost,
            path=unit_cost.make(leaf),
        ).value
        == 3
    )

    assert make_transfer_path(
        **args,
        structure=history_generator,
        path=history_generator.make(leaf),
    ).value == _tree_set(
        """
        (
            [&kind=loss,host=1L,contents='{"a","b","c"}',segment='{"a","b","c"}'],
            (
                [&kind=loss,host=1R,contents='{"c"}',segment='{"c"}'],
                x[&kind=extant,host=7,contents='{"a","b"}']
            )
            [&
                kind=diverge,
                host=1R,
                contents='{"a","b","c"}',
                segment='{"a","b"}',
                transfer=True,
                cut=True,
                result=1
            ]
        )[&kind=codiverge,host=1,contents='{"a","b","c"}'];
        """,
        """
        (
            [&kind=loss,host=1L,contents='{"a","b","c"}',segment='{"a","b","c"}'],
            (
                [&kind=loss,host=1R,contents='{"a","b","c"}',segment='{"a","b","c"}'],
                x[&kind=extant,host=7,contents='{"a","b"}']
            )
            [&
                kind=diverge,
                host=1R,
                contents='{"a","b","c"}',
                segment='{"a","b"}',
                transfer=True,
                result=1
            ]
        )[&kind=codiverge,host=1,contents='{"a","b","c"}'];
        """,
    )

    # Above unsampled towards sampled: No loss before transfer
    leaf = Extant(name="x", host="7", contents=frozenset("ab"))
    args = {
        "start_host": "1L",
        "end_host": "7",
        "start_contents": frozenset("abc"),
        "end_contents": frozenset("ab"),
        "host_index": host_index,
    }

    assert (
        make_transfer_path(
            **args,
            structure=unit_cost,
            path=unit_cost.make(leaf),
        ).value
        == 2
    )

    assert make_transfer_path(
        **args,
        structure=history_generator,
        path=history_generator.make(leaf),
    ).value == _tree_set(
        """
        (
            [&kind=loss,host=2,contents='{"a","b","c"}',segment='{"a","b","c"}'],
            (
                [&kind=extant,host=1LU,contents='{"c"}'],
                x[&kind=extant,host=7,contents='{"a","b"}']
            )
            [&
                kind=diverge,
                host=1LU,
                contents='{"a","b","c"}',
                segment='{"a","b"}',
                transfer=True,
                cut=True,
                result=1
            ]
        )[&kind=codiverge,host=1L,contents='{"a","b","c"}'];
        """,
        """
        (
            [&kind=loss,host=2,contents='{"a","b","c"}',segment='{"a","b","c"}'],
            (
                [&kind=extant,host=1LU,contents='{"a","b","c"}'],
                x[&kind=extant,host=7,contents='{"a","b"}']
            )
            [&
                kind=diverge,
                host=1LU,
                contents='{"a","b","c"}',
                segment='{"a","b"}',
                transfer=True,
                result=1
            ]
        )[&kind=codiverge,host=1L,contents='{"a","b","c"}'];
        """,
    )


def test_make_path():
    # Towards ancestor: No solution
    leaf = Extant(name="x", host="1", contents=frozenset("abc"))
    args = {
        "start_host": "6",
        "end_host": "1",
        "start_contents": frozenset("abc"),
        "end_contents": frozenset("abc"),
        "host_index": host_index,
    }

    assert (
        make_path(**args, structure=unit_cost, path=unit_cost.make(leaf)).value == inf
    )

    assert (
        make_path(
            **args, structure=history_generator, path=history_generator.make(leaf)
        ).value
        == frozenset()
    )

    # Missing extra contents with equal remainder: No solution
    leaf = Extant(name="x", host="6", contents=frozenset({"a", "b", EXTRA_CONTENTS}))
    args = {
        "start_host": "6",
        "end_host": "6",
        "start_contents": frozenset("ab"),
        "end_contents": frozenset({"a", "b", EXTRA_CONTENTS}),
        "host_index": host_index,
    }

    assert (
        make_path(**args, structure=unit_cost, path=unit_cost.make(leaf)).value == inf
    )

    assert (
        make_path(
            **args, structure=history_generator, path=history_generator.make(leaf)
        ).value
        == frozenset()
    )

    # Missing extra contents with contents to lose: No events needed
    leaf = Extant(name="x", host="6", contents=frozenset({"a", EXTRA_CONTENTS}))
    args = {
        "start_host": "6",
        "end_host": "6",
        "start_contents": frozenset("ab"),
        "end_contents": frozenset({"a", EXTRA_CONTENTS}),
        "host_index": host_index,
    }

    assert make_path(**args, structure=unit_cost, path=unit_cost.make(leaf)).value == 0

    assert make_path(
        **args, structure=history_generator, path=history_generator.make(leaf)
    ).value == frozenset({history_builder.make(leaf)})

    # Kept extra contents with contents to lose: No events needed
    leaf = Extant(name="x", host="6", contents=frozenset({"a", EXTRA_CONTENTS}))
    args = {
        "start_host": "6",
        "end_host": "6",
        "start_contents": frozenset({"a", "b", EXTRA_CONTENTS}),
        "end_contents": frozenset({"a", EXTRA_CONTENTS}),
        "host_index": host_index,
    }

    assert make_path(**args, structure=unit_cost, path=unit_cost.make(leaf)).value == 0

    assert make_path(
        **args, structure=history_generator, path=history_generator.make(leaf)
    ).value == frozenset({history_builder.make(leaf)})

    # Gained contents in same host: Gain events
    leaf = Extant(name="x", host="6", contents=frozenset("abcd"))
    args = {
        "start_host": "6",
        "end_host": "6",
        "start_contents": frozenset("ab"),
        "end_contents": frozenset("abcd"),
        "host_index": host_index,
    }

    assert make_path(**args, structure=unit_cost, path=unit_cost.make(leaf)).value == 0

    assert make_path(
        **args, structure=history_generator, path=history_generator.make(leaf)
    ).value == _tree_set(
        """
        (x[&kind=extant,host=6,contents='{"a","b","c","d"}'])
        [&kind=gain,host=6,contents='{"a","b"}',gained='{"c","d"}'];
        """,
    )

    # Lost contents in same host: Loss event
    leaf = Extant(name="x", host="6", contents=frozenset("ab"))
    args = {
        "start_host": "6",
        "end_host": "6",
        "start_contents": frozenset({"a", "b", EXTRA_CONTENTS}),
        "end_contents": frozenset("ab"),
        "host_index": host_index,
    }

    assert make_path(**args, structure=unit_cost, path=unit_cost.make(leaf)).value == 1

    assert make_path(
        **args, structure=history_generator, path=history_generator.make(leaf)
    ).value == _tree_set(
        """
        (x[&kind=extant,host=6,contents='{"a","b"}'])
        [&kind=loss,host=6,contents='{"a","b","__extra__"}',segment='{"__extra__"}'];
        """,
    )

    # Gained and lost contents in same host: Loss then gains
    leaf = Extant(name="x", host="6", contents=frozenset("abc"))
    args = {
        "start_host": "6",
        "end_host": "6",
        "start_contents": frozenset({"a", "b", EXTRA_CONTENTS}),
        "end_contents": frozenset("abc"),
        "host_index": host_index,
    }

    assert make_path(**args, structure=unit_cost, path=unit_cost.make(leaf)).value == 1

    assert make_path(
        **args, structure=history_generator, path=history_generator.make(leaf)
    ).value == _tree_set(
        """
        (
            (x[&kind=extant,host=6,contents='{"a","b","c"}'])
            [&kind=gain,host=6,contents='{"a","b"}',gained='{"c"}']
        )
        [&kind=loss,host=6,contents='{"a","b","__extra__"}',segment='{"__extra__"}'];
        """,
    )

    # Gained and lost contents in different hosts: Switch, (loss), then gains
    leaf = Extant(name="x", host="6", contents=frozenset("abc"))
    args = {
        "start_host": "2R",
        "end_host": "6",
        "start_contents": frozenset({"a", "b", EXTRA_CONTENTS}),
        "end_contents": frozenset("abc"),
        "host_index": host_index,
    }

    assert make_path(**args, structure=unit_cost, path=unit_cost.make(leaf)).value == 2

    assert make_path(
        **args, structure=history_generator, path=history_generator.make(leaf)
    ).value == _tree_set(
        """
        (
            [&kind=extant,host=2RU,contents='{"a","b","__extra__"}'],
            (
                [&
                    kind=loss,
                    host=5R,
                    contents='{"a","b","__extra__"}',
                    segment='{"a","b","__extra__"}'
                ],
                (
                    [&kind=extant,host=5LU,contents='{"a","b","__extra__"}'],
                    (
                        (x[&kind=extant,host=6,contents='{"a","b","c"}'])
                        [&kind=gain,host=6,contents='{"a","b"}',gained='{"c"}']
                    )[&kind=loss,host=6,contents='{"a","b","__extra__"}',segment='{"__extra__"}']
                )
                [&kind=codiverge,host=5L,contents='{"a","b","__extra__"}']
            )
            [&kind=codiverge,host=5,contents='{"a","b","__extra__"}']
        )
        [&kind=codiverge,host=2R,contents='{"a","b","__extra__"}'];
        """,
        """
        (
            [&
                kind=loss,
                host=5,
                contents='{"a","b","__extra__"}',
                segment='{"a","b","__extra__"}'
            ],
            (
                [&kind=extant,host=2RU,contents='{"a","b","__extra__"}'],
                (x[&kind=extant,host=6,contents='{"a","b","c"}'])
                [&kind=gain,host=6,contents='{"a","b"}',gained='{"c"}']
            )
            [&
                kind=diverge,
                host=2RU,
                contents='{"a","b","__extra__"}',
                segment='{"a","b"}',
                transfer=True,
                result=1
            ]
        )
        [&kind=codiverge,host=2R,contents='{"a","b","__extra__"}'];
        """,
        """
        (
            [&
                kind=loss,
                host=5,
                contents='{"a","b","__extra__"}',
                segment='{"a","b","__extra__"}'
            ],
            (
                [&kind=extant,host=2RU,contents='{"__extra__"}'],
                (x[&kind=extant,host=6,contents='{"a","b","c"}'])
                [&kind=gain,host=6,contents='{"a","b"}',gained='{"c"}']
            )
            [&
                kind=diverge,
                host=2RU,
                contents='{"a","b","__extra__"}',
                segment='{"a","b"}',
                transfer=True,
                cut=True,
                result=1
            ]
        )
        [&kind=codiverge,host=2R,contents='{"a","b","__extra__"}'];
        """,
    )
