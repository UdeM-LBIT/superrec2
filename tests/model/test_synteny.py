from ete3 import Tree
from superrec2.model.synteny import (
    sort_synteny,
    format_synteny,
    parse_synteny_mapping,
    serialize_synteny_mapping,
)


tree = Tree("((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8;", format=1)


def test_sort_synteny():
    assert sort_synteny("dates") == ["a", "d", "e", "s", "t"]
    assert format_synteny("dates", width=21) == "d, a, t, e, s"
    assert sort_synteny(set("dates")) == ["a", "d", "e", "s", "t"]
    assert format_synteny(set("dates"), width=21) == "a, d, e, s, t"
    assert sort_synteny(["gene1", "gene7", "gene10", "gene7a1"]) == [
        "gene1",
        "gene7",
        "gene7a1",
        "gene10",
    ]
    assert (
        format_synteny(["gene1", "gene7", "gene10", "gene7a1"], width=21)
        == "gene1, gene7,\ngene10, gene7a1"
    )
    assert sort_synteny({"gene1", "gene7", "gene10", "gene7a1"}) == [
        "gene1",
        "gene7",
        "gene7a1",
        "gene10",
    ]
    assert (
        format_synteny({"gene1", "gene7", "gene10", "gene7a1"}, width=21)
        == "gene1, gene7,\ngene7a1, gene10"
    )


def test_parse_synteny_mapping():
    assert parse_synteny_mapping(
        tree,
        {
            "x_3": "abcdef",
            "z_2": ["cas1", "a", "b", "c", "cas3"],
        },
    ) == {
        tree & "x_3": "abcdef",
        tree & "z_2": ["cas1", "a", "b", "c", "cas3"],
    }


def test_serialize_synteny_mapping():
    assert serialize_synteny_mapping(
        {
            tree & "x_3": ["a", "b", "c", "d", "e", "f"],
            tree & "z_2": ["cas1", "a", "b", "c", "cas3"],
        }
    ) == {
        "x_3": ["a", "b", "c", "d", "e", "f"],
        "z_2": ["cas1", "a", "b", "c", "cas3"],
    }
    assert serialize_synteny_mapping(
        {
            tree & "x_3": ["a", "b", "c", "d", "e", "f"],
            tree & "z_2": {"cas1", "a", "b", "cas10", "cas2"},
        }
    ) == {
        "x_3": ["a", "b", "c", "d", "e", "f"],
        "z_2": ["a", "b", "cas1", "cas2", "cas10"],
    }
