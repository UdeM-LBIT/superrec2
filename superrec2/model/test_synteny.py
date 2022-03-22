import unittest
from ete3 import Tree
from .synteny import (
    sort_synteny,
    format_synteny,
    parse_synteny_mapping,
    serialize_synteny_mapping,
)


class TestModelSynteny(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tree = Tree("((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8;", format=1)

    def test_sort_synteny(self):
        self.assertEqual(
            sort_synteny("dates"),
            ["a", "d", "e", "s", "t"],
        )
        self.assertEqual(
            format_synteny("dates", width=21),
            "d, a, t, e, s",
        )
        self.assertEqual(
            sort_synteny(set("dates")),
            ["a", "d", "e", "s", "t"],
        )
        self.assertEqual(
            format_synteny(set("dates"), width=21),
            "a, d, e, s, t",
        )
        self.assertEqual(
            sort_synteny(["gene1", "gene7", "gene10", "gene7a1"]),
            ["gene1", "gene7", "gene7a1", "gene10"],
        )
        self.assertEqual(
            format_synteny(["gene1", "gene7", "gene10", "gene7a1"], width=21),
            "gene1, gene7,\ngene10, gene7a1",
        )
        self.assertEqual(
            sort_synteny({"gene1", "gene7", "gene10", "gene7a1"}),
            ["gene1", "gene7", "gene7a1", "gene10"],
        )
        self.assertEqual(
            format_synteny({"gene1", "gene7", "gene10", "gene7a1"}, width=21),
            "gene1, gene7,\ngene7a1, gene10",
        )

    def test_parse_synteny_mapping(self):
        self.assertEqual(
            parse_synteny_mapping(
                self.tree,
                {
                    "x_3": "abcdef",
                    "z_2": ["cas1", "a", "b", "c", "cas3"],
                },
            ),
            {
                self.tree & "x_3": "abcdef",
                self.tree & "z_2": ["cas1", "a", "b", "c", "cas3"],
            },
        )

    def test_serialize_synteny_mapping(self):
        self.assertEqual(
            serialize_synteny_mapping(
                {
                    self.tree & "x_3": ["a", "b", "c", "d", "e", "f"],
                    self.tree & "z_2": ["cas1", "a", "b", "c", "cas3"],
                }
            ),
            {
                "x_3": ["a", "b", "c", "d", "e", "f"],
                "z_2": ["cas1", "a", "b", "c", "cas3"],
            },
        )
        self.assertEqual(
            serialize_synteny_mapping(
                {
                    self.tree & "x_3": ["a", "b", "c", "d", "e", "f"],
                    self.tree & "z_2": {"cas1", "a", "b", "cas10", "cas2"},
                }
            ),
            {
                "x_3": ["a", "b", "c", "d", "e", "f"],
                "z_2": ["a", "b", "cas1", "cas2", "cas10"],
            },
        )
