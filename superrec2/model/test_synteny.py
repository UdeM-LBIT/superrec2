import unittest
from ete3 import Tree
from .synteny import (
    parse_synteny, serialize_synteny,
    parse_synteny_mapping, serialize_synteny_mapping
)


class TestModelSynteny(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tree = Tree("((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8;", format=1)
        cls.comma_tree = cls.tree.copy()
        (cls.comma_tree & "x_3").name = "x,3"
        (cls.comma_tree & "y_1").name = "y,1"
        (cls.comma_tree & "y_2").name = "y,2"
        (cls.comma_tree & "y_3").name = "y,3"
        (cls.comma_tree & "z_2").name = "z,2"

    def test_parse_synteny(self):
        self.assertEqual(
            parse_synteny(r"abcdef"),
            ["a", "b", "c", "d", "e", "f"],
        )
        self.assertEqual(
            parse_synteny(r"(cas6)(cas8)(cas7)(cas5)(cas3)(cas4)(cas1)(cas2)"),
            ["cas6", "cas8", "cas7", "cas5", "cas3", "cas4", "cas1", "cas2"],
        )
        self.assertEqual(
            parse_synteny(r"(g1) (g2) (g3)"),
            ["g1", "g2", "g3"],
        )
        self.assertEqual(
            parse_synteny(r"(cas3')(esc\(ape)(test\\gene)"),
            ["cas3'", "esc(ape", r"test\gene"],
        )
        self.assertEqual(
            parse_synteny(r"  ( spaces ) (	are	)   ( ignored   ) "),
            ["spaces", "are", "ignored"],
        )
        self.assertEqual(
            parse_synteny(r"(balanced(paren))(work)(corr(ect)ly)"),
            ["balanced(paren)", "work", "corr(ect)ly"],
        )
        self.assertEqual(
            parse_synteny(r"a\((b)(unclosed"),
            ["a", "(", "b"],
        )

    def test_serialize_synteny(self):
        self.assertEqual(
            serialize_synteny(["a", "b", "c", "d", "e", "f"]),
            r"abcdef"
        )
        self.assertEqual(
            serialize_synteny(["a", "(", ")", "\\"]),
            r"a\(\)\\",
        )
        self.assertEqual(
            serialize_synteny(["cas11", "cas7", "cas5"]),
            r"(cas11)(cas7)(cas5)",
        )

    def test_parse_synteny_mapping(self):
        self.assertEqual(
            parse_synteny_mapping(
                self.tree,
                "x_3 : abcdef, z_2 : (cas1)abc(cas3),",
            ),
            {
                self.tree & "x_3": ["a", "b", "c", "d", "e", "f"],
                self.tree & "z_2": ["cas1", "a", "b", "c", "cas3"],
            }
        )
        self.assertEqual(
            parse_synteny_mapping(
                self.tree,
                "x_3:abcdef,z_2:(cas1)abc(cas3)",
            ),
            {
                self.tree & "x_3": ["a", "b", "c", "d", "e", "f"],
                self.tree & "z_2": ["cas1", "a", "b", "c", "cas3"],
            }
        )
        self.assertEqual(
            parse_synteny_mapping(
                self.tree,
                "x_3 : abcdef, #z_2 : (cas1)abc(cas3),",
            ),
            {
                self.tree & "x_3": ["a", "b", "c", "d", "e", "f"],
            }
        )
        self.assertEqual(
            parse_synteny_mapping(
                self.comma_tree,
                "x\,3 : abc\,def, z\,2 : a(abc\, def)f",
            ),
            {
                self.comma_tree & "x,3": ["a", "b", "c", ",", "d", "e", "f"],
                self.comma_tree & "z,2": ["a", "abc, def", "f"],
            }
        )

    def test_serialize_synteny_mapping(self):
        self.assertEqual(
            serialize_synteny_mapping({
                self.tree & "x_3": ["a", "b", "c", "d", "e", "f"],
                self.tree & "z_2": ["cas1", "a", "b", "c", "cas3"],
            }),
            r"x_3:abcdef,z_2:(cas1)abc(cas3)",
        )
        self.assertEqual(
            serialize_synteny_mapping({
                self.comma_tree & "x,3": ["a", "b", "c", ",", "d", "e", "f"],
                self.comma_tree & "z,2": ["a", "abc, def", "f"],
            }),
            r"x\,3:abc\,def,z\,2:a(abc\, def)f",
        )
