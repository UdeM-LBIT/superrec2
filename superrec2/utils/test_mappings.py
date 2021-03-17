from ete3 import Tree
import unittest
from .mappings import invert_mapping


class TestUtilsMappings(unittest.TestCase):
    def test_invert_injective(self):
        mapping = {"A": 1, "Z": 0, "B": 2, "C": 3}
        invert = invert_mapping(mapping)
        self.assertEqual(invert, {1: ["A"], 2: ["B"], 3: ["C"], 0: ["Z"]})

    def test_invert_surjective(self):
        mapping = {"A": 1, "B": 1, "C": 2, "D": 1, "E": 2, "Z": 3}
        invert = invert_mapping(mapping)
        self.assertEqual(invert, {1: ["A", "B", "D"], 2: ["C", "E"], 3: ["Z"]})

    def test_invert_outside_domain(self):
        mapping = {"A": 1, "B": 2}
        invert = invert_mapping(mapping)
        self.assertEqual(invert[3], [])

    def test_invert_treenodes(self):
        tree_1 = Tree("((x_1,y_1)2,(((x_2,y_2)5,z_1)4,z_2)3)1;", format=8)
        tree_2 = Tree("((x,y)xy,z)xyz;", format=8)
        mapping = {
            tree_1 & "x_1": tree_2 & "x",
            tree_1 & "x_2": tree_2 & "x",
            tree_1 & "y_1": tree_2 & "y",
            tree_1 & "y_2": tree_2 & "y",
            tree_1 & "z_1": tree_2 & "z",
            tree_1 & "z_2": tree_2 & "z",
            tree_1 & "1": tree_2,
            tree_1 & "2": tree_2 & "xy",
            tree_1 & "3": tree_2,
            tree_1 & "4": tree_2 & "xy",
            tree_1 & "5": tree_2 & "xy",
        }
        invert = invert_mapping(mapping)
        self.assertEqual(
            invert,
            {
                tree_2: [tree_1 & "1", tree_1 & "3"],
                tree_2 & "xy": [tree_1 & "2", tree_1 & "4", tree_1 & "5"],
                tree_2 & "x": [tree_1 & "x_1", tree_1 & "x_2"],
                tree_2 & "y": [tree_1 & "y_1", tree_1 & "y_2"],
                tree_2 & "z": [tree_1 & "z_1", tree_1 & "z_2"],
            },
        )
