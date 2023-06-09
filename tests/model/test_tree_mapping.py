import unittest
from ete3 import Tree
from superrec2.model.tree_mapping import (
    parse_tree_mapping,
    get_species_mapping,
    serialize_tree_mapping,
)


class TestReconciliationTools(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gene_tree = Tree(
            """
            (
                ((x_1,z_1)3,(w_1,w_2)4)2,
                (
                    (
                        (x_2,y_4)7,
                        ((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8
                    )6,
                    (w_3,(z_3,(t_1,t_2)14)13)12
                )5
            )1;
        """,
            format=1,
        )
        cls.species_tree = Tree("(((X,Y)XY,Z)XYZ,(W,T)WT)XYZWT;", format=1)

    def test_parse_tree_mapping(self):
        self.assertEqual(
            parse_tree_mapping(
                self.gene_tree,
                self.species_tree,
                {
                    "1": "XYZWT",
                    "2": "XYZ",
                    "3": "XYZ",
                    "4": "W",
                    "5": "XYZWT",
                    "6": "XYZ",
                    "7": "XY",
                    "8": "XYZ",
                    "9": "XY",
                    "10": "Y",
                    "11": "Y",
                    "12": "WT",
                    "13": "T",
                    "14": "T",
                },
            ),
            {
                self.gene_tree & "1": self.species_tree & "XYZWT",
                self.gene_tree & "2": self.species_tree & "XYZ",
                self.gene_tree & "3": self.species_tree & "XYZ",
                self.gene_tree & "4": self.species_tree & "W",
                self.gene_tree & "5": self.species_tree & "XYZWT",
                self.gene_tree & "6": self.species_tree & "XYZ",
                self.gene_tree & "7": self.species_tree & "XY",
                self.gene_tree & "8": self.species_tree & "XYZ",
                self.gene_tree & "9": self.species_tree & "XY",
                self.gene_tree & "10": self.species_tree & "Y",
                self.gene_tree & "11": self.species_tree & "Y",
                self.gene_tree & "12": self.species_tree & "WT",
                self.gene_tree & "13": self.species_tree & "T",
                self.gene_tree & "14": self.species_tree & "T",
            },
        )

    def test_get_species_mapping(self):
        self.assertEqual(
            get_species_mapping(self.gene_tree, self.species_tree),
            {
                self.gene_tree & "x_1": self.species_tree & "X",
                self.gene_tree & "x_2": self.species_tree & "X",
                self.gene_tree & "x_3": self.species_tree & "X",
                self.gene_tree & "y_1": self.species_tree & "Y",
                self.gene_tree & "y_2": self.species_tree & "Y",
                self.gene_tree & "y_3": self.species_tree & "Y",
                self.gene_tree & "y_4": self.species_tree & "Y",
                self.gene_tree & "z_1": self.species_tree & "Z",
                self.gene_tree & "z_2": self.species_tree & "Z",
                self.gene_tree & "z_3": self.species_tree & "Z",
                self.gene_tree & "w_1": self.species_tree & "W",
                self.gene_tree & "w_2": self.species_tree & "W",
                self.gene_tree & "w_3": self.species_tree & "W",
                self.gene_tree & "t_1": self.species_tree & "T",
                self.gene_tree & "t_2": self.species_tree & "T",
            },
        )

    def test_serialize_tree_mapping(self):
        self.assertEqual(
            serialize_tree_mapping(
                {
                    self.gene_tree & "x_1": self.species_tree & "X",
                    self.gene_tree & "2": self.species_tree & "XYZ",
                    self.gene_tree & "x_2": self.species_tree & "X",
                    self.gene_tree & "3": self.species_tree & "XYZ",
                }
            ),
            {
                "2": "XYZ",
                "3": "XYZ",
                "x_1": "X",
                "x_2": "X",
            },
        )
