import unittest
from ete3 import PhyloTree
from .lowest_common_ancestor import LowestCommonAncestor
from .reconciliation import (
    CostType,
    Event,
    get_cost,
    get_event,
    get_species_name,
    parse_reconciliation,
    reconcile_leaves,
)

class TestReconciliationTools(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gene_tree = PhyloTree("""
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
        """, sp_naming_function=get_species_name, format=1)

        cls.species_tree = PhyloTree(
            "(((X,Y)XY,Z)XYZ,(W,T)WT)XYZWT;",
            format=1
        )
        cls.species_lca = LowestCommonAncestor(cls.species_tree)

        cls.rec_leaves = {
            cls.gene_tree & "x_1": cls.species_tree & "X",
            cls.gene_tree & "x_2": cls.species_tree & "X",
            cls.gene_tree & "x_3": cls.species_tree & "X",
            cls.gene_tree & "y_1": cls.species_tree & "Y",
            cls.gene_tree & "y_2": cls.species_tree & "Y",
            cls.gene_tree & "y_3": cls.species_tree & "Y",
            cls.gene_tree & "y_4": cls.species_tree & "Y",
            cls.gene_tree & "z_1": cls.species_tree & "Z",
            cls.gene_tree & "z_2": cls.species_tree & "Z",
            cls.gene_tree & "z_3": cls.species_tree & "Z",
            cls.gene_tree & "w_1": cls.species_tree & "W",
            cls.gene_tree & "w_2": cls.species_tree & "W",
            cls.gene_tree & "w_3": cls.species_tree & "W",
            cls.gene_tree & "t_1": cls.species_tree & "T",
            cls.gene_tree & "t_2": cls.species_tree & "T",
        }

        cls.rec_internal = {
            cls.gene_tree & "1": cls.species_tree & "XYZWT",
            cls.gene_tree & "2": cls.species_tree & "XYZ",
            cls.gene_tree & "3": cls.species_tree & "XYZ",
            cls.gene_tree & "4": cls.species_tree & "W",
            cls.gene_tree & "5": cls.species_tree & "XYZWT",
            cls.gene_tree & "6": cls.species_tree & "XYZ",
            cls.gene_tree & "7": cls.species_tree & "XY",
            cls.gene_tree & "8": cls.species_tree & "XYZ",
            cls.gene_tree & "9": cls.species_tree & "XY",
            cls.gene_tree & "10": cls.species_tree & "Y",
            cls.gene_tree & "11": cls.species_tree & "Y",
            cls.gene_tree & "12": cls.species_tree & "WT",
            cls.gene_tree & "13": cls.species_tree & "T",
            cls.gene_tree & "14": cls.species_tree & "T",
        }

        cls.rec = {
            **cls.rec_leaves,
            **cls.rec_internal,
        }

    def test_species_name(self):
        self.assertEqual(get_species_name("x_1"), "X")
        self.assertEqual(get_species_name("x_123"), "X")
        self.assertEqual(get_species_name("x__123"), "X")
        self.assertEqual(get_species_name("xyz_123"), "XYZ")
        self.assertEqual(get_species_name("1_x"), "1")
        self.assertEqual(get_species_name("X_X"), "X")

    def test_reconcile_leaves(self):
        self.assertEqual(
            reconcile_leaves(self.gene_tree, self.species_tree),
            self.rec_leaves,
        )

    def test_parse(self):
        self.assertEqual(
            parse_reconciliation(
                self.gene_tree, self.species_tree,
                "1:XYZWT,2:XYZ,3:XYZ,4:W,5:XYZWT,6:XYZ,7:XY,"
                "8:XYZ,9:XY,10:Y,11:Y,12:WT,13:T,14:T"
            ),
            self.rec_internal,
        )

    def test_get_event(self):
        expected_events = {
            "1": Event.Duplication,
            "2": Event.HorizontalGeneTransfer,
            "3": Event.Speciation,
            "4": Event.Duplication,
            "5": Event.Speciation,
            "6": Event.Duplication,
            "7": Event.Speciation,
            "8": Event.Speciation,
            "9": Event.Speciation,
            "10": Event.Duplication,
            "11": Event.Duplication,
            "12": Event.Speciation,
            "13": Event.HorizontalGeneTransfer,
            "14": Event.Duplication,
        }

        for name, event in expected_events.items():
            self.assertEqual(get_event(
                self.gene_tree & name,
                self.species_lca, self.rec,
            ), event)

    def test_get_cost(self):
        expected_costs = {
            (1, 0, 0): 6,
            (0, 1, 0): 2,
            (0, 0, 1): 3,
            (1, 1, 1): 11,
        }

        for (dup, hgt, loss), value in expected_costs.items():
            self.assertEqual(get_cost(
                self.gene_tree,
                self.species_lca,
                self.rec,
                {
                    CostType.Duplication: dup,
                    CostType.HorizontalGeneTransfer: hgt,
                    CostType.Loss: loss,
                }
            ), value)
