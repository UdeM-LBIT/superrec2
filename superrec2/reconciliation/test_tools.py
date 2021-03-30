import unittest
from ete3 import PhyloTree
from ..utils.trees import LowestCommonAncestor
from .tools import (
    CostType,
    Event,
    get_reconciliation_cost,
    get_event,
    get_species_name,
    parse_reconciliation,
    serialize_reconciliation,
    parse_labeling,
    serialize_labeling,
    reconcile_leaves,
    reconcile_all,
)


class TestReconciliationTools(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gene_tree = PhyloTree(
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
            sp_naming_function=get_species_name,
            format=1,
        )

        cls.species_tree = PhyloTree("(((X,Y)XY,Z)XYZ,(W,T)WT)XYZWT;", format=1)
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

    def test_get_event(self):
        expected_events = {
            "1": Event.DUPLICATION,
            "2": Event.HORIZONTAL_GENE_TRANSFER,
            "3": Event.SPECIATION,
            "4": Event.DUPLICATION,
            "5": Event.SPECIATION,
            "6": Event.DUPLICATION,
            "7": Event.SPECIATION,
            "8": Event.SPECIATION,
            "9": Event.SPECIATION,
            "10": Event.DUPLICATION,
            "11": Event.DUPLICATION,
            "12": Event.SPECIATION,
            "13": Event.HORIZONTAL_GENE_TRANSFER,
            "14": Event.DUPLICATION,
        }

        for name, event in expected_events.items():
            self.assertEqual(
                get_event(
                    self.gene_tree & name,
                    self.species_lca,
                    self.rec,
                ),
                event,
            )

    def test_get_reconciliation_cost(self):
        expected_costs = {
            (1, 0, 0): 6,
            (0, 1, 0): 2,
            (0, 0, 1): 3,
            (1, 1, 1): 11,
        }

        for (dup, hgt, loss), value in expected_costs.items():
            self.assertEqual(
                get_reconciliation_cost(
                    self.gene_tree,
                    self.species_lca,
                    self.rec,
                    {
                        CostType.DUPLICATION: dup,
                        CostType.HORIZONTAL_GENE_TRANSFER: hgt,
                        CostType.FULL_LOSS: loss,
                    },
                ),
                value,
            )

    def test_reconcile_all(self):
        gene_tree = PhyloTree(
            "((x_1,x_2)2,(y_1,z_1)3)1;",
            sp_naming_function=get_species_name,
            format=1,
        )

        gene_1 = gene_tree & "1"
        gene_2 = gene_tree & "2"
        gene_3 = gene_tree & "3"

        species_tree = PhyloTree("(X,(Y,Z)YZ)XYZ;", format=1)

        species_x = species_tree & "X"
        species_y = species_tree & "Y"
        species_z = species_tree & "Z"
        species_yz = species_tree & "YZ"
        species_xyz = species_tree & "XYZ"

        species_lca = LowestCommonAncestor(species_tree)

        recs = list(reconcile_all(gene_tree, species_lca))
        leaves = reconcile_leaves(gene_tree, species_tree)

        # Check that all valid reconciliations are generated
        self.assertEqual(len(recs), 16)
        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_yz,
                gene_1: species_xyz,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_yz,
                gene_1: species_x,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_yz,
                gene_1: species_yz,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_xyz,
                gene_1: species_xyz,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_y,
                gene_1: species_yz,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_y,
                gene_1: species_xyz,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_y,
                gene_1: species_x,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_y,
                gene_1: species_y,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_z,
                gene_1: species_yz,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_z,
                gene_1: species_xyz,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_z,
                gene_1: species_x,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_x,
                gene_3: species_z,
                gene_1: species_z,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_xyz,
                gene_3: species_yz,
                gene_1: species_xyz,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_xyz,
                gene_3: species_xyz,
                gene_1: species_xyz,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_xyz,
                gene_3: species_y,
                gene_1: species_xyz,
            },
            recs,
        )

        self.assertIn(
            {
                **leaves,
                gene_2: species_xyz,
                gene_3: species_z,
                gene_1: species_xyz,
            },
            recs,
        )

        # Check that all the generated reconciliations are valid
        for rec in recs:
            for gene in rec:
                self.assertNotEqual(
                    get_event(gene, species_lca, rec),
                    Event.INVALID,
                )

    def test_serialize_reconciliation(self):
        serialized = (
            "1:XYZWT,2:XYZ,3:XYZ,4:W,5:XYZWT,6:XYZ,7:XY,"
            "8:XYZ,9:XY,10:Y,11:Y,12:WT,13:T,14:T"
        )

        self.assertEqual(
            parse_reconciliation(self.gene_tree, self.species_tree, serialized),
            self.rec_internal,
        )
        self.assertEqual(
            serialize_reconciliation(self.rec_internal),
            serialized,
        )

    def test_serialize_labeling(self):
        labeling = {
            self.gene_tree & "x_1": ["a"],
            self.gene_tree & "x_2": ["a", "b", "c", "d", "e", "f"],
            self.gene_tree & "x_3": ["d", "e", "f", "a", "b", "c"],
        }
        serialized = "x_1:a,x_2:abcdef,x_3:defabc"

        self.assertEqual(parse_labeling(self.gene_tree, serialized), labeling)
        self.assertEqual(serialize_labeling(labeling), serialized)
