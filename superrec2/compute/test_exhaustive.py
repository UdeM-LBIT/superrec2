import unittest
from ete3 import Tree
from ..utils.trees import LowestCommonAncestor
from ..model.tree_mapping import get_species_mapping
from ..model.reconciliation import (
    ReconciliationInput,
    ReconciliationOutput,
    NodeEvent,
)
from .exhaustive import generate_all, exhaustive_any, exhaustive_all


class TestComputeExhaustive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gene_tree = Tree("((x_1,x_2)2,(y_1,z_1)3)1;", format=1)
        cls.species_tree = Tree("(X,(Y,Z)YZ)XYZ;", format=1)
        cls.leaf_gene_species = get_species_mapping(
            cls.gene_tree, cls.species_tree
        )

        cls.species_lca = LowestCommonAncestor(cls.species_tree)
        cls.rec_input = ReconciliationInput(
            cls.gene_tree,
            cls.species_lca,
            cls.leaf_gene_species,
        )

    def test_generate_all(self):
        rec_outputs = list(generate_all(self.rec_input))

        # Check that all valid reconciliations are generated
        self.assertEqual(len(rec_outputs), 16)
        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "YZ",
                self.gene_tree & "1": self.species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "YZ",
                self.gene_tree & "1": self.species_tree & "X",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "YZ",
                self.gene_tree & "1": self.species_tree & "YZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "XYZ",
                self.gene_tree & "1": self.species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "Y",
                self.gene_tree & "1": self.species_tree & "YZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "Y",
                self.gene_tree & "1": self.species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "Y",
                self.gene_tree & "1": self.species_tree & "X",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "Y",
                self.gene_tree & "1": self.species_tree & "Y",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "Z",
                self.gene_tree & "1": self.species_tree & "YZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "Z",
                self.gene_tree & "1": self.species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "Z",
                self.gene_tree & "1": self.species_tree & "X",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "X",
                self.gene_tree & "3": self.species_tree & "Z",
                self.gene_tree & "1": self.species_tree & "Z",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "XYZ",
                self.gene_tree & "3": self.species_tree & "YZ",
                self.gene_tree & "1": self.species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "XYZ",
                self.gene_tree & "3": self.species_tree & "XYZ",
                self.gene_tree & "1": self.species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "XYZ",
                self.gene_tree & "3": self.species_tree & "Y",
                self.gene_tree & "1": self.species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "2": self.species_tree & "XYZ",
                self.gene_tree & "3": self.species_tree & "Z",
                self.gene_tree & "1": self.species_tree & "XYZ",
            }),
            rec_outputs,
        )

        # Check that all the generated reconciliations are valid
        for rec in rec_outputs:
            for gene in self.gene_tree:
                self.assertNotEqual(
                    rec.node_event(gene),
                    NodeEvent.INVALID,
                )

    def test_exhaustive_any_all(self):
        all_outputs = list(generate_all(self.rec_input))
        min_cost = min(output.cost() for output in all_outputs)
        any_min = exhaustive_any(self.rec_input)
        all_min = exhaustive_all(self.rec_input)

        self.assertIn(any_min, all_outputs)
        self.assertIn(any_min, all_min)
        self.assertEqual(any_min.cost(), min_cost)

        for a_min in all_min:
            self.assertEqual(a_min.cost(), min_cost)
            self.assertIn(a_min, all_outputs)
