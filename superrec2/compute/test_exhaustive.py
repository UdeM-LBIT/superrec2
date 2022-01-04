import unittest
from ete3 import Tree
from ..utils.trees import LowestCommonAncestor
from ..model.tree_mapping import get_species_mapping
from ..model.reconciliation import (
    ReconciliationInput,
    ReconciliationOutput,
    NodeEvent,
)
from .exhaustive import reconcile_all


class TestComputeExhaustive(unittest.TestCase):
    def test_reconcile_all(self):
        gene_tree = Tree("((x_1,x_2)2,(y_1,z_1)3)1;", format=1)
        species_tree = Tree("(X,(Y,Z)YZ)XYZ;", format=1)
        leaf_gene_species = get_species_mapping(gene_tree, species_tree)

        species_lca = LowestCommonAncestor(species_tree)
        rec_input = ReconciliationInput(
            gene_tree,
            species_lca,
            leaf_gene_species,
        )
        rec_outputs = list(reconcile_all(rec_input))

        # Check that all valid reconciliations are generated
        self.assertEqual(len(rec_outputs), 16)
        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "YZ",
                gene_tree & "1": species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "YZ",
                gene_tree & "1": species_tree & "X",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "YZ",
                gene_tree & "1": species_tree & "YZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "XYZ",
                gene_tree & "1": species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Y",
                gene_tree & "1": species_tree & "YZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Y",
                gene_tree & "1": species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Y",
                gene_tree & "1": species_tree & "X",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Y",
                gene_tree & "1": species_tree & "Y",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Z",
                gene_tree & "1": species_tree & "YZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Z",
                gene_tree & "1": species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Z",
                gene_tree & "1": species_tree & "X",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Z",
                gene_tree & "1": species_tree & "Z",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "XYZ",
                gene_tree & "3": species_tree & "YZ",
                gene_tree & "1": species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "XYZ",
                gene_tree & "3": species_tree & "XYZ",
                gene_tree & "1": species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "XYZ",
                gene_tree & "3": species_tree & "Y",
                gene_tree & "1": species_tree & "XYZ",
            }),
            rec_outputs,
        )

        self.assertIn(
            ReconciliationOutput(rec_input, {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "XYZ",
                gene_tree & "3": species_tree & "Z",
                gene_tree & "1": species_tree & "XYZ",
            }),
            rec_outputs,
        )

        # Check that all the generated reconciliations are valid
        for rec in rec_outputs:
            for gene in gene_tree:
                self.assertNotEqual(
                    rec.node_event(gene),
                    NodeEvent.INVALID,
                )
