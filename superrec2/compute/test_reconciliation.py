import unittest
from infinity import inf
from ete3 import Tree
from ..utils.trees import LowestCommonAncestor
from ..model.tree_mapping import get_species_mapping
from ..model.reconciliation import (
    ReconciliationInput,
    ReconciliationOutput,
    NodeEvent,
    EdgeEvent,
)
from .exhaustive import reconcile_all
from .reconciliation import reconcile_lca, reconcile_thl_any, reconcile_thl_all


class TestComputeReconciliation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gene_tree = Tree(
            "((x_1,(x_2,(y_1,z_1)5)4)3,(y_2,z_2)2)1;",
            format=1,
        )

        cls.species_tree = Tree("(X,(Y,Z)YZ)XYZ;", format=1)
        cls.species_lca = LowestCommonAncestor(cls.species_tree)
        cls.leaf_gene_species = get_species_mapping(
            cls.gene_tree, cls.species_tree
        )

        cls.rec_input = ReconciliationInput(
            cls.gene_tree,
            cls.species_lca,
            cls.leaf_gene_species,
        )

        cls.all_outputs = list(reconcile_all(cls.rec_input))

    def test_output_count(self):
        self.assertEqual(len(self.all_outputs), 199)

    def test_reconcile_lca(self):
        self.rec_input.costs[NodeEvent.DUPLICATION] = 1
        self.rec_input.costs[NodeEvent.HORIZONTAL_TRANSFER] = inf
        self.rec_input.costs[EdgeEvent.FULL_LOSS] = 1

        result = reconcile_lca(self.rec_input)

        # Check that the expected result is returned
        self.assertEqual(
            result.object_species,
            {
                **self.leaf_gene_species,
                self.gene_tree & "1": self.species_tree & "XYZ",
                self.gene_tree & "2": self.species_tree & "YZ",
                self.gene_tree & "3": self.species_tree & "XYZ",
                self.gene_tree & "4": self.species_tree & "XYZ",
                self.gene_tree & "5": self.species_tree & "YZ",
            },
        )

        self.assertEqual(result.cost(), 4)
        self.assertIn(result, self.all_outputs)

        # Check optimality
        for possible_rec in self.all_outputs:
            cost = possible_rec.cost()
            self.assertTrue(cost >= result.cost())

            if cost == result.cost():
                self.assertEqual(possible_rec, result)

    def test_reconcile_thl(self):
        self.rec_input.costs[NodeEvent.DUPLICATION] = 1
        self.rec_input.costs[NodeEvent.HORIZONTAL_TRANSFER] = 1
        self.rec_input.costs[EdgeEvent.FULL_LOSS] = 1

        any_result = reconcile_thl_any(self.rec_input)
        results = list(reconcile_thl_all(self.rec_input))

        # Check that any is part of all results
        self.assertIn(any_result, results)

        # Check that all expected results are returned
        self.assertEqual(len(results), 2)
        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "1": self.species_tree & "XYZ",
                self.gene_tree & "2": self.species_tree & "YZ",
                self.gene_tree & "3": self.species_tree & "X",
                self.gene_tree & "4": self.species_tree & "X",
                self.gene_tree & "5": self.species_tree & "YZ",
            }),
            results,
        )

        self.assertIn(
            ReconciliationOutput(self.rec_input, {
                **self.leaf_gene_species,
                self.gene_tree & "1": self.species_tree & "XYZ",
                self.gene_tree & "2": self.species_tree & "YZ",
                self.gene_tree & "3": self.species_tree & "X",
                self.gene_tree & "4": self.species_tree & "YZ",
                self.gene_tree & "5": self.species_tree & "YZ",
            }),
            results,
        )

        # Check that all results have the same expected cost
        for result in results:
            self.assertEqual(result.cost(), 2)

        for result in results:
            self.assertIn(result, self.all_outputs)

        # Check optimality
        for possible_rec in self.all_outputs:
            cost = possible_rec.cost()
            self.assertTrue(cost >= 2)

            if cost == 2:
                self.assertIn(possible_rec, results)
