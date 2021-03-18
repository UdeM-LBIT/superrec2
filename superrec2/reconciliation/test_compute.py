import unittest
from infinity import inf
from ete3 import PhyloTree
from superrec2 import reconciliation as rec
from superrec2.utils.lowest_common_ancestor import LowestCommonAncestor


class TestReconciliationCompute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gene_tree = PhyloTree(
            "((x_1,(x_2,(y_1,z_1)5)4)3,(y_2,z_2)2)1;",
            sp_naming_function=rec.tools.get_species_name,
            format=1,
        )

        cls.species_tree = PhyloTree("(X,(Y,Z)YZ)XYZ;", format=1)
        cls.species_lca = LowestCommonAncestor(cls.species_tree)

        cls.all_recs = list(
            rec.tools.reconcile_all(
                cls.gene_tree, cls.species_tree, cls.species_lca
            )
        )

    def test_all_recs(self):
        self.assertEqual(len(self.all_recs), 199)

    def test_reconcile_lca(self):
        costs = {
            rec.tools.CostType.DUPLICATION: 1,
            rec.tools.CostType.HORIZONTAL_GENE_TRANSFER: inf,
            rec.tools.CostType.LOSS: 1,
        }

        result = rec.compute.reconcile_lca(self.gene_tree, self.species_lca)

        # Check that the expected result is returned
        self.assertEqual(
            result,
            {
                **rec.tools.reconcile_leaves(self.gene_tree, self.species_tree),
                self.gene_tree & "1": self.species_tree & "XYZ",
                self.gene_tree & "2": self.species_tree & "YZ",
                self.gene_tree & "3": self.species_tree & "XYZ",
                self.gene_tree & "4": self.species_tree & "XYZ",
                self.gene_tree & "5": self.species_tree & "YZ",
            },
        )

        result_cost = rec.tools.get_reconciliation_cost(
            self.gene_tree,
            self.species_lca,
            result,
            costs,
        )

        self.assertEqual(result_cost, 4)
        self.assertIn(result, self.all_recs)

        # Check optimality
        for possible_rec in self.all_recs:
            cost = rec.tools.get_reconciliation_cost(
                self.gene_tree, self.species_lca, possible_rec, costs
            )
            self.assertTrue(cost >= result_cost)

            if cost == result_cost:
                self.assertEqual(possible_rec, result)

    def test_reconcile_thl(self):
        costs = {
            rec.tools.CostType.DUPLICATION: 1,
            rec.tools.CostType.HORIZONTAL_GENE_TRANSFER: 1,
            rec.tools.CostType.LOSS: 1,
        }

        result_cost, results = rec.compute.reconcile_thl(
            self.gene_tree, self.species_lca, costs
        )

        # Check that all results have the same expected cost
        self.assertEqual(result_cost, 2)

        for result in results:
            self.assertEqual(
                rec.tools.get_reconciliation_cost(
                    self.gene_tree,
                    self.species_lca,
                    result,
                    costs,
                ),
                result_cost,
            )

        # Check that all expected results are returned
        self.assertEqual(len(results), 2)
        self.assertIn(
            {
                **rec.tools.reconcile_leaves(self.gene_tree, self.species_tree),
                self.gene_tree & "1": self.species_tree & "XYZ",
                self.gene_tree & "2": self.species_tree & "YZ",
                self.gene_tree & "3": self.species_tree & "X",
                self.gene_tree & "4": self.species_tree & "X",
                self.gene_tree & "5": self.species_tree & "YZ",
            },
            results,
        )

        self.assertIn(
            {
                **rec.tools.reconcile_leaves(self.gene_tree, self.species_tree),
                self.gene_tree & "1": self.species_tree & "XYZ",
                self.gene_tree & "2": self.species_tree & "YZ",
                self.gene_tree & "3": self.species_tree & "X",
                self.gene_tree & "4": self.species_tree & "YZ",
                self.gene_tree & "5": self.species_tree & "YZ",
            },
            results,
        )

        for result in results:
            self.assertIn(result, self.all_recs)

        # Check optimality
        for possible_rec in self.all_recs:
            cost = rec.tools.get_reconciliation_cost(
                self.gene_tree, self.species_lca, possible_rec, costs
            )
            self.assertTrue(cost >= result_cost)

            if cost == result_cost:
                self.assertIn(possible_rec, results)
