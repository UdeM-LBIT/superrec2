import unittest
from ete3 import PhyloTree
from superrec2 import reconciliation as rec
from superrec2.utils.lowest_common_ancestor import LowestCommonAncestor


class TestReconciliationCompute(unittest.TestCase):
    def test_reconcile_lca(self):
        gene_tree = PhyloTree(
            "((x_1,(x_2,(y_1,z_1)5)4)3,(y_2,z_2)2)1;",
            sp_naming_function=rec.tools.get_species_name, format=1
        )

        species_tree = PhyloTree("(X,(Y,Z)YZ)XYZ;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        result = rec.compute.reconcile_lca(gene_tree, species_lca)

        self.assertEqual(
            result,
            {
                **rec.tools.reconcile_leaves(gene_tree, species_tree),
                gene_tree & "1": species_tree & "XYZ",
                gene_tree & "2": species_tree & "YZ",
                gene_tree & "3": species_tree & "XYZ",
                gene_tree & "4": species_tree & "XYZ",
                gene_tree & "5": species_tree & "YZ",
            }
        )

        self.assertEqual(
            rec.tools.get_reconciliation_cost(
                gene_tree,
                species_lca,
                result,
                {
                    rec.tools.CostType.Duplication: 1,
                    rec.tools.CostType.HorizontalGeneTransfer: 1,
                    rec.tools.CostType.Loss: 1,
                }
            ),
            4
        )

    def test_reconcile_thl(self):
        gene_tree = PhyloTree(
            "((x_1,(x_2,(y_1,z_1)5)4)3,(y_2,z_2)2)1;",
            sp_naming_function=rec.tools.get_species_name, format=1
        )

        species_tree = PhyloTree("(X,(Y,Z)YZ)XYZ;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        costs = {
            rec.tools.CostType.Duplication: 1,
            rec.tools.CostType.HorizontalGeneTransfer: 1,
            rec.tools.CostType.Loss: 1,
        }
        result_cost, results = rec.compute.reconcile_thl(
            gene_tree,
            species_lca,
            costs
        )

        self.assertEqual(result_cost, 2)
        self.assertEqual(len(results), 2)

        for result in results:
            self.assertEqual(rec.tools.get_reconciliation_cost(
                gene_tree,
                species_lca,
                result,
                costs,
            ), result_cost)

        self.assertIn(
            {
                **rec.tools.reconcile_leaves(gene_tree, species_tree),
                gene_tree & "1": species_tree & "XYZ",
                gene_tree & "2": species_tree & "YZ",
                gene_tree & "3": species_tree & "X",
                gene_tree & "4": species_tree & "X",
                gene_tree & "5": species_tree & "YZ",
            },
            results
        )

        self.assertIn(
            {
                **rec.tools.reconcile_leaves(gene_tree, species_tree),
                gene_tree & "1": species_tree & "XYZ",
                gene_tree & "2": species_tree & "YZ",
                gene_tree & "3": species_tree & "X",
                gene_tree & "4": species_tree & "YZ",
                gene_tree & "5": species_tree & "YZ",
            },
            results
        )
