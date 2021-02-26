import unittest
from ete3 import PhyloTree
from infinity import inf
from superrec2 import reconciliation as rec
from superrec2.utils.lowest_common_ancestor import LowestCommonAncestor


class TestReconciliationSyntenies(unittest.TestCase):
    def assertLabelingEquals(
        self,
        gene_tree,
        species_lca,
        rec_result,
        input_labeling,
        expected_cost,
        expected_labelings,
    ):
        output = rec.syntenies.label_ancestral_syntenies(
            gene_tree, species_lca, rec_result, input_labeling
        )

        self.assertEqual(output[0], expected_cost)
        self.assertEqual(len(output[1]), len(expected_labelings))

        for labeling in expected_labelings:
            self.assertEqual(
                rec.syntenies.get_labeling_cost(
                    gene_tree, species_lca, rec_result, labeling
                ),
                expected_cost,
            )
            self.assertIn(labeling, output[1])

    def test_speciations(self):
        gene_tree = PhyloTree(
            "((x_1,y_1)2,z_1)1;",
            format=1, sp_naming_function=rec.tools.get_species_name
        )

        species_tree = PhyloTree("((X,Y)XY,Z)XYZ;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        rec_result = rec.compute.reconcile_lca(gene_tree, species_lca)

        # Test 1: Single optimal solution, conserved root
        input_1 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("de"),
            gene_tree & "z_1": list("abcde"),
        }

        self.assertLabelingEquals(
            gene_tree, species_lca, rec_result, input_1,
            expected_cost=2,
            expected_labelings=[
                {
                    **input_1,
                    gene_tree & "2": list("abcde"),
                    gene_tree & "1": list("abcde"),
                }
            ]
        )

        # Test 2: Single optimal solution, altered root
        input_2 = {
            gene_tree & "x_1": list("abd"),
            gene_tree & "y_1": list("bde"),
            gene_tree & "z_1": list("abcde"),
        }

        self.assertLabelingEquals(
            gene_tree, species_lca, rec_result, input_2,
            expected_cost=3,
            expected_labelings=[
                {
                    **input_2,
                    gene_tree & "2": list("abde"),
                    gene_tree & "1": list("abcde"),
                }
            ]
        )

        # Test 3: Multiple optimal solutions
        input_3 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("cd"),
            gene_tree & "z_1": list("be"),
        }

        self.assertLabelingEquals(
            gene_tree, species_lca, rec_result, input_3,
            expected_cost=4,
            expected_labelings=[
                {
                    **input_3,
                    gene_tree & "2": list("cdab"),
                    gene_tree & "1": list("cdabe"),
                },
                {
                    **input_3,
                    gene_tree & "2": list("abecd"),
                    gene_tree & "1": list("abecd"),
                },
            ]
        )

        # Test 4: No solution
        input_4 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("cd"),
            gene_tree & "z_1": list("ba"),
        }

        self.assertLabelingEquals(
            gene_tree, species_lca, rec_result, input_4,
            expected_cost=inf,
            expected_labelings=[],
        )
