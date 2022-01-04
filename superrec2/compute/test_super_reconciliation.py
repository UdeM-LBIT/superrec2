import unittest
from ete3 import Tree
from infinity import inf
from ..model.reconciliation import (
    ReconciliationInput,
    SuperReconciliationInput,
    SuperReconciliationOutput,
    get_default_cost,
)
from ..utils.trees import LowestCommonAncestor
from ..model.tree_mapping import get_species_mapping
from .reconciliation import reconcile_lca
from .super_reconciliation import label_syntenies_any, label_syntenies_all


class TestComputeSuperReconciliation(unittest.TestCase):
    def assertLabelingEquals(
        self,
        srec_input,
        rec_output,
        expected_cost,
        expected_labelings,
    ):
        any_l = label_syntenies_any(srec_input, rec_output)
        all_l = list(label_syntenies_all(srec_input, rec_output))
        expected_results = [
            SuperReconciliationOutput(
                input=srec_input,
                object_species=rec_output.object_species,
                syntenies=expected_labeling,
            )
            for expected_labeling in expected_labelings
        ]

        self.assertCountEqual(all_l, expected_results)

        if expected_results:
            self.assertEqual(any_l.cost(), expected_cost)
            self.assertIn(any_l, all_l)

            for labeling in expected_results:
                self.assertEqual(labeling.cost(), expected_cost)
        else:
            self.assertIsNone(any_l)

    def test_speciations(self):
        gene_tree = Tree("((x_1,y_1)2,z_1)1;", format=1)
        species_tree = Tree("((X,Y)XY,Z)XYZ;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        leaf_gene_species = get_species_mapping(gene_tree, species_tree)

        rec_input = ReconciliationInput(
            gene_tree,
            species_lca,
            leaf_gene_species,
        )
        rec_output = reconcile_lca(rec_input)

        # Test 1: Single optimal solution, conserved root
        input_1 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("de"),
            gene_tree & "z_1": list("abcde"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_1,
            ),
            rec_output,
            expected_cost=2,
            expected_labelings=[
                {
                    **input_1,
                    gene_tree & "2": list("abcde"),
                    gene_tree & "1": list("abcde"),
                },
            ],
        )

        # Test 2: Single optimal solution, altered root
        input_2 = {
            gene_tree & "x_1": list("abd"),
            gene_tree & "y_1": list("bde"),
            gene_tree & "z_1": list("abcde"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_2,
            ),
            rec_output,
            expected_cost=3,
            expected_labelings=[
                {
                    **input_2,
                    gene_tree & "2": list("abde"),
                    gene_tree & "1": list("abcde"),
                },
            ],
        )

        # Test 3: Multiple optimal solutions
        input_3 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("cd"),
            gene_tree & "z_1": list("be"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_3,
            ),
            rec_output,
            expected_cost=4,
            expected_labelings=[
                {
                    **input_3,
                    gene_tree & "2": list("cdab"),
                    gene_tree & "1": list("cdabe"),
                },
                {
                    **input_3,
                    gene_tree & "2": list("cdabe"),
                    gene_tree & "1": list("cdabe"),
                },
                {
                    **input_3,
                    gene_tree & "2": list("abecd"),
                    gene_tree & "1": list("abecd"),
                },
            ],
        )

        # Test 4: No solution
        input_4 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("cd"),
            gene_tree & "z_1": list("ba"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_4,
            ),
            rec_output,
            expected_cost=inf,
            expected_labelings=[],
        )

        # Test 5: Single gene
        input_5 = {
            gene_tree & "x_1": list("a"),
            gene_tree & "y_1": list("a"),
            gene_tree & "z_1": list("b"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_5,
            ),
            rec_output,
            expected_cost=2,
            expected_labelings=[
                {
                    **input_5,
                    gene_tree & "2": list("a"),
                    gene_tree & "1": list("ab"),
                },
                {
                    **input_5,
                    gene_tree & "2": list("a"),
                    gene_tree & "1": list("ba"),
                },
            ],
        )

    def test_duplications(self):
        gene_tree = Tree("((x_1,x_2)2,y_1)1;", format=1)
        species_tree = Tree("(X,Y)XY;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        leaf_gene_species = get_species_mapping(gene_tree, species_tree)

        rec_input = ReconciliationInput(
            gene_tree,
            species_lca,
            leaf_gene_species,
        )
        rec_output = reconcile_lca(rec_input)

        # Test 1: Single optimal solution, conserved root
        input_1 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "x_2": list("de"),
            gene_tree & "y_1": list("abcde"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_1,
            ),
            rec_output,
            expected_cost=2,
            expected_labelings=[
                {
                    **input_1,
                    gene_tree & "2": list("abcde"),
                    gene_tree & "1": list("abcde"),
                }
            ],
        )

        # Test 2: Single optimal solution, altered root
        input_2 = {
            gene_tree & "x_1": list("abd"),
            gene_tree & "x_2": list("bde"),
            gene_tree & "y_1": list("abcde"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_2,
            ),
            rec_output,
            expected_cost=3,
            expected_labelings=[
                {
                    **input_2,
                    gene_tree & "2": list("abde"),
                    gene_tree & "1": list("abcde"),
                }
            ],
        )

        # Test 3: Single optimal solution, conserved root (exchanged children)
        input_3 = {
            gene_tree & "x_1": list("de"),
            gene_tree & "x_2": list("ab"),
            gene_tree & "y_1": list("abcde"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_3,
            ),
            rec_output,
            expected_cost=2,
            expected_labelings=[
                {
                    **input_3,
                    gene_tree & "2": list("abcde"),
                    gene_tree & "1": list("abcde"),
                }
            ],
        )

        # Test 4: Single optimal solution, altered root (exchanged children)
        input_4 = {
            gene_tree & "x_1": list("bde"),
            gene_tree & "x_2": list("abd"),
            gene_tree & "y_1": list("abcde"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_4,
            ),
            rec_output,
            expected_cost=3,
            expected_labelings=[
                {
                    **input_4,
                    gene_tree & "2": list("abde"),
                    gene_tree & "1": list("abcde"),
                }
            ],
        )

        # Test 5: Multiple possible root labelings to examine
        input_5 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "x_2": list("cd"),
            gene_tree & "y_1": list("be"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_5,
            ),
            rec_output,
            expected_cost=3,
            expected_labelings=[
                {
                    **input_5,
                    gene_tree & "2": list("cdabe"),
                    gene_tree & "1": list("cdabe"),
                },
            ],
        )

        # Test 6: No solution
        input_6 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "x_2": list("cd"),
            gene_tree & "y_1": list("ba"),
        }

        self.assertLabelingEquals(
            SuperReconciliationInput(
                gene_tree,
                species_lca,
                leaf_gene_species,
                get_default_cost(),
                input_6,
            ),
            rec_output,
            expected_cost=inf,
            expected_labelings=[],
        )
