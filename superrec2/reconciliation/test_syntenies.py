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
                rec.tools.get_labeling_cost(
                    gene_tree, species_lca, rec_result, labeling
                ),
                expected_cost,
            )
            self.assertIn(labeling, output[1])

    def test_speciations(self):
        gene_tree = PhyloTree(
            "((x_1,y_1)2,z_1)1;",
            format=1,
            sp_naming_function=rec.tools.get_species_name,
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
            gene_tree,
            species_lca,
            rec_result,
            input_1,
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
            gene_tree & "y_1": list("bde"),
            gene_tree & "z_1": list("abcde"),
        }

        self.assertLabelingEquals(
            gene_tree,
            species_lca,
            rec_result,
            input_2,
            expected_cost=3,
            expected_labelings=[
                {
                    **input_2,
                    gene_tree & "2": list("abde"),
                    gene_tree & "1": list("abcde"),
                }
            ],
        )

        # Test 3: Multiple optimal solutions
        input_3 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("cd"),
            gene_tree & "z_1": list("be"),
        }

        self.assertLabelingEquals(
            gene_tree,
            species_lca,
            rec_result,
            input_3,
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
            ],
        )

        # Test 4: No solution
        input_4 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("cd"),
            gene_tree & "z_1": list("ba"),
        }

        self.assertLabelingEquals(
            gene_tree,
            species_lca,
            rec_result,
            input_4,
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
            gene_tree,
            species_lca,
            rec_result,
            input_5,
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
        gene_tree = PhyloTree(
            "((x_1,x_2)2,y_1)1;",
            format=1,
            sp_naming_function=rec.tools.get_species_name,
        )

        species_tree = PhyloTree("(X,Y)XY;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        rec_result = rec.compute.reconcile_lca(gene_tree, species_lca)

        # Test 1: Single optimal solution, conserved root
        input_1 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "x_2": list("de"),
            gene_tree & "y_1": list("abcde"),
        }

        self.assertLabelingEquals(
            gene_tree,
            species_lca,
            rec_result,
            input_1,
            expected_cost=1,
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
            gene_tree,
            species_lca,
            rec_result,
            input_2,
            expected_cost=2,
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
            gene_tree,
            species_lca,
            rec_result,
            input_3,
            expected_cost=1,
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
            gene_tree,
            species_lca,
            rec_result,
            input_4,
            expected_cost=2,
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
            gene_tree,
            species_lca,
            rec_result,
            input_5,
            expected_cost=2,
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
            gene_tree,
            species_lca,
            rec_result,
            input_6,
            expected_cost=inf,
            expected_labelings=[],
        )

    def test_transfer(self):
        gene_tree = PhyloTree(
            "((x_1,y_1)2,y_2)1;",
            format=1,
            sp_naming_function=rec.tools.get_species_name,
        )

        species_tree = PhyloTree("(X,Y)XY;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        rec_cost, rec_results = rec.compute.reconcile_thl(
            gene_tree,
            species_lca,
            costs={
                rec.tools.CostType.Duplication: 1,
                rec.tools.CostType.HorizontalGeneTransfer: 1,
                rec.tools.CostType.Loss: 1,
            },
        )

        self.assertEqual(rec_cost, 1)
        self.assertEqual(len(rec_results), 1)
        rec_result = rec_results[0]

        # Test 1: Single optimal solution, conserved root
        input_1 = {
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("de"),
            gene_tree & "y_2": list("abcde"),
        }

        self.assertLabelingEquals(
            gene_tree,
            species_lca,
            rec_result,
            input_1,
            expected_cost=1,
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
            gene_tree & "y_1": list("bde"),
            gene_tree & "y_2": list("abcde"),
        }

        self.assertLabelingEquals(
            gene_tree,
            species_lca,
            rec_result,
            input_2,
            expected_cost=2,
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
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("de"),
            gene_tree & "y_2": list("abcde"),
        }

        self.assertLabelingEquals(
            gene_tree,
            species_lca,
            rec_result,
            input_3,
            expected_cost=1,
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
            gene_tree & "x_1": list("abd"),
            gene_tree & "y_1": list("bde"),
            gene_tree & "y_2": list("abcde"),
        }

        self.assertLabelingEquals(
            gene_tree,
            species_lca,
            rec_result,
            input_4,
            expected_cost=2,
            expected_labelings=[
                {
                    **input_4,
                    gene_tree & "2": list("abde"),
                    gene_tree & "1": list("abcde"),
                }
            ],
        )

    def test_all(self):
        gene_tree = PhyloTree(
            """(
                ((x_1,z_1)3,(w_1,w_2)4)2,
                (
                    ((x_2,y_4)7,((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8)6,
                    (w_3,(z_3,(t_1,t_2)14)13)12
                )5
            )1;""",
            format=1,
            sp_naming_function=rec.tools.get_species_name,
        )

        species_tree = PhyloTree("(((X,Y)XY,Z)XYZ,(W,T)WT)XYZWT;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        rec_cost, rec_results = rec.compute.reconcile_thl(
            gene_tree,
            species_lca,
            costs={
                rec.tools.CostType.Duplication: 1,
                rec.tools.CostType.HorizontalGeneTransfer: 1,
                rec.tools.CostType.Loss: 1,
            },
        )

        self.assertEqual(rec_cost, 9)
        self.assertEqual(len(rec_results), 2)
        rec_result = rec_results[0]

        leaf_labeling = {
            gene_tree & "x_1": list("abcd"),
            gene_tree & "x_2": list("defg"),
            gene_tree & "x_3": list("cdef"),
            gene_tree & "y_1": list("cef"),
            gene_tree & "y_2": list("cde"),
            gene_tree & "y_3": list("cde"),
            gene_tree & "y_4": list("def"),
            gene_tree & "z_1": list("abcd"),
            gene_tree & "z_2": list("cef"),
            gene_tree & "z_3": list("defg"),
            gene_tree & "w_1": list("ab"),
            gene_tree & "w_2": list("abc"),
            gene_tree & "w_3": list("defg"),
            gene_tree & "t_1": list("def"),
            gene_tree & "t_2": list("defg"),
        }

        self.assertLabelingEquals(
            gene_tree,
            species_lca,
            rec_result,
            leaf_labeling,
            expected_cost=7,
            expected_labelings=[
                {
                    **leaf_labeling,
                    gene_tree & "1": list("abcdefg"),
                    gene_tree & "2": list("abcd"),
                    gene_tree & "3": list("abcd"),
                    gene_tree & "4": list("abc"),
                    gene_tree & "5": list("abcdefg"),
                    gene_tree & "6": list("cdefg"),
                    gene_tree & "7": list("defg"),
                    gene_tree & "8": list("cdef"),
                    gene_tree & "9": list("cdef"),
                    gene_tree & "10": list("cdef"),
                    gene_tree & "11": list("cde"),
                    gene_tree & "12": list("defg"),
                    gene_tree & "13": list("defg"),
                    gene_tree & "14": list("defg"),
                }
            ],
        )
