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
from ..utils.dynamic_programming import RetentionPolicy
from ..model.tree_mapping import get_species_mapping
from .reconciliation import reconcile_lca
from .unordered_super_reconciliation import (
    usreconcile_base_uspfs,
    usreconcile_extended_uspfs,
)


class TestComputeUnorderedSuperReconciliation(unittest.TestCase):
    def assertResults(
        self,
        algo,
        srec_input,
        expected_cost,
        expected_results,
    ):
        any_l = algo(srec_input, RetentionPolicy.ANY)
        all_l = list(algo(srec_input, RetentionPolicy.ALL))
        self.assertCountEqual(all_l, expected_results)

        if expected_results:
            self.assertTrue(any_l)
            self.assertEqual(min(any_l).cost(), expected_cost)
            self.assertIn(min(any_l), all_l)

            for labeling in expected_results:
                self.assertEqual(labeling.cost(), expected_cost)
        else:
            self.assertFalse(any_l)

    def test_speciations(self):
        gene_tree = Tree("((x_1,y_1)2,z_1)1;", format=1)
        species_tree = Tree("((X,Y)XY,Z)XYZ;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        leaf_gene_species = get_species_mapping(gene_tree, species_tree)

        # Test 1: Only LCA syntenies
        input_1 = SuperReconciliationInput(
            gene_tree,
            species_lca,
            leaf_gene_species,
            get_default_cost(),
            leaf_syntenies={
                gene_tree & "x_1": list("a"),
                gene_tree & "y_1": list("a"),
                gene_tree & "z_1": list("ab"),
            },
        )

        self.assertResults(
            usreconcile_base_uspfs,
            input_1,
            expected_cost=1,
            expected_results=[
                SuperReconciliationOutput(
                    input_1,
                    object_species=reconcile_lca(input_1).object_species,
                    syntenies={
                        **input_1.leaf_syntenies,
                        gene_tree & "2": list("a"),
                        gene_tree & "1": list("ab"),
                    },
                    ordered=False,
                ),
            ],
        )

        # Test 2: Inherit synteny
        input_2 = SuperReconciliationInput(
            gene_tree,
            species_lca,
            leaf_gene_species,
            get_default_cost(),
            leaf_syntenies={
                gene_tree & "x_1": list("ab"),
                gene_tree & "y_1": list("ac"),
                gene_tree & "z_1": list("bcd"),
            },
        )

        self.assertResults(
            usreconcile_base_uspfs,
            input_2,
            expected_cost=3,
            expected_results=[
                SuperReconciliationOutput(
                    input_2,
                    object_species=reconcile_lca(input_2).object_species,
                    syntenies={
                        **input_2.leaf_syntenies,
                        gene_tree & "2": list("abcd"),
                        gene_tree & "1": list("abcd"),
                    },
                    ordered=False,
                ),
            ],
        )

    def test_duplications(self):
        gene_tree = Tree("((x_1,x_2)2,y_1)1;", format=1)
        species_tree = Tree("(X,Y)XY;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        leaf_gene_species = get_species_mapping(gene_tree, species_tree)

        # Test 1: LCA or INH syntenies
        input_1 = SuperReconciliationInput(
            gene_tree,
            species_lca,
            leaf_gene_species,
            get_default_cost(),
            leaf_syntenies={
                gene_tree & "x_1": list("a"),
                gene_tree & "x_2": list("a"),
                gene_tree & "y_1": list("ab"),
            },
        )

        self.assertResults(
            usreconcile_base_uspfs,
            input_1,
            expected_cost=2,
            expected_results=[
                SuperReconciliationOutput(
                    input_1,
                    object_species=reconcile_lca(input_1).object_species,
                    syntenies={
                        **input_1.leaf_syntenies,
                        gene_tree & "2": list("a"),
                        gene_tree & "1": list("ab"),
                    },
                    ordered=False,
                ),
                SuperReconciliationOutput(
                    input_1,
                    object_species=reconcile_lca(input_1).object_species,
                    syntenies={
                        **input_1.leaf_syntenies,
                        gene_tree & "2": list("ab"),
                        gene_tree & "1": list("ab"),
                    },
                    ordered=False,
                ),
            ],
        )

        # Test 2: Inherit synteny
        input_2 = SuperReconciliationInput(
            gene_tree,
            species_lca,
            leaf_gene_species,
            get_default_cost(),
            leaf_syntenies={
                gene_tree & "x_1": list("ab"),
                gene_tree & "x_2": list("ac"),
                gene_tree & "y_1": list("bcd"),
            },
        )

        self.assertResults(
            usreconcile_base_uspfs,
            input_2,
            expected_cost=3,
            expected_results=[
                SuperReconciliationOutput(
                    input_2,
                    object_species=reconcile_lca(input_2).object_species,
                    syntenies={
                        **input_2.leaf_syntenies,
                        gene_tree & "2": list("abcd"),
                        gene_tree & "1": list("abcd"),
                    },
                    ordered=False,
                ),
            ],
        )

    def test_transfers(self):
        gene_tree = Tree("((x_1,y_1)1,(x_2,y_2)2)3;", format=1)
        species_tree = Tree("(X,Y)XY;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        leaf_gene_species = get_species_mapping(gene_tree, species_tree)

        # Test 1: costlier reconciliation to save segmental losses
        input_1 = SuperReconciliationInput(
            gene_tree,
            species_lca,
            leaf_gene_species,
            get_default_cost(),
            leaf_syntenies={
                gene_tree & "x_1": list("abc"),
                gene_tree & "y_1": list("ab"),
                gene_tree & "x_2": list("bc"),
                gene_tree & "y_2": list("abc"),
            },
        )

        self.assertResults(
            usreconcile_base_uspfs,
            input_1,
            expected_cost=3,
            expected_results=[
                SuperReconciliationOutput(
                    input_1,
                    object_species=reconcile_lca(input_1).object_species,
                    syntenies={
                        **input_1.leaf_syntenies,
                        gene_tree & "1": list("abc"),
                        gene_tree & "2": list("abc"),
                        gene_tree & "3": list("abc"),
                    },
                    ordered=False,
                ),
            ],
        )

        self.assertResults(
            usreconcile_extended_uspfs,
            input_1,
            expected_cost=2,
            expected_results=[
                SuperReconciliationOutput(
                    input_1,
                    object_species={
                        **leaf_gene_species,
                        gene_tree & "1": species_tree & "X",
                        gene_tree & "2": species_tree & "Y",
                        gene_tree & "3": species_tree & "XY",
                    },
                    syntenies={
                        **input_1.leaf_syntenies,
                        gene_tree & "1": list("abc"),
                        gene_tree & "2": list("abc"),
                        gene_tree & "3": list("abc"),
                    },
                    ordered=False,
                ),
            ],
        )
