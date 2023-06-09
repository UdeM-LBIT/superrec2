import unittest
from ete3 import Tree
from infinity import inf
from superrec2.model.reconciliation import (
    ReconciliationInput,
    SuperReconciliationInput,
    SuperReconciliationOutput,
    get_default_cost,
)
from superrec2.utils.trees import LowestCommonAncestor
from superrec2.utils.dynamic_programming import RetentionPolicy
from superrec2.model.tree_mapping import get_species_mapping
from superrec2.compute.reconciliation import reconcile_lca
from superrec2.compute.unordered_super_reconciliation import (
    _compute_gain_sets,
    _compute_lca_sets,
    usreconcile_base_uspfs,
    usreconcile_extended_uspfs,
)


class TestComputeUnorderedSuperReconciliation(unittest.TestCase):
    def test_gain_lca_sets(self):
        gene_tree = Tree("((((a,b)1,c)2,d)3,(e,f)4)5;", format=1)
        leaf_syntenies = {
            gene_tree & "a": "ab",
            gene_tree & "b": "ac",
            gene_tree & "c": "ad",
            gene_tree & "d": "bd",
            gene_tree & "e": "cde",
            gene_tree & "f": "ef",
        }
        s_input = SuperReconciliationInput(
            gene_tree,
            Tree(),
            {},
            get_default_cost(),
            leaf_syntenies,
        )

        gain_sets = _compute_gain_sets(s_input)
        self.assertEqual(
            gain_sets,
            {
                gene_tree & "a": set(),
                gene_tree & "b": set(),
                gene_tree & "c": set(),
                gene_tree & "d": set(),
                gene_tree & "e": set(),
                gene_tree & "f": set("f"),
                gene_tree & "1": set(),
                gene_tree & "2": set("a"),
                gene_tree & "3": set("b"),
                gene_tree & "4": set("e"),
                gene_tree & "5": set("cd"),
            },
        )

        lca_sets = _compute_lca_sets(s_input, gain_sets)
        self.assertEqual(
            lca_sets,
            {
                gene_tree & "a": set("ab"),
                gene_tree & "b": set("ac"),
                gene_tree & "c": set("ad"),
                gene_tree & "d": set("bd"),
                gene_tree & "e": set("cde"),
                gene_tree & "f": set("ef"),
                gene_tree & "1": set("abc"),
                gene_tree & "2": set("abcd"),
                gene_tree & "3": set("bcd"),
                gene_tree & "4": set("cde"),
                gene_tree & "5": set("cd"),
            },
        )

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
            expected_cost=0,
            expected_results=[
                SuperReconciliationOutput(
                    input_1,
                    object_species=reconcile_lca(input_1).object_species,
                    syntenies={
                        **input_1.leaf_syntenies,
                        gene_tree & "2": list("a"),
                        gene_tree & "1": list("a"),
                    },
                    ordered=False,
                ),
            ],
        )

        # Test 2: Gains and losses
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
            expected_cost=2,
            expected_results=[
                SuperReconciliationOutput(
                    input_2,
                    object_species=reconcile_lca(input_2).object_species,
                    syntenies={
                        **input_2.leaf_syntenies,
                        gene_tree & "2": list("abc"),
                        gene_tree & "1": list("bc"),
                    },
                    ordered=False,
                ),
            ],
        )

        # Test 3: Gains and losses with inheritance
        gene_tree = Tree("((((x_1,y_1)1,z_1)2,t_1)3,(w_1,v_1)4)5;", format=1)
        species_tree = Tree("((((X,Y)XY,Z)XYZ,T)XYZT,(W,V)WV)XYZTVW;", format=1)
        species_lca = LowestCommonAncestor(species_tree)
        leaf_gene_species = get_species_mapping(gene_tree, species_tree)

        input_3 = SuperReconciliationInput(
            gene_tree,
            species_lca,
            leaf_gene_species,
            get_default_cost(),
            leaf_syntenies={
                gene_tree & "x_1": list("abx"),
                gene_tree & "y_1": list("acx"),
                gene_tree & "z_1": list("ad"),
                gene_tree & "t_1": list("bd"),
                gene_tree & "w_1": list("cde"),
                gene_tree & "v_1": list("ef"),
            },
        )

        self.assertResults(
            usreconcile_base_uspfs,
            input_3,
            expected_cost=5,
            expected_results=[
                SuperReconciliationOutput(
                    input_3,
                    object_species=reconcile_lca(input_3).object_species,
                    syntenies={
                        **input_3.leaf_syntenies,
                        gene_tree & "1": list("abcdx"),
                        gene_tree & "2": list("abcd"),
                        gene_tree & "3": list("bcd"),
                        gene_tree & "4": list("cde"),
                        gene_tree & "5": list("cd"),
                    },
                    ordered=False,
                ),
            ],
        )

    def test_duplications(self):
        gene_tree = Tree("(((x_1,x_2)3,y_1)2,z_1)1;", format=1)
        species_tree = Tree("((X,Y)XY,Z);", format=1)
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
                gene_tree & "z_1": list("ab"),
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
                        gene_tree & "3": list("a"),
                        gene_tree & "2": list("ab"),
                        gene_tree & "1": list("ab"),
                    },
                    ordered=False,
                ),
                SuperReconciliationOutput(
                    input_1,
                    object_species=reconcile_lca(input_1).object_species,
                    syntenies={
                        **input_1.leaf_syntenies,
                        gene_tree & "3": list("ab"),
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
                gene_tree & "z_1": list("abcd"),
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
                        gene_tree & "3": list("abcd"),
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
