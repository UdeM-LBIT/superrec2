from ete3 import Tree
from infinity import inf
from superrec2.model.reconciliation import (
    SuperReconciliationInput,
    SuperReconciliationOutput,
    get_default_cost,
)
from superrec2.utils.trees import LowestCommonAncestor
from superrec2.utils.dynamic_programming import RetentionPolicy
from superrec2.model.tree_mapping import get_species_mapping
from superrec2.compute.reconciliation import reconcile_lca
from superrec2.compute.super_reconciliation import (
    sreconcile_base_spfs,
    sreconcile_extended_spfs,
)


def _assert_count_equal(list1, list2):
    assert len(list1) == len(list2)

    for element in list1:
        assert element in list2

    for element in list2:
        assert element in list1


def _assert_results(
    algo,
    srec_input,
    expected_cost,
    expected_results,
):
    any_l = algo(srec_input, RetentionPolicy.ANY)
    all_l = list(algo(srec_input, RetentionPolicy.ALL))
    _assert_count_equal(all_l, expected_results)

    if expected_results:
        assert any_l
        assert min(any_l).cost() == expected_cost
        assert min(any_l) in all_l

        for labeling in expected_results:
            assert labeling.cost() == expected_cost
    else:
        assert not any_l


def test_speciations():
    gene_tree = Tree("((x_1,y_1)2,z_1)1;", format=1)
    species_tree = Tree("((X,Y)XY,Z)XYZ;", format=1)
    species_lca = LowestCommonAncestor(species_tree)
    leaf_gene_species = get_species_mapping(gene_tree, species_tree)

    # Test 1: Single optimal solution, conserved root
    input_1 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("de"),
            gene_tree & "z_1": list("abcde"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_1,
        expected_cost=2,
        expected_results=[
            SuperReconciliationOutput(
                input_1,
                object_species=reconcile_lca(input_1).object_species,
                syntenies={
                    **input_1.leaf_syntenies,
                    gene_tree & "2": list("abcde"),
                    gene_tree & "1": list("abcde"),
                },
                ordered=True,
            ),
        ],
    )

    # Test 2: Single optimal solution, altered root
    input_2 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("abd"),
            gene_tree & "y_1": list("bde"),
            gene_tree & "z_1": list("abcde"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_2,
        expected_cost=3,
        expected_results=[
            SuperReconciliationOutput(
                input_2,
                object_species=reconcile_lca(input_2).object_species,
                syntenies={
                    **input_2.leaf_syntenies,
                    gene_tree & "2": list("abde"),
                    gene_tree & "1": list("abcde"),
                },
                ordered=True,
            ),
        ],
    )

    # Test 3: Multiple optimal solutions
    input_3 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("cd"),
            gene_tree & "z_1": list("be"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_3,
        expected_cost=4,
        expected_results=[
            SuperReconciliationOutput(
                input_3,
                object_species=reconcile_lca(input_3).object_species,
                syntenies={
                    **input_3.leaf_syntenies,
                    gene_tree & "2": list("cdab"),
                    gene_tree & "1": list("cdabe"),
                },
                ordered=True,
            ),
            SuperReconciliationOutput(
                input_3,
                object_species=reconcile_lca(input_3).object_species,
                syntenies={
                    **input_3.leaf_syntenies,
                    gene_tree & "2": list("cdabe"),
                    gene_tree & "1": list("cdabe"),
                },
                ordered=True,
            ),
            SuperReconciliationOutput(
                input_3,
                object_species=reconcile_lca(input_3).object_species,
                syntenies={
                    **input_3.leaf_syntenies,
                    gene_tree & "2": list("abecd"),
                    gene_tree & "1": list("abecd"),
                },
                ordered=True,
            ),
        ],
    )

    # Test 4: No solution
    input_4 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("ab"),
            gene_tree & "y_1": list("cd"),
            gene_tree & "z_1": list("ba"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_4,
        expected_cost=inf,
        expected_results=[],
    )

    # Test 5: Single gene
    input_5 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("a"),
            gene_tree & "y_1": list("a"),
            gene_tree & "z_1": list("b"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_5,
        expected_cost=2,
        expected_results=[
            SuperReconciliationOutput(
                input_5,
                object_species=reconcile_lca(input_5).object_species,
                syntenies={
                    **input_5.leaf_syntenies,
                    gene_tree & "2": list("a"),
                    gene_tree & "1": list("ab"),
                },
                ordered=True,
            ),
            SuperReconciliationOutput(
                input_5,
                object_species=reconcile_lca(input_5).object_species,
                syntenies={
                    **input_5.leaf_syntenies,
                    gene_tree & "2": list("a"),
                    gene_tree & "1": list("ba"),
                },
                ordered=True,
            ),
        ],
    )


def test_duplications():
    gene_tree = Tree("((x_1,x_2)2,y_1)1;", format=1)
    species_tree = Tree("(X,Y)XY;", format=1)
    species_lca = LowestCommonAncestor(species_tree)
    leaf_gene_species = get_species_mapping(gene_tree, species_tree)

    # Test 1: Single optimal solution, conserved root
    input_1 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("ab"),
            gene_tree & "x_2": list("de"),
            gene_tree & "y_1": list("abcde"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_1,
        expected_cost=2,
        expected_results=[
            SuperReconciliationOutput(
                input_1,
                object_species=reconcile_lca(input_1).object_species,
                syntenies={
                    **input_1.leaf_syntenies,
                    gene_tree & "2": list("abcde"),
                    gene_tree & "1": list("abcde"),
                },
                ordered=True,
            ),
        ],
    )

    # Test 2: Single optimal solution, altered root
    input_2 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("abd"),
            gene_tree & "x_2": list("bde"),
            gene_tree & "y_1": list("abcde"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_2,
        expected_cost=3,
        expected_results=[
            SuperReconciliationOutput(
                input_2,
                object_species=reconcile_lca(input_2).object_species,
                syntenies={
                    **input_2.leaf_syntenies,
                    gene_tree & "2": list("abde"),
                    gene_tree & "1": list("abcde"),
                },
                ordered=True,
            ),
        ],
    )

    # Test 3: Single optimal solution, conserved root (exchanged children)
    input_3 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("de"),
            gene_tree & "x_2": list("ab"),
            gene_tree & "y_1": list("abcde"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_3,
        expected_cost=2,
        expected_results=[
            SuperReconciliationOutput(
                input_3,
                object_species=reconcile_lca(input_3).object_species,
                syntenies={
                    **input_3.leaf_syntenies,
                    gene_tree & "2": list("abcde"),
                    gene_tree & "1": list("abcde"),
                },
                ordered=True,
            ),
        ],
    )

    # Test 4: Single optimal solution, altered root (exchanged children)
    input_4 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("bde"),
            gene_tree & "x_2": list("abd"),
            gene_tree & "y_1": list("abcde"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_4,
        expected_cost=3,
        expected_results=[
            SuperReconciliationOutput(
                input_4,
                object_species=reconcile_lca(input_4).object_species,
                syntenies={
                    **input_4.leaf_syntenies,
                    gene_tree & "2": list("abde"),
                    gene_tree & "1": list("abcde"),
                },
                ordered=True,
            ),
        ],
    )

    # Test 5: Multiple possible root labelings to examine
    input_5 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("ab"),
            gene_tree & "x_2": list("cd"),
            gene_tree & "y_1": list("be"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_5,
        expected_cost=3,
        expected_results=[
            SuperReconciliationOutput(
                input_5,
                object_species=reconcile_lca(input_5).object_species,
                syntenies={
                    **input_5.leaf_syntenies,
                    gene_tree & "2": list("cdabe"),
                    gene_tree & "1": list("cdabe"),
                },
                ordered=True,
            ),
        ],
    )

    # Test 6: No solution
    input_6 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("ab"),
            gene_tree & "x_2": list("cd"),
            gene_tree & "y_1": list("ba"),
        },
    )

    _assert_results(
        sreconcile_base_spfs,
        input_6,
        expected_cost=inf,
        expected_results=[],
    )


def test_transfers():
    gene_tree = Tree("((x_1,y_1)1,(x_2,y_2)2)3;", format=1)
    species_tree = Tree("(X,Y)XY;", format=1)
    species_lca = LowestCommonAncestor(species_tree)
    leaf_gene_species = get_species_mapping(gene_tree, species_tree)

    input_1 = SuperReconciliationInput(
        gene_tree,
        species_lca,
        leaf_gene_species,
        get_default_cost(),
        leaf_syntenies={
            gene_tree & "x_1": list("a"),
            gene_tree & "x_2": list("a"),
            gene_tree & "y_1": list("b"),
            gene_tree & "y_2": list("b"),
        },
    )

    _assert_results(
        sreconcile_extended_spfs,
        input_1,
        expected_cost=4,
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
                    gene_tree & "1": list("ab"),
                    gene_tree & "2": list("ab"),
                    gene_tree & "3": list("ab"),
                },
                ordered=True,
            ),
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
                    gene_tree & "1": list("ba"),
                    gene_tree & "2": list("ba"),
                    gene_tree & "3": list("ba"),
                },
                ordered=True,
            ),
            SuperReconciliationOutput(
                input_1,
                object_species={
                    **leaf_gene_species,
                    gene_tree & "1": species_tree & "Y",
                    gene_tree & "2": species_tree & "X",
                    gene_tree & "3": species_tree & "XY",
                },
                syntenies={
                    **input_1.leaf_syntenies,
                    gene_tree & "1": list("ab"),
                    gene_tree & "2": list("ab"),
                    gene_tree & "3": list("ab"),
                },
                ordered=True,
            ),
            SuperReconciliationOutput(
                input_1,
                object_species={
                    **leaf_gene_species,
                    gene_tree & "1": species_tree & "Y",
                    gene_tree & "2": species_tree & "X",
                    gene_tree & "3": species_tree & "XY",
                },
                syntenies={
                    **input_1.leaf_syntenies,
                    gene_tree & "1": list("ba"),
                    gene_tree & "2": list("ba"),
                    gene_tree & "3": list("ba"),
                },
                ordered=True,
            ),
        ],
    )
