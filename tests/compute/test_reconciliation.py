from infinity import inf
from ete3 import Tree
from superrec2.utils.trees import LowestCommonAncestor
from superrec2.utils.dynamic_programming import RetentionPolicy
from superrec2.model.tree_mapping import get_species_mapping
from superrec2.model.reconciliation import (
    ReconciliationInput,
    ReconciliationOutput,
    NodeEvent,
    EdgeEvent,
)
from superrec2.compute.exhaustive import generate_all
from superrec2.compute.reconciliation import reconcile_lca, reconcile_thl


gene_tree = Tree(
    "((x_1,(x_2,(y_1,z_1)5)4)3,(y_2,z_2)2)1;",
    format=1,
)

species_tree = Tree("(X,(Y,Z)YZ)XYZ;", format=1)
species_lca = LowestCommonAncestor(species_tree)
leaf_gene_species = get_species_mapping(gene_tree, species_tree)

rec_input = ReconciliationInput(
    gene_tree,
    species_lca,
    leaf_gene_species,
)

all_outputs = list(generate_all(rec_input))


def test_output_count():
    assert len(all_outputs) == 199


def test_reconcile_lca():
    rec_input.costs[NodeEvent.DUPLICATION] = 1
    rec_input.costs[NodeEvent.HORIZONTAL_TRANSFER] = inf
    rec_input.costs[EdgeEvent.FULL_LOSS] = 1

    result = reconcile_lca(rec_input)

    # Check that the expected result is returned
    assert result.object_species == {
        **leaf_gene_species,
        gene_tree & "1": species_tree & "XYZ",
        gene_tree & "2": species_tree & "YZ",
        gene_tree & "3": species_tree & "XYZ",
        gene_tree & "4": species_tree & "XYZ",
        gene_tree & "5": species_tree & "YZ",
    }

    assert result.cost() == 4
    assert result in all_outputs

    # Check optimality
    for possible_rec in all_outputs:
        cost = possible_rec.cost()
        assert cost >= result.cost()

        if cost == result.cost():
            assert possible_rec == result


def test_reconcile_thl():
    rec_input.costs[NodeEvent.DUPLICATION] = 1
    rec_input.costs[NodeEvent.HORIZONTAL_TRANSFER] = 1
    rec_input.costs[EdgeEvent.FULL_LOSS] = 1

    any_result = reconcile_thl(rec_input, RetentionPolicy.ANY)
    results = list(reconcile_thl(rec_input, RetentionPolicy.ALL))

    # Check that any is part of all results
    assert min(any_result) in results

    # Check that all expected results are returned
    assert len(results) == 2
    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "1": species_tree & "XYZ",
                gene_tree & "2": species_tree & "YZ",
                gene_tree & "3": species_tree & "X",
                gene_tree & "4": species_tree & "X",
                gene_tree & "5": species_tree & "YZ",
            },
        )
        in results
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "1": species_tree & "XYZ",
                gene_tree & "2": species_tree & "YZ",
                gene_tree & "3": species_tree & "X",
                gene_tree & "4": species_tree & "YZ",
                gene_tree & "5": species_tree & "YZ",
            },
        )
        in results
    )

    # Check that all results have the same expected cost
    for result in results:
        assert result.cost() == 2

    for result in results:
        assert result in all_outputs

    # Check optimality
    for possible_rec in all_outputs:
        cost = possible_rec.cost()
        assert cost >= 2

        if cost == 2:
            assert possible_rec in results
