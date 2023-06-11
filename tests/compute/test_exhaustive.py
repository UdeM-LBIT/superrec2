from ete3 import Tree
from superrec2.utils.dynamic_programming import RetentionPolicy
from superrec2.utils.trees import LowestCommonAncestor
from superrec2.model.tree_mapping import get_species_mapping
from superrec2.model.reconciliation import (
    ReconciliationInput,
    ReconciliationOutput,
    NodeEvent,
)
from superrec2.compute.exhaustive import generate_all, reconcile_exhaustive


gene_tree = Tree("((x_1,x_2)2,(y_1,z_1)3)1;", format=1)
species_tree = Tree("(X,(Y,Z)YZ)XYZ;", format=1)
leaf_gene_species = get_species_mapping(gene_tree, species_tree)

species_lca = LowestCommonAncestor(species_tree)
rec_input = ReconciliationInput(
    gene_tree,
    species_lca,
    leaf_gene_species,
)


def test_generate_all():
    rec_outputs = list(generate_all(rec_input))

    # Check that all valid reconciliations are generated
    assert len(rec_outputs) == 16
    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "YZ",
                gene_tree & "1": species_tree & "XYZ",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "YZ",
                gene_tree & "1": species_tree & "X",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "YZ",
                gene_tree & "1": species_tree & "YZ",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "XYZ",
                gene_tree & "1": species_tree & "XYZ",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Y",
                gene_tree & "1": species_tree & "YZ",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Y",
                gene_tree & "1": species_tree & "XYZ",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Y",
                gene_tree & "1": species_tree & "X",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Y",
                gene_tree & "1": species_tree & "Y",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Z",
                gene_tree & "1": species_tree & "YZ",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Z",
                gene_tree & "1": species_tree & "XYZ",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Z",
                gene_tree & "1": species_tree & "X",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "X",
                gene_tree & "3": species_tree & "Z",
                gene_tree & "1": species_tree & "Z",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "XYZ",
                gene_tree & "3": species_tree & "YZ",
                gene_tree & "1": species_tree & "XYZ",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "XYZ",
                gene_tree & "3": species_tree & "XYZ",
                gene_tree & "1": species_tree & "XYZ",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "XYZ",
                gene_tree & "3": species_tree & "Y",
                gene_tree & "1": species_tree & "XYZ",
            },
        )
        in rec_outputs
    )

    assert (
        ReconciliationOutput(
            rec_input,
            {
                **leaf_gene_species,
                gene_tree & "2": species_tree & "XYZ",
                gene_tree & "3": species_tree & "Z",
                gene_tree & "1": species_tree & "XYZ",
            },
        )
        in rec_outputs
    )

    # Check that all the generated reconciliations are valid
    for rec in rec_outputs:
        for gene in gene_tree:
            assert rec.node_event(gene) != NodeEvent.INVALID


def test_exhaustive():
    all_outputs = list(generate_all(rec_input))
    min_cost = min(output.cost() for output in all_outputs)
    any_min = reconcile_exhaustive(rec_input, RetentionPolicy.ANY)
    all_min = reconcile_exhaustive(rec_input, RetentionPolicy.ALL)

    assert min(any_min) in all_outputs
    assert min(any_min) in all_min
    assert min(any_min).cost() == min_cost

    for a_min in all_min:
        assert a_min.cost() == min_cost
        assert a_min in all_outputs
