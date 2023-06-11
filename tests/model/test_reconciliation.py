from ete3 import Tree
from superrec2.utils.trees import LowestCommonAncestor
from superrec2.model.reconciliation import (
    ReconciliationInput,
    SuperReconciliationInput,
    ReconciliationOutput,
    SuperReconciliationOutput,
    EdgeEvent,
    NodeEvent,
)


gene_tree = Tree(
    """
    (
        ((x_1,z_1)3,(w_1,w_2)4)2,
        (
            ((x_2,y_4)7, ((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8)6,
            (w_3,(z_3,(t_1,t_2)14)13)12
        )5
    )1;
""",
    format=1,
)

species_tree = Tree("(((X,Y)XY,Z)XYZ,(W,T)WT)XYZWT;", format=1)
species_lca = LowestCommonAncestor(species_tree)

leaf_object_species = {
    gene_tree & "x_1": species_tree & "X",
    gene_tree & "x_2": species_tree & "X",
    gene_tree & "x_3": species_tree & "X",
    gene_tree & "y_1": species_tree & "Y",
    gene_tree & "y_2": species_tree & "Y",
    gene_tree & "y_3": species_tree & "Y",
    gene_tree & "y_4": species_tree & "Y",
    gene_tree & "z_1": species_tree & "Z",
    gene_tree & "z_2": species_tree & "Z",
    gene_tree & "z_3": species_tree & "Z",
    gene_tree & "w_1": species_tree & "W",
    gene_tree & "w_2": species_tree & "W",
    gene_tree & "w_3": species_tree & "W",
    gene_tree & "t_1": species_tree & "T",
    gene_tree & "t_2": species_tree & "T",
}

object_species = {
    **leaf_object_species,
    gene_tree & "1": species_tree & "XYZWT",
    gene_tree & "2": species_tree & "XYZ",
    gene_tree & "3": species_tree & "XYZ",
    gene_tree & "4": species_tree & "W",
    gene_tree & "5": species_tree & "XYZWT",
    gene_tree & "6": species_tree & "XYZ",
    gene_tree & "7": species_tree & "XY",
    gene_tree & "8": species_tree & "XYZ",
    gene_tree & "9": species_tree & "XY",
    gene_tree & "10": species_tree & "Y",
    gene_tree & "11": species_tree & "Y",
    gene_tree & "12": species_tree & "WT",
    gene_tree & "13": species_tree & "T",
    gene_tree & "14": species_tree & "T",
}

leaf_syntenies = {
    gene_tree & "x_1": "abcd",
    gene_tree & "z_1": "abcd",
    gene_tree & "w_1": "ab",
    gene_tree & "w_2": "abc",
    gene_tree & "x_2": "defg",
    gene_tree & "y_4": "def",
    gene_tree & "x_3": "cdef",
    gene_tree & "y_1": "ce",
    gene_tree & "y_2": "cde",
    gene_tree & "y_3": "cde",
    gene_tree & "z_2": "cef",
    gene_tree & "w_3": "defg",
    gene_tree & "z_3": "defg",
    gene_tree & "t_1": "def",
    gene_tree & "t_2": "defg",
}

syntenies = {
    **leaf_syntenies,
    gene_tree & "1": "abcdefg",
    gene_tree & "2": "abcd",
    gene_tree & "3": "abcd",
    gene_tree & "4": "abc",
    gene_tree & "5": "abcdefg",
    gene_tree & "6": "abcdefg",
    gene_tree & "7": "defg",
    gene_tree & "8": "cdef",
    gene_tree & "9": "cdef",
    gene_tree & "10": "cdef",
    gene_tree & "11": "cde",
    gene_tree & "12": "defg",
    gene_tree & "13": "defg",
    gene_tree & "14": "defg",
}

rec_input = ReconciliationInput(
    gene_tree,
    species_lca,
    leaf_object_species,
)

srec_input = SuperReconciliationInput(
    gene_tree,
    species_lca,
    leaf_object_species,
    leaf_syntenies,
)

rec_output = ReconciliationOutput(
    rec_input,
    object_species,
)

srec_output = SuperReconciliationOutput(
    rec_input, object_species, syntenies, ordered=True
)

usrec_output = SuperReconciliationOutput(
    rec_input, object_species, syntenies, ordered=False
)


def test_node_event():
    expected_events = {
        "1": NodeEvent.DUPLICATION,
        "2": NodeEvent.HORIZONTAL_TRANSFER,
        "3": NodeEvent.SPECIATION,
        "4": NodeEvent.DUPLICATION,
        "5": NodeEvent.SPECIATION,
        "6": NodeEvent.DUPLICATION,
        "7": NodeEvent.SPECIATION,
        "8": NodeEvent.SPECIATION,
        "9": NodeEvent.SPECIATION,
        "10": NodeEvent.DUPLICATION,
        "11": NodeEvent.DUPLICATION,
        "12": NodeEvent.SPECIATION,
        "13": NodeEvent.HORIZONTAL_TRANSFER,
        "14": NodeEvent.DUPLICATION,
    }

    for name, event in expected_events.items():
        assert rec_output.node_event(gene_tree & name) == event


def test_labeling_cost():
    assert srec_output.labeling_cost() == 6
    assert usrec_output.labeling_cost() == 5


def test_reconciliation_cost():
    expected_costs = {
        (1, 0, 0): 6,
        (0, 1, 0): 2,
        (0, 0, 1): 3,
        (1, 1, 1): 11,
    }

    for (dup, hgt, loss), value in expected_costs.items():
        rec_input.costs[NodeEvent.DUPLICATION] = dup
        rec_input.costs[NodeEvent.HORIZONTAL_TRANSFER] = hgt
        rec_input.costs[EdgeEvent.FULL_LOSS] = loss
        assert rec_output.cost() == value
