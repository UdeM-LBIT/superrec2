from ete3 import Tree
from superrec2.model.tree_mapping import (
    parse_tree_mapping,
    get_species_mapping,
    serialize_tree_mapping,
)


gene_tree = Tree(
    """
    (
        ((x_1,z_1)3,(w_1,w_2)4)2,
        (
            ((x_2,y_4)7,((x_3,(y_1,(y_2,y_3)11)10)9,z_2)8)6,
            (w_3,(z_3,(t_1,t_2)14)13)12
        )5
    )1;
""",
    format=1,
)
species_tree = Tree("(((X,Y)XY,Z)XYZ,(W,T)WT)XYZWT;", format=1)


def test_parse_tree_mapping():
    assert parse_tree_mapping(
        gene_tree,
        species_tree,
        {
            "1": "XYZWT",
            "2": "XYZ",
            "3": "XYZ",
            "4": "W",
            "5": "XYZWT",
            "6": "XYZ",
            "7": "XY",
            "8": "XYZ",
            "9": "XY",
            "10": "Y",
            "11": "Y",
            "12": "WT",
            "13": "T",
            "14": "T",
        },
    ) == {
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


def test_get_species_mapping():
    assert get_species_mapping(gene_tree, species_tree) == {
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


def test_serialize_tree_mapping():
    assert serialize_tree_mapping(
        {
            gene_tree & "x_1": species_tree & "X",
            gene_tree & "2": species_tree & "XYZ",
            gene_tree & "x_2": species_tree & "X",
            gene_tree & "3": species_tree & "XYZ",
        }
    ) == {
        "2": "XYZ",
        "3": "XYZ",
        "x_1": "X",
        "x_2": "X",
    }
