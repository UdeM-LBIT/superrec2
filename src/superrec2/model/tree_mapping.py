"""Representation and parsing of mappings between tree nodes."""
from typing import Dict, Mapping
from ete3 import Tree, TreeNode


TreeMapping = Mapping[TreeNode, TreeNode]


def parse_tree_mapping(
    from_tree: Tree, to_tree: Tree, data: Mapping[str, str]
) -> TreeMapping:
    """
    Convert a plain mapping of tree node names to a mapping between two trees.

    :param from_tree: first tree
    :param to_tree: second tree
    :param data: plain mapping to convert from
    :returns: parsed mapping
    """
    return {
        from_tree & from_node: to_tree & to_node for from_node, to_node in data.items()
    }


def get_species_mapping(tree: Tree, species_tree: Tree) -> TreeMapping:
    """
    Extract a mapping of a tree onto a species tree from node names.

    Each node whose name matches the `<species>_<suffix>` format, where
    `<species>` is the name of a node in :param:`species_tree` and `_<suffix>`
    is an optional arbitrary suffix, gets mapped to that species. Matching
    is case-insensitive. If the species name contains underscores, the first
    matching prefix that is followed by an underscore is considered. Other
    nodes are left unmapped.

    :param tree: tree of objects
    :param species_tree: target species tree
    :returns: extracted mapping
    """
    result = {}
    species_ignorecase = {}

    for node in species_tree:
        if node.name:
            species_ignorecase[node.name.lower()] = node

    for node in tree:
        parts = node.name.split("_")

        for i in range(1, len(parts)):
            prefix = "_".join(parts[:i]).lower()

            if prefix in species_ignorecase:
                result[node] = species_ignorecase[prefix]
                break

    return result


def serialize_tree_mapping(mapping: TreeMapping) -> Dict[str, str]:
    """
    Convert a mapping between two trees to a plain dictionary.

    :param mapping: mapping to serialize
    :returns: representation that can be used to serialize the mapping
    """
    return {from_node.name: to_node.name for from_node, to_node in mapping.items()}
