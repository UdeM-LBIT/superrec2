import re
from typing import Dict, Mapping
from ete3 import PhyloTree, Tree, TreeNode
from ..utils.escape import escape, unescape


TreeMapping = Mapping[TreeNode, TreeNode]

# Set of characters that must be escaped from names
ESCAPE_CHARS = r"\,"
MAPPING_SPLIT_REGEX = re.compile(r"(?<!\\),")


def parse_tree_mapping(from_tree: Tree, to_tree: Tree, data: str) \
-> TreeMapping:
    """
    Parse a string representation of a mapping between two trees.

    :param from_tree: first tree
    :param to_tree: second tree
    :param data: string to parse
    :returns: parsed mapping
    """
    result: Dict[TreeNode, TreeNode] = {}

    for pair in MAPPING_SPLIT_REGEX.split(data):
        if pair.strip():
            from_node, to_node = map(
                lambda x: unescape(x.strip()),
                pair.split(":")
            )

            if not from_node.startswith("#"):
                result[from_tree & from_node] = to_tree & to_node

    return result


def get_species_mapping(tree: PhyloTree, species_tree: PhyloTree) \
-> TreeMapping:
    """
    Extract a mapping of a tree onto a species tree from node names.

    Only nodes that actually have species mapping information get included.

    :param tree: tree containing mapping information
    :param species_tree: target species tree
    :returns: extracted mapping
    """
    return {
        node: species_tree & node.species
        for node in tree
        if node.species
    }


def serialize_tree_mapping(mapping: TreeMapping) -> str:
    """
    Serialize a mapping between two trees.

    :param mapping: mapping to serialize
    :returns: serialized representation
    """
    return ",".join(
        escape(from_node.name, ESCAPE_CHARS)
        + ":" + escape(to_node.name, ESCAPE_CHARS)
        for from_node, to_node in mapping.items()
    )
