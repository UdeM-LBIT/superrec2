import re
from typing import Dict, Mapping
from ete3 import Tree, TreeNode
from ..utils.escape import escape, unescape


TreeMapping = Mapping[TreeNode, TreeNode]

# Set of characters that must be escaped from names
ESCAPE_CHARS = r"\,"
MAPPING_SPLIT_REGEX = re.compile(r"(?<!\\),")


def parse_tree_mapping(
    from_tree: Tree, to_tree: Tree, data: str
) -> TreeMapping:
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
                lambda x: unescape(x.strip()), pair.split(":")
            )

            if not from_node.startswith("#"):
                result[from_tree & from_node] = to_tree & to_node

    return result


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


def serialize_tree_mapping(mapping: TreeMapping) -> str:
    """
    Serialize a mapping between two trees.

    :param mapping: mapping to serialize
    :returns: serialized representation
    """
    return ",".join(
        escape(from_name, ESCAPE_CHARS) + ":" + escape(to_name, ESCAPE_CHARS)
        for from_name, to_name in sorted(
            (from_node.name, to_node.name)
            for from_node, to_node in mapping.items()
        )
    )
