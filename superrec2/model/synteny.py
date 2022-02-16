import re
from typing import Dict, Iterable, List, Mapping, Sequence, Set
from ete3 import Tree, TreeNode


GeneFamily = str
Synteny = Iterable[GeneFamily]
OrderedSynteny = Sequence[GeneFamily]
UnorderedSynteny = Set[GeneFamily]

SyntenyMapping = Mapping[TreeNode, Synteny]
OrderedSyntenyMapping = Mapping[TreeNode, OrderedSynteny]
UnorderedSyntenyMapping = Mapping[TreeNode, UnorderedSynteny]

# Regex to match digit groups in a string
DIGITS = re.compile(r"([0-9]+)")


def sort_synteny(synteny: Synteny) -> OrderedSynteny:
    """
    Sort objects of a synteny in a canonical order.

    :param synteny: original synteny, with or without a defined order
    :returns: synteny in canonical order
    """

    def key(obj: str):
        parts = DIGITS.split(obj)
        return [int(part) if part.isdigit() else part for part in parts]

    return sorted(synteny, key=key)


def format_synteny(synteny: Synteny) -> str:
    """
    Convert a synteny to a human-readable string.

    :param synteny: input synteny
    :returns: string representation
    """
    sequence = sort_synteny(synteny) if isinstance(synteny, set) else synteny

    if all(len(family) == 1 for family in synteny):
        return "".join(sequence)
    else:
        return ",".join(sequence)


def parse_synteny_mapping(
    tree: Tree, data: Dict[str, Synteny]
) -> SyntenyMapping:
    """
    Convert a plain mapping of syntenies to tree node names to a mapping
    to node objects.

    :param tree: labeled tree
    :param data: plain mapping to map from
    :returns: parsed mapping
    """
    return {tree & node: synteny for node, synteny in data.items()}


def serialize_synteny_mapping(
    mapping: SyntenyMapping,
) -> Dict[str, OrderedSynteny]:
    """
    Convert a mapping of tree nodes to syntenies to a plain mapping.

    :param mapping: mapping to serialize
    :returns: representation that can be used to serialize the mapping
    """
    return {
        node.name: (
            sort_synteny(synteny) if isinstance(synteny, set) else list(synteny)
        )
        for node, synteny in mapping.items()
    }
