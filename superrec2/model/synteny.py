import re
from typing import Dict, Iterable, List, Mapping, Sequence, Set
from textwrap import wrap
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


def _wrap_badness(lines: List[str]) -> int:
    max_len = max(len(line) for line in lines)
    return sum((len(line) - max_len) ** 2 for line in lines)


def format_synteny(synteny: Synteny, width: int = 21) -> str:
    """
    Convert a synteny to a human-readable string.

    :param synteny: input synteny
    :param width: produce lines not exceeding the given character width
    :returns: string representation
    """
    source = ", ".join(
        sort_synteny(synteny) if isinstance(synteny, set) else synteny
    )

    best_result = wrap(source, width)
    best_badness = _wrap_badness(best_result)
    line_count = len(best_result)

    # Reduce wrap width looking for the most balanced solution
    # while keeping the same number of lines
    while width > 1:
        width -= 1
        next_result = wrap(source, width)

        if len(next_result) != line_count:
            break

        next_badness = _wrap_badness(next_result)

        if next_badness < best_badness:
            best_result = next_result
            best_badness = next_badness

    return "\n".join(best_result)


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
