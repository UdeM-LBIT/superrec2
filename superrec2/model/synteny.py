import re
from typing import Dict, List, Mapping, Sequence
from ete3 import Tree, TreeNode
from ..utils.escape import escape, unescape


Synteny = Sequence[str]
SyntenyMapping = Mapping[TreeNode, Synteny]

# Set of characters that must be escaped from names
ESCAPE_CHARS = r"\,()"


def parse_synteny(data: str) -> Synteny:
    """
    Parse a string representation of a synteny.

    :param data: string to parse
    :returns: parsed synteny
    :example:
    >>> parse_synteny("abcdef")
    ["a", "b", "c", "d", "e", "f"]
    >>> parse_synteny("(gene1)(gene2)xyz(gene3)")
    ["gene1", "gene2", "x", "y", "z", "gene3"]
    """
    depth = 0
    escaped = False
    item = ""
    result = []

    for char in data:
        if not escaped and char == "\\":
            escaped = True
            continue

        if depth == 0:
            if escaped:
                result.append(char)
            elif char == "(":
                depth += 1
            elif not str.isspace(char):
                result.append(char)
        else:
            if escaped:
                item += char
            elif char == "(":
                item += char
                depth += 1
            elif char == ")":
                depth -= 1

                if depth == 0:
                    result.append(item.strip())
                    item = ""
                else:
                    item += char
            else:
                item += char

        escaped = False

    return result


def serialize_synteny(synteny: Synteny) -> str:
    """
    Serialize a synteny.

    :param mapping: synteny to serialize
    :returns: serialized representation
    :example:
    >>> serialize_synteny(["a", "b", "c", "d", "e", "f"])
    "abcdef"
    >>> serialize_synteny(["gene1", "gene2", "x", "y", "z", "gene3"])
    "(gene1)(gene2)xyz(gene3)"
    """
    return "".join(
        f"({escape(item, ESCAPE_CHARS)})"
        if len(item) >= 2
        else escape(item, ESCAPE_CHARS)
        for item in synteny
    )


MAPPING_SPLIT_REGEX = re.compile(r"(?<!\\),")


def parse_synteny_mapping(tree: Tree, data: str) -> SyntenyMapping:
    """
    Parse a string representation of a synteny mapping.

    :param tree: labeled tree
    :param data: string to parse
    :returns: parsed mapping
    :example:
    >>> parse_synteny_mapping(Tree("(x, y);"), "x: abcd, y: ade")
    {Tree node "x": ["a", "b", "c", "d"],
     Tree node "y": ["a", "d", "e"]}
    """
    result: Dict[TreeNode, List[str]] = {}

    for pair in MAPPING_SPLIT_REGEX.split(data):
        if pair.strip():
            node, synteny = map(lambda x: x.strip(), pair.split(":"))

            if not node.startswith("#"):
                result[tree & unescape(node)] = parse_synteny(synteny)

    return result


def serialize_synteny_mapping(mapping: SyntenyMapping) -> str:
    """
    Serialize a synteny mapping.

    :param mapping: mapping to serialize
    :returns: serialized representation
    >>> tree = Tree("(x, y);")
    >>> serialize_synteny_mapping({
          tree & "x": ["a", "b", "c", "d"],
          tree & "y": ["a", "d", "e"]
        })
    "x:abcd,y:ade"
    """
    return ",".join(
        f"{escape(name, ESCAPE_CHARS)}:{serialize_synteny(synteny)}"
        for name, synteny in sorted(
            (node.name, synteny) for node, synteny in mapping.items()
        )
    )
