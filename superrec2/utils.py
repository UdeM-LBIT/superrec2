from collections import defaultdict
from ete3 import TreeNode
import graphlib
from typing import DefaultDict, Iterable, List, Mapping, Sequence, TypeVar

A = TypeVar('A')
B = TypeVar('B')


def invert_mapping(mapping: Mapping[A, B]) -> DefaultDict[B, List[A]]:
    """
    Compute the reciprocal of a mapping.

    :param mapping: initial mapping to invert
    :returns: for each value of the input mapping, a list of
        keys for which that value is used
    """
    result = defaultdict(list)

    for key, value in mapping.items():
        result[value].append(key)

    return result


def sort_tree_nodes(nodes: Sequence[TreeNode]) -> Iterable[TreeNode]:
    """
    Sort a set of tree nodes so that no node comes before its children.

    :param nodes: set of nodes to sort
    :returns: topologically sorted set of nodes
    """
    toposort = graphlib.TopologicalSorter()
    nodes_set = set(nodes)

    for node in nodes:
        toposort.add(node)

        for child in node.children:
            if child in nodes_set:
                toposort.add(node, child)

    return toposort.static_order()
