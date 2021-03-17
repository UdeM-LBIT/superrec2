from ete3 import TreeNode
from collections import deque
from typing import (
    Deque,
    Dict,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    TypeVar,
)


Node = TypeVar("Node")


def toposort(graph: Mapping[Node, Set[Node]]) -> Optional[List[Node]]:
    """
    Sort nodes of a graph in topological order.

    See also :meth:`toposort_all` to get a list of all topological orderings.

    :param graph: a dict mapping nodes of the graph to their successor nodes
    :returns: a topological ordering of the graph if such an ordering exists,
        None otherwise (i.e., when the graph has a directed cycle)
    """
    starts: Deque[Node] = deque(graph)
    indeg: Dict[Node, int] = {node: 0 for node in graph}
    result: List[Node] = []

    for node, succs in graph.items():
        for succ in succs:
            if indeg[succ] == 0:
                starts.remove(succ)
            indeg[succ] += 1

    while starts:
        node_from = starts.popleft()
        result.append(node_from)

        for node_to in graph[node_from]:
            indeg[node_to] -= 1

            if indeg[node_to] == 0:
                starts.append(node_to)

    if len(result) == len(graph):
        return result

    return None


def tree_nodes_toposort(nodes: Sequence[TreeNode]) -> Optional[List[TreeNode]]:
    """
    Sort a set of tree nodes so that no node comes before its children.

    :param nodes: set of nodes to sort
    :returns: topologically sorted set of nodes
    """
    subgraph: Dict[TreeNode, Set[TreeNode]] = {node: set() for node in nodes}

    for node in nodes:
        for child in node.children:
            if child in subgraph:
                subgraph[child].add(node)

    return toposort(subgraph)


def _toposort_all_bt(
    starts: Set[Node],
    graph: Mapping[Node, Set[Node]],
    indeg: Dict[Node, int],
) -> List[List[Node]]:
    """
    Enumerate all topological orderings of a sub-graph.

    :param starts: nodes of the graph that have no more predecessors
    :param graph: original graph
    :param indeg: number of remaining predecessors of each node
    :returns: list of valid orderings for the considered sub-graph
    """
    if not starts:
        return [[]]

    results = []

    for node_from in starts:
        next_starts = set(starts)
        next_starts.remove(node_from)

        for node_to in graph[node_from]:
            indeg[node_to] -= 1

            if indeg[node_to] == 0:
                next_starts.add(node_to)

        for subresult in _toposort_all_bt(next_starts, graph, indeg):
            subresult.append(node_from)
            results.append(subresult)

        for node_to in graph[node_from]:
            indeg[node_to] += 1

    return results


def toposort_all(graph: Mapping[Node, Set[Node]]) -> List[List[Node]]:
    """
    Enumerate all topological orderings of the vertices of a graph.

    :param graph: a dict mapping nodes of the graph to their successor nodes
    :returns: list of valid orderings
    """
    starts: Set[Node] = set(graph)
    indeg: Dict[Node, int] = {node: 0 for node in graph}

    for node, succs in graph.items():
        for succ in succs:
            starts.discard(succ)
            indeg[succ] += 1

    results = _toposort_all_bt(starts, graph, indeg)

    for subresult in results:
        if len(subresult) != len(graph):
            return []
        subresult.reverse()

    return results
