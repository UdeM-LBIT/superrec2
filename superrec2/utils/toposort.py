"""Compute a topological ordering or all possible such orderings."""
from collections import deque
from typing import (
    Deque,
    Dict,
    List,
    Mapping,
    Optional,
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

    for succs in graph.values():
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

    for succs in graph.values():
        for succ in succs:
            starts.discard(succ)
            indeg[succ] += 1

    results = _toposort_all_bt(starts, graph, indeg)

    for subresult in results:
        if len(subresult) != len(graph):
            return []
        subresult.reverse()

    return results


def find_cycle(graph: Mapping[Node, Set[Node]]) -> Optional[List[Node]]:
    """
    Find a directed cycle in a graph.

    :param graph: a dict mapping nodes of the graph to their successor nodes
    :returns: a directed cycle, if there is any
    """
    stack: List[Node] = []
    parents: Dict[Node, Node] = {}

    cycle_start = None
    initial = next(iter(graph.keys()))
    stack.append(initial)
    parents[initial] = initial

    while stack and cycle_start is None:
        current = stack.pop()

        for neighbor in graph[current]:
            if neighbor in parents:
                parents[neighbor] = current
                cycle_start = neighbor
                break

            parents[neighbor] = current
            stack.append(neighbor)

    if cycle_start is None:
        return None

    cycle: List[Node] = [cycle_start]
    current = parents[cycle_start]

    while current not in (cycle_start, parents[current]):
        cycle.append(current)
        current = parents[current]

    return cycle
