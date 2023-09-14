from math import inf
from dataclasses import dataclass
from graphlib import CycleError
from typing import Generic, Iterable, TypeVar


T = TypeVar("T")


@dataclass(frozen=True)
class Edge(Generic[T]):
    start: T
    end: T
    weight: float


def _smallest_rotation(sequence: list[T]) -> list[T]:
    """Rotate a list until the smallest element is in front."""
    smallest = min(sequence)
    index = sequence.index(smallest)
    return sequence[index:] + sequence[:index]


def shortest_paths(
    source: T, nodes: Iterable[T], edges: Iterable[Edge[T]]
) -> tuple[dict[T, float], dict[T, T]]:
    """
    Compute shortest paths from a single source.

    :param source: source node
    :param nodes: list of all nodes
    :param edges: list of directed and weighted edges
    :returns: shortest distance from the source to each node
        and predecessor of each node on the path
    :raises CycleError: if the graph contains a negative-weight cycle
    """
    distance = {node: inf for node in nodes}
    predecessor = {}
    distance[source] = 0

    relaxed = True

    for _ in range(len(distance) - 1):
        relaxed = False

        for edge in edges:
            candidate = distance[edge.start] + edge.weight

            if candidate < distance[edge.end]:
                distance[edge.end] = candidate
                predecessor[edge.end] = edge.start
                relaxed = True

        if not relaxed:
            break
    else:
        for edge in edges:
            if distance[edge.start] + edge.weight < distance[edge.end]:
                predecessor[edge.end] = edge.start
                seen = set([edge.end])
                start = edge.start

                while start not in seen:
                    seen.add(start)
                    start = predecessor[start]

                cycle = [start]
                item = predecessor[start]

                while item != start:
                    cycle = [item] + cycle
                    item = predecessor[item]

                raise CycleError(
                    "negative-weight cycle exists",
                    _smallest_rotation(cycle),
                )

    return distance, predecessor
