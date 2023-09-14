import pytest
from superrec2.model.graph import CycleError, Edge, shortest_paths


def test_shortest_paths():
    nodes = ["s", "t", "x", "y", "z"]
    edges = [
        Edge(start="s", end="t", weight=6),
        Edge(start="s", end="y", weight=7),
        Edge(start="t", end="x", weight=5),
        Edge(start="t", end="y", weight=8),
        Edge(start="t", end="z", weight=-4),
        Edge(start="x", end="t", weight=-2),
        Edge(start="y", end="x", weight=-3),
        Edge(start="y", end="z", weight=9),
        Edge(start="z", end="s", weight=2),
        Edge(start="z", end="x", weight=7),
    ]
    assert shortest_paths("s", nodes, edges) == (
        {
            "s": 0,
            "t": 2,
            "x": 4,
            "y": 7,
            "z": -2,
        },
        {
            "x": "y",
            "y": "s",
            "z": "t",
            "t": "x",
        },
    )

    edges = [
        Edge(start="s", end="t", weight=6),
        Edge(start="s", end="y", weight=7),
        Edge(start="t", end="x", weight=5),
        Edge(start="t", end="y", weight=-4),
        Edge(start="t", end="z", weight=8),
        Edge(start="x", end="t", weight=-2),
        Edge(start="y", end="x", weight=-3),
        Edge(start="y", end="z", weight=9),
        Edge(start="z", end="s", weight=2),
        Edge(start="z", end="x", weight=7),
    ]

    with pytest.raises(CycleError, match="negative-weight cycle exists") as err:
        shortest_paths("s", nodes, edges)

    assert err.value.args[1] == ["t", "y", "x"]
