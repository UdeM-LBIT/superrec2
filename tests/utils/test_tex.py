from superrec2.utils.tex import measure_tikz
from superrec2.utils.geometry import Rect, Position, Size


def test_measure_tikz():
    assert measure_tikz(
        [
            r"\node {abcdef};",
            r"\node[label={right:uvwxyz}] {abcdef};",
        ],
        preamble=r"\usepackage{tikz}",
    ) == [
        Rect(Position(x=-17.50298, y=-6.91296), Size(w=35.00596, h=13.82592)),
        Rect(Position(x=-17.50298, y=-6.91296), Size(w=75.1319, h=13.82592)),
    ]
