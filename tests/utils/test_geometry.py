import pytest
from superrec2.utils.geometry import Position, Size, Rect


def test_position():
    pos1 = Position(8, 6)
    pos2 = Position(2, -8)

    assert pos1.x == 8
    assert pos1.y == 6
    assert Position.zero().x == 0
    assert Position.zero().y == 0

    assert pos1 + Position.zero() == pos1
    assert pos1 + pos2 == Position(10, -2)
    assert pos1 - pos2 == Position(6, 14)
    assert -pos1 == Position(-8, -6)
    assert pos1 * 2 == Position(16, 12)
    assert pos1 / 2 == Position(4, 3)

    assert str(pos1) == "8,6"
    assert str(pos2) == "2,-8"
    assert str(-pos2) == "-2,8"
    assert str(-pos1) == "-8,-6"

    pos3 = Position(8.5555, 6.5555)
    assert str(pos3) == "8.5555,6.5555"
    assert f"{pos3:1}" == "8.6,6.6"

    assert pos1.meet_hv(pos2) == Position(2, 6)
    assert pos1.meet_vh(pos2) == Position(8, -8)
    assert pos2.meet_hv(pos1) == pos1.meet_vh(pos2)
    assert pos2.meet_vh(pos1) == pos1.meet_hv(pos2)

    assert sum((pos1, pos2, pos3), start=Position.zero()) == pos1 + pos2 + pos3


def test_rect():
    rect1 = Rect(Position(8, 6), Size(10, 20))
    rect2 = Rect(Position(2, -8), Size(5, 3))

    assert rect1 + Position(1, 1) == Rect(Position(9, 7), Size(10, 20))
    assert rect1 - Position(1, 1) == Rect(Position(7, 5), Size(10, 20))
    assert rect1 + Position.zero() == rect1
    assert rect1 - Position.zero() == rect1

    assert Rect.zero() == Rect(Position(0, 0), Size(0, 0))

    with pytest.raises(ValueError, match=r"fit\(\) arg is an empty sequence"):
        Rect.fit(())

    assert Rect.fit((Position(8, 6), Position(2, -8), Position(3, 12))) == Rect(
        position=Position(2, -8),
        size=Size(6, 20),
    )

    assert Rect.fit((rect1, rect2)) == Rect(
        position=Position(2, -8),
        size=Size(16, 34),
    )

    assert rect1.top_left() == Position(8, 6)
    assert rect1.top() == Position(13, 6)
    assert rect1.top_right() == Position(18, 6)
    assert rect1.right() == Position(18, 16)
    assert rect1.bottom_right() == Position(18, 26)
    assert rect1.bottom() == Position(13, 26)
    assert rect1.bottom_left() == Position(8, 26)
    assert rect1.left() == Position(8, 16)
    assert rect1.center() == Position(13, 16)
    assert rect2.center() == Position(4.5, -6.5)
