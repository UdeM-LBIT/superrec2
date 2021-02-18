from typing import NamedTuple

class Position(NamedTuple):
    x: int
    y: int

    def __add__(self, pos: "Position") -> "Position":
        return Position(self.x + pos.x, self.y + pos.y)

    def __sub__(self, pos: "Position") -> "Position":
        return Position(self.x - pos.x, self.y - pos.y)

    def __str__(self) -> str:
        return f"{self.x},{self.y}"


class Size(NamedTuple):
    w: int
    h: int


class Rect(NamedTuple):
    x: int
    y: int
    w: int
    h: int

    def make_from(position: Position, size: Size) -> "Rect":
        return Rect(
            x=position.x,
            y=position.y,
            w=size.w,
            h=size.h,
        )

    def top_left(self) -> Position:
        return Position(self.x, self.y)

    def top_right(self) -> Position:
        return Position(self.x + self.w, self.y)

    def bottom_right(self) -> Position:
        return Position(self.x + self.w, self.y + self.h)

    def bottom_left(self) -> Position:
        return Position(self.x, self.y + self.h)
