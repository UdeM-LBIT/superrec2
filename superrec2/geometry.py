from typing import NamedTuple

class Position(NamedTuple):
    """Vector or point on the 2D plane."""
    x: int
    y: int

    def __add__(self, pos: "Position") -> "Position":
        """Add two vectors together."""
        return Position(self.x + pos.x, self.y + pos.y)

    def __sub__(self, pos: "Position") -> "Position":
        """Subtract a vector from another."""
        return Position(self.x - pos.x, self.y - pos.y)

    def __str__(self) -> str:
        """Print vector coordinates."""
        return f"{self.x},{self.y}"


class Size(NamedTuple):
    """Dimension on the 2D plane."""
    w: int
    h: int


class Rect(NamedTuple):
    """Axis-aligned 2D rectangle."""
    x: int
    y: int
    w: int
    h: int

    def make_from(position: Position, size: Size) -> "Rect":
        """
        Create a rectangle from a position and a dimension.

        :param position: position of the rectangle’s top left corner
        :param size: overall size of the rectangle
        """
        return Rect(
            x=position.x,
            y=position.y,
            w=size.w,
            h=size.h,
        )

    def top_left(self) -> Position:
        """Position of the upper left corner."""
        return Position(self.x, self.y)

    def top_right(self) -> Position:
        """Position of the upper left corner."""
        return Position(self.x + self.w, self.y)

    def bottom_right(self) -> Position:
        """Position of the lower right corner."""
        return Position(self.x + self.w, self.y + self.h)

    def bottom_left(self) -> Position:
        """Position of the lower left corner."""
        return Position(self.x, self.y + self.h)
