"""2D geometry primitives."""
from typing import NamedTuple


class Position(NamedTuple):
    """Vector or point on the 2D plane."""

    x: float
    y: float

    def __add__(self, pos: tuple) -> "Position":
        """Add two vectors together."""
        return Position(self.x + pos[0], self.y + pos[1])

    def __sub__(self, pos: tuple) -> "Position":
        """Subtract a vector from another."""
        return Position(self.x - pos[0], self.y - pos[1])

    def __str__(self) -> str:
        """Print vector coordinates."""
        return f"{self.x},{self.y}"

    def __format__(self, spec):
        """Format the vector as a coordinate tuple, possibly rounding it."""
        if spec:
            digits = int(spec)
            return f"{round(self.x, digits)},{round(self.y, digits)}"

        return str(self)

    def meet_hv(self, pos: tuple) -> "Position":
        """
        Intersect a horizontal line from this position with a vertical line
        from another position.
        """
        return Position(pos[0], self.y)

    def meet_vh(self, pos: tuple) -> "Position":
        """
        Intersect a vertical line from this position with a horizontal line
        from another position.
        """
        return Position(self.x, pos[1])


class Size(NamedTuple):
    """Dimension on the 2D plane."""

    w: float
    h: float


class Rect(NamedTuple):
    """Axis-aligned 2D rectangle."""

    x: float
    y: float
    w: float
    h: float

    @staticmethod
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

    def __add__(self, pos: tuple) -> "Rect":
        """Shift the rectangle by adding the given vector."""
        return Rect(x=self.x + pos[0], y=self.y + pos[1], w=self.w, h=self.h)

    def __sub__(self, pos: tuple) -> "Rect":
        """Shift the rectangle by subtracting the given vector."""
        return Rect(x=self.x - pos[0], y=self.y - pos[1], w=self.w, h=self.h)

    def top_left(self) -> Position:
        """Position of the upper left corner."""
        return Position(self.x, self.y)

    def top(self) -> Position:
        """Position of the upper edge’s center."""
        return Position(self.x + self.w / 2, self.y)

    def top_right(self) -> Position:
        """Position of the upper left corner."""
        return Position(self.x + self.w, self.y)

    def right(self) -> Position:
        """Position of the right edge’s center."""
        return Position(self.x + self.w, self.y + self.h / 2)

    def bottom_right(self) -> Position:
        """Position of the lower right corner."""
        return Position(self.x + self.w, self.y + self.h)

    def bottom(self) -> Position:
        """Position of the lower edge’s center."""
        return Position(self.x + self.w / 2, self.y + self.h)

    def bottom_left(self) -> Position:
        """Position of the lower left corner."""
        return Position(self.x, self.y + self.h)

    def left(self) -> Position:
        """Position of the left edge’s center."""
        return Position(self.x, self.y + self.h / 2)

    def center(self) -> Position:
        """Position of the center."""
        return Position(self.x + self.w / 2, self.y + self.h / 2)
