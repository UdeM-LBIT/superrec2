"""2D geometry primitives."""
from typing import Self, Sequence
from math import inf
from dataclasses import dataclass


@dataclass(frozen=True)
class Position:
    """Vector or point on the 2D plane."""

    x: float
    y: float

    @classmethod
    def zero(cls) -> Self:
        return cls(0, 0)

    def __add__(self, pos: Self) -> Self:
        """Add two positions together."""
        return Position(self.x + pos.x, self.y + pos.y)

    def __sub__(self, pos: Self) -> Self:
        """Subtract a position from another."""
        return Position(self.x - pos.x, self.y - pos.y)

    def __mul__(self, scalar: float) -> Self:
        """Multiply this position by a scalar."""
        return Position(self.x * scalar, self.y * scalar)

    def __truediv__(self, scalar: float) -> Self:
        """Divide this position by a scalar."""
        return Position(self.x / scalar, self.y / scalar)

    def __neg__(self) -> Self:
        """Flip this position about the X and Y axes."""
        return Position(-self.x, -self.y)

    def __str__(self) -> str:
        """Print position coordinates."""
        return f"{self.x},{self.y}"

    def __format__(self, spec):
        """Format the position as a coordinate tuple, possibly rounding it."""
        if spec:
            digits = int(spec)
            return f"{round(self.x, digits)},{round(self.y, digits)}"

        return str(self)

    def meet_hv(self, pos: Self) -> Self:
        """
        Intersect a horizontal line from this position with a vertical line
        from another position.
        """
        return Position(pos.x, self.y)

    def meet_vh(self, pos: Self) -> Self:
        """
        Intersect a vertical line from this position with a horizontal line
        from another position.
        """
        return Position(self.x, pos.y)


@dataclass(frozen=True)
class Size:
    """Dimension on the 2D plane."""

    w: float
    h: float

    @classmethod
    def zero(cls) -> Self:
        return cls(0, 0)


@dataclass(frozen=True)
class Rect:
    """Axis-aligned 2D rectangle."""

    position: Position
    size: Size

    @classmethod
    def zero(cls) -> Self:
        return cls(position=Position.zero(), size=Size.zero())

    @classmethod
    def fit(cls, items: Sequence[Position | Self]) -> Self:
        """Compute the smallest rectangle fitting a set of positions or rectangles."""
        min_x = min_y = inf
        max_x = max_y = -inf

        for item in items:
            if isinstance(item, Position):
                positions = (item,)
            else:
                positions = (item.top_left(), item.bottom_right())

            for position in positions:
                min_x = min(min_x, position.x)
                max_x = max(max_x, position.x)
                min_y = min(min_y, position.y)
                max_y = max(max_y, position.y)

        if min_x == min_y == inf and max_x == max_y == -inf:
            raise ValueError("fit() arg is an empty sequence")

        return cls(
            position=Position(min_x, min_y),
            size=Size(max_x - min_x, max_y - min_y),
        )

    def __add__(self, shift: Position) -> Self:
        """Shift the rectangle by adding the given vector."""
        return Rect(position=self.position + shift, size=self.size)

    def __sub__(self, shift: Position) -> Self:
        """Shift the rectangle by subtracting the given vector."""
        return Rect(position=self.position - shift, size=self.size)

    def grow(self, w: float, h: float | None = None) -> Self:
        """
        Add a given amount of padding around the rectangle.

        :param w: horizontal padding added on both the left and right sides
        :param h: vertical padding added on both the top and bottom sides;
            if omitted, use the same amount as horizontal padding
        """
        h = h if h is not None else w
        return Rect(
            position=self.position - Position(w, h),
            size=Size(self.size.w + 2 * w, self.size.h + 2 * h),
        )

    def top_left(self) -> Position:
        """Position of the upper left corner."""
        return self.position

    def top(self) -> Position:
        """Position of the upper edge’s center."""
        return Position(self.position.x + self.size.w / 2, self.position.y)

    def top_right(self) -> Position:
        """Position of the upper left corner."""
        return Position(self.position.x + self.size.w, self.position.y)

    def right(self) -> Position:
        """Position of the right edge’s center."""
        return Position(
            self.position.x + self.size.w, self.position.y + self.size.h / 2
        )

    def bottom_right(self) -> Position:
        """Position of the lower right corner."""
        return Position(self.position.x + self.size.w, self.position.y + self.size.h)

    def bottom(self) -> Position:
        """Position of the lower edge’s center."""
        return Position(
            self.position.x + self.size.w / 2, self.position.y + self.size.h
        )

    def bottom_left(self) -> Position:
        """Position of the lower left corner."""
        return Position(self.position.x, self.position.y + self.size.h)

    def left(self) -> Position:
        """Position of the left edge’s center."""
        return Position(self.position.x, self.position.y + self.size.h / 2)

    def center(self) -> Position:
        """Position of the center."""
        return Position(
            self.position.x + self.size.w / 2,
            self.position.y + self.size.h / 2,
        )
