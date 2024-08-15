"""Structures used for holding history layout information."""

from enum import Enum, auto
from dataclasses import dataclass, field
from typing import Mapping, Self
from sowing.node import Node
from sowing.indexed import IndexedTree
from itertools import chain
from ..model.history import Associate, Host, Event
from ..utils.geometry import Position, Rect


class Waypoint(Event):
    """Null event inserted into histories for layout purposes."""

    def outcomes(self, host_index: IndexedTree[Host, None]) -> tuple[Associate, ...]:
        return ()

    def validate(
        self,
        host_index: IndexedTree[Host, None],
        children: tuple[Associate, ...],
    ) -> None:
        Event.validate(self, host_index, children)


class Orientation(Enum):
    """Tree layout orientation."""

    # Root at the top, leaves at the bottom
    Vertical = auto()

    # Root on the left, leaves on the right
    Horizontal = auto()


@dataclass(frozen=True)
class DrawParams:
    """Parameters for laying out and drawing a reconciliation."""

    ## Drawing parameters

    # Horizontal unit scale (unit-less parameters are multiples of this unit)
    x_unit: str = "1pt"

    # Vertical unit scale (unit-less parameters are multiples of this unit)
    y_unit: str = "1pt"

    # Rounding radius of corners of the host outline
    host_border_radius: str = "2pt"

    # Thickness of the lines around the host outline
    host_border_thickness: str = ".01pt"

    # Thickness of the lines of the associate tree
    branch_thickness: str = "0.5pt"

    # Rounding radius of branches in the associate tree
    branch_border_radius: str = "1.5pt"

    # Draw rectangles around computed areas to aid in debugging
    debug: bool = False

    ## Layout parameters

    # Layout orientation (direction of the main axis along which the tree grows)
    orientation: Orientation = Orientation.Vertical

    # Minimum space between two event nodes
    events_spacing: float = 3

    # Space to allocate in the host outline around events
    events_host_padding: float = 3

    # Maximum width (in ems) before the species labels will be
    # wrapped (leave as None to disable wrapping)
    host_label_width: int | None = 12

    # Maximum width (in ems) before the event node labels will be
    # wrapped (leave as None to disable wrapping)
    event_label_width: int | None = 6

    # Space to leave between two sibling hosts
    subtree_spacing: float = 4

    # Vertical space between each epoch level
    epoch_spacing: float = 4

    # Space between extant gene names and the end of host outlines
    host_leaf_spacing: float = 1

    # Size of the filled circles that represent extant associates
    extant_diameter: float = 3

    # Minimum size of the hollow circles that represent speciation events
    speciation_size: float = 8

    # Minimum size of the squares that represent duplication events
    duplication_size: float = 8

    # Minimum size of the diamonds that represent transfer events
    transfer_size: float = 8


@dataclass
class EventLayout:
    """Layout information for a single event node in a history."""

    # Children of this event in the same host
    in_children: list[Node[Event, None]]

    # Children of this event in descending hosts
    desc_children: list[Node[Event, None]]

    # Children of this event in parallel hosts
    side_children: list[Node[Event, None]]

    # True if the transfer edges to side children can be drawn horizontally
    side_horizontal: bool = False

    # Area spanned by this event’s node
    area: Rect = Rect.zero()

    # Position of the event anchor point
    anchor: Position = Position.zero()

    @property
    def children(self) -> list[Node[Event, None]]:
        return self.in_children + self.desc_children + self.side_children

    @property
    def forking(self) -> bool:
        return bool(self.desc_children)

    @property
    def leaf(self) -> bool:
        return not self.children or not (self.in_children or self.side_horizontal)

    def __iadd__(self, shift: Position) -> Self:
        """Shift the event position by adding the given vector."""
        self.area += shift
        self.anchor += shift
        return self

    def __isub__(self, shift: Position) -> Self:
        """Shift the event position by subtracting the given vector."""
        self += -shift
        return self


@dataclass
class HostLayout:
    """
    Layout information for a host holding one or more event nodes.

    Visual representation of the layout of a host node:

                           events
                            area
                           ╭─────╮
                              <──── anchors
              ╭            ║  │  ║
        trunk │            ║  │  ║
              ╰            ║ ┌d┐ ║
              ╭ ╔══════════╝ │ │ ╚══════════╗
         fork │ ║ ┌──────────s─│──────────┐ ║
              │ ║ │ ┌──────────s────────┐ │ ║
              ╰ ║ │ │ ╔═══════════════╗ │ │ ║
                ║ ┊ ┊ ║               ║ ┊ ┊ ║
               (to left              (to right
                subtree)              subtree)
    """

    # Drawing parameters
    params: DrawParams

    # Information about the host
    host: Host

    # Children of this host with their layout information
    children: Mapping[str, Self] = field(default_factory=dict)

    # Events inside this host, with individual layout information
    events: Mapping[Node[Event, None], EventLayout] = field(default_factory=dict)

    # Anchor points to which branches can be connected
    anchors: Mapping[Node[Event, None], Position] = field(default_factory=dict)

    @property
    def area(self) -> Rect:
        """Overall area of the host, including its children."""
        return Rect.fit(
            chain(
                (self.trunk_area,),
                (layout.area for layout in self.children.values()),
            )
        )

    @property
    def fork_area(self) -> Rect:
        """Area spanned by the forking and leaf events in this host."""
        return Rect.fit(
            (
                layout.area
                for layout in self.events.values()
                if layout.forking or layout.leaf
            )
        ).grow(self.params.events_host_padding)

    @property
    def trunk_area(self) -> Rect:
        """Area spanned by all events in this host."""
        inside = [
            self.fork_area,
            self.fork_area.top() + Position(0, -self.params.epoch_spacing),
        ]

        if any(
            not layout.forking and not layout.leaf for layout in self.events.values()
        ):
            inside.append(
                Rect.fit(
                    (
                        layout.area
                        for layout in self.events.values()
                        if not layout.forking and not layout.leaf
                    )
                ).grow(w=self.params.events_host_padding, h=0)
            )

        return Rect.fit(inside)

    def __iadd__(self, shift: Position) -> Self:
        """Shift the host position by adding the given vector."""
        for child in self.children:
            self.children[child] += shift

        for event in self.events:
            self.events[event] += shift

        for anchor in self.anchors:
            self.anchors[event] += shift

        return self

    def __isub__(self, shift: Position) -> Self:
        """Shift the host position by subtracting the given vector."""
        self += -shift
        return self


Layout = Mapping[str, HostLayout]
