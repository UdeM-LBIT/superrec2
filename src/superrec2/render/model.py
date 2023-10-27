"""Structures used for holding history layout information."""
from enum import Enum, auto
from dataclasses import dataclass, field
from typing import Mapping
from sowing.node import Node
from ..model.history import Host, Event
from ..utils.geometry import Position, Rect


@dataclass
class EventLayout:
    """Layout information for a single event node in a history."""

    # Whether this event forks towards children hosts
    forking: bool

    # Children of this event in the same host
    in_children: list[Node[Event, None]]

    # Children of this event in a different host
    out_children: list[Node[Event, None]]

    # Area spanned by this event
    area: Rect = Rect.zero()


@dataclass
class HostLayout:
    """
    Layout information for a host holding one or more event nodes.

    Visual representation of the layout of a host node:

                           events
                            area
                           ╭─────╮
                              ╱──── anchors
                           ║  │  ║
                           ║ ┌d┐ ║
              ╭ ╔══════════╝ │ │ ╚══════════╗
         fork │ ║ ┌──────────s─│──────────┐ ║
              │ ║ │ ┌──────────s────────┐ │ ║
              ╰ ║ │ │ ╔═══════════════╗ │ │ ║
                ║ ┊ ┊ ║               ║ ┊ ┊ ║
               (to left              (to right
                subtree)              subtree)
    """

    # Information about the host
    host: Host

    # Overall area spanned by this host node and its descendants
    area: Rect = Rect.zero()

    # Area spanned by the event nodes inside this host node
    events_area: Rect = Rect.zero()

    # Area spanned by the event nodes that fork towards other hosts
    fork_events_area: Rect = Rect.zero()

    # Children of this host
    children: list[str] = field(default_factory=list)

    # Events inside this host, with individual layout information
    events: Mapping[Node[Event, None], EventLayout] = field(default_factory=dict)

    # Anchor points to which branches can be connected
    anchors: Mapping[Node[Event, None], Position] = field(default_factory=dict)


Layout = Mapping[str, HostLayout]


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
    host_border_radius: str = "3pt"

    # Thickness of the lines around the host outline
    host_border_thickness: str = ".01pt"

    # Rounding radius of branches in the associate tree
    branch_border_radius: str = "1.5pt"

    ## Layout parameters

    # Layout orientation (direction in which the tree grows)
    orientation: Orientation = Orientation.Vertical

    # Minimum space between two event nodes
    events_spacing: float = 4

    # Space to allocate in the host outline around events
    events_host_padding: float = 3

    # Maximum width (in ems) before the species labels will be
    # wrapped (leave as None to disable wrapping)
    host_label_width: int | None = 12

    # Maximum width (in ems) before the event node labels will be
    # wrapped (leave as None to disable wrapping)
    event_label_width: int | None = 6

    ###

    # Minimum space between the outline of the species tree
    # and one of the gene branches it contains
    species_branch_padding: float = 4

    # Minimum space between two gene branches in the tree
    gene_branch_spacing: float = 5

    # Minimum space between the two subtrees of an internal species
    min_subtree_spacing: float = 12

    # Vertical space between each epoch level
    epoch_spacing: float = 6

    # Thickness of the lines that make up the inner gene tree
    branch_thickness: str = "0.5pt"

    # Reserved space around branches of the gene tree when two lines cross
    branch_outer_thickness: str = "4pt"

    # Space between extant gene names and the end of species branches
    host_leaf_spacing: float = 1

    # Distance of the species labels from the species leaves
    host_label_spacing: float = 4

    # Size of the filled circles that represent extant genes
    extant_gene_diameter: float = 3

    # Size of the crosses that represent lost genes
    loss_size: float = 3

    # Minimum size of the hollow circles that represent speciation events
    speciation_size: float = 8

    # Minimum size of the squares that represent duplication events
    duplication_size: float = 8

    # Minimum size of the sideways squares that represent transfer events
    transfer_size: float = 8
