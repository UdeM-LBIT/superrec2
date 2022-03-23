"""Structures used for holding reconciliation layout information."""
from enum import Enum, auto
from typing import Mapping, NamedTuple, Optional, Union
from ete3 import TreeNode
from ..model.reconciliation import Event
from ..utils.geometry import Position, Rect


class Orientation(Enum):
    """Tree layout orientation."""

    # Root at the top, leaves at the bottom
    VERTICAL = auto()

    # Root on the left, leaves on the right
    HORIZONTAL = auto()


class Branch(NamedTuple):
    """Branch in a gene tree embedded in a species tree."""

    # Type of branch, i.e. the event that triggered the gene branching
    kind: Event

    # Size and position of the branch node
    rect: Rect

    # Anchor point for the parent node edge
    anchor_parent: Position

    # Anchor point for the left child edge
    anchor_left: Position

    # Anchor point for the right child edge
    anchor_right: Position

    # Anchor point for the unique child edge (for transfers)
    anchor_child: Position

    # Label of the branch node or of the extant gene
    name: str = ""

    # Color of the node and its incoming edge (HTML notation)
    color: str = "000000"

    # Children genes
    # * If this is a full loss branch, either `left` or `right` will be None,
    #   which corresponds to the lost gene copy.
    # * If this is a horizontal gene transfer branch, `left` is the conserved
    #   copy and `right` is the transfered copy.
    # * If this is a leaf branch, both children will be None
    left: Optional[TreeNode] = None
    right: Optional[TreeNode] = None


class PseudoGene:  # pylint:disable=too-few-public-methods
    """Objects used as virtual nodes for lost genes."""


GeneAnchor = Union[TreeNode, PseudoGene]


class SubtreeLayout(NamedTuple):
    """
    Layout information about a subtree of the species tree.

    Each subtree, rooted at a node of the species tree, is
    made up of three main parts:

    - Its left and right subtrees (if not a leaf node).
    - A fork that joins its left and right subtrees and contains
      the gene copies that are sent to the child species.
    - A trunk that protrudes above (or to the left of) the center of the fork
      and contains the branching nodes. At the end of the trunk are anchor
      points that link the subtree to its parent subtree.

    Visual representation of those elements, if laid out top-down:

                            trunk
                           ╭─────╮
                              ╱──── anchor
                           ║  │  ║
                           ║ ┌d┐ ║
              ╭ ╔══════════╝ │ │ ╚══════════╗
         fork │ ║ ┌──────────s─│──────────┐ ║
    thickness │ ║ │ ┌──────────s────────┐ │ ║
              ╰ ║ │ │ ╔═════ ^ ═══════╗ │ │ ║
                ║ ┊ ┊ ║      |        ║ ┊ ┊ ║
               (to left    branch    (to right
                subtree)    nodes     subtree)
    """

    # Size and position of this subtree. This rect includes the
    # fork and the trunk and also the rects of the child subtrees
    rect: Rect

    # Size and position of the trunk
    trunk: Rect

    # Number of gene branches in the fork
    fork_thickness: int

    # Positions of the anchor points
    anchors: Mapping[GeneAnchor, Position]

    # Branching nodes in this subtree
    branches: Mapping[GeneAnchor, Branch]


Layout = Mapping[TreeNode, SubtreeLayout]


class DrawParams(NamedTuple):
    """Parameters for laying out and drawing a reconciliation."""

    # Layout orientation (direction in which the tree grows)
    orientation: Orientation = Orientation.VERTICAL

    # Horizontal unit scale (unit-less parameters are multiples of this unit)
    x_unit: str = "1pt"

    # Vertical unit scale (unit-less parameters are multiples of this unit)
    y_unit: str = "1pt"

    # Minimum space between the outline of the species tree
    # and one of the gene branches it contains
    species_branch_padding: float = 4

    # Minimum space between two gene branches in the tree
    gene_branch_spacing: float = 5

    # Space above trunks
    trunk_overhead: float = 10

    # Minimum space between the two subtrees of an internal species
    min_subtree_spacing: float = 12

    # Vertical space between the fork of the parent node
    # and its children subtrees
    level_spacing: float = 4

    # Thickness of the lines around the outer species tree
    species_border_thickness: str = "1pt"

    # Rounding radius of corners of the species outline
    species_border_rounding: str = "4pt"

    # Thickness of the lines that make up the inner gene tree
    branch_thickness: str = "0.5pt"

    # Reserved space around branches of the gene tree when two lines cross
    branch_outer_thickness: str = "4pt"

    # Space between extant gene names and the end of species branches
    species_leaf_spacing: float = 1

    # Distance of the species labels from the species leaves
    species_label_spacing: float = 10

    # Maximum width (in number of chars) before the species labels will be
    # wrapped (leave as None to disable wrapping)
    species_label_width: Optional[int] = 21

    # Size of the filled circles that represent extant genes
    extant_gene_diameter: float = 3

    # Size of the crosses that represent lost genes
    loss_size: float = 3

    # Maximum width (in number of chars) before the event node labels will be
    # wrapped (leave as None to disable wrapping)
    event_label_width: Optional[int] = 18

    # Minimum size of the hollow circles that represent speciation events
    speciation_size: float = 8

    # Minimum size of the squares that represent duplication events
    duplication_size: float = 8

    # Minimum size of the sideways squares that represent transfer events
    transfer_size: float = 8
