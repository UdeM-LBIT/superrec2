"""Automatically draw reconciliations with optional labels."""
from enum import Enum, auto
import textwrap
from typing import Any, Dict, List, Mapping, NamedTuple, Optional, Union
from ete3 import Tree, TreeNode
from .utils.geometry import Position, Rect, Size
from .utils.trees import LowestCommonAncestor
from .utils import tex
from .model.reconciliation import (
    NodeEvent,
    EdgeEvent,
    ReconciliationOutput,
    SuperReconciliationOutput,
)
from .model.tree_mapping import TreeMapping
from .model.synteny import Synteny, SyntenyMapping


class Orientation(Enum):
    """Tree layout orientation."""

    VERTICAL = auto()
    HORIZONTAL = auto()


class Branch(NamedTuple):
    """Branch in a gene tree embedded in a species tree."""

    # Type of branch, i.e. the event that triggered the gene branching
    kind: Union[NodeEvent, EdgeEvent]

    # Size and position of the branch node
    rect: Rect

    # Label of the branch node or of the extant gene
    name: str = ""

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
    """Parameters for drawing a reconciliation."""

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
    min_subtree_spacing: float = 25

    # Vertical space between the fork of the parent node
    # and its children subtrees
    level_spacing: float = 0

    # Thickness of the lines around the outer species tree
    species_border_thickness: str = "1pt"

    # Thickness of the lines that make up the inner gene tree
    branch_thickness: str = "0.5pt"

    # Reserved space around branches of the gene tree when two lines cross
    branch_outer_thickness: str = "4pt"

    # Space between extant gene names and the end of species branches
    species_leaf_spacing: float = 1

    # Distance of the species labels from the species leaves
    species_label_spacing: float = 10

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

    def get_tikz_definitions(self):
        """Get TikZ definitions matching a set of drawing parameters."""
        if self.orientation == Orientation.VERTICAL:
            leaf_label_style = textwrap.dedent(
                """
                [font={{\strut}},
                    fill=white,
                    inner xsep=0pt, inner ysep=2pt,
                    outer xsep=0pt, outer ysep=0pt]
                below:#1
                """
            )
        else:
            leaf_label_style = textwrap.dedent(
                """
                [fill=white,
                    inner xsep=4pt, inner ysep=0pt,
                    outer xsep=0pt, outer ysep=0pt]
                right:#1
                """
            )

        return textwrap.dedent(
            rf"""
            \tikzset{{
                x={{{self.x_unit}}},
                y={{-{self.y_unit}}},
                species border/.style={{
                    line width={{{self.species_border_thickness}}},
                    shorten <={{-{self.species_border_thickness} / 2 + 0.05pt}},
                    shorten >={{-{self.species_border_thickness} / 2 + 0.05pt}},
                }},
                species label/.style={{
                    font=\bfseries,
                    midway,
                    {
                        "yshift=-" if self.orientation == Orientation.VERTICAL
                        else "xshift="
                    }{self.species_label_spacing}
                }},
                branch/.style={{
                    line width={{{self.branch_thickness}}},
                    preaction={{
                        draw, white, -{{}},
                        line width={{{self.branch_outer_thickness}}},
                        shorten <={{{self.branch_thickness}}},
                        shorten >={{{self.branch_thickness}}},
                    }},
                }},
                transfer branch/.style={{
                    branch,
                    -Stealth,
                }},
                loss/.style={{
                    draw, cross out, thick,
                    line width={{{self.branch_thickness}}},
                    inner sep=0pt,
                    outer sep=0pt,
                    minimum width={{{self.loss_size}}},
                    minimum height={{{self.loss_size}}},
                }},
                extant gene/.style={{
                    circle, fill,
                    outer sep=0pt, inner sep=0pt,
                    minimum size={{{self.extant_gene_diameter}}},
                    label={{
                        {textwrap.indent(leaf_label_style, " " * 24).strip()}
                    }},
                }},
                branch node/.style={{
                    draw, fill=white,
                    outer sep=0pt, inner xsep=0pt, inner ysep=2pt,
                    line width={{{self.branch_thickness}}},
                }},
                speciation/.style={{
                    branch node, rounded rectangle,
                    inner xsep=4pt,
                    minimum width={{{self.speciation_size}}},
                    minimum height={{{self.speciation_size}}},
                }},
                duplication/.style={{
                    branch node, rectangle,
                    inner xsep=4pt,
                    minimum width={{{self.duplication_size}}},
                    minimum height={{{self.duplication_size}}},
                }},
                horizontal gene transfer/.style={{
                    branch node, signal, signal to=east and west,
                    minimum width={{{self.transfer_size}}},
                    minimum height={{{self.transfer_size}}},
                }},
            }}"""
        ).lstrip()


def _add_losses(
    layout_state: Dict[TreeNode, Dict],
    gene: TreeNode,
    start_species: TreeNode,
    end_species: TreeNode,
) -> GeneAnchor:
    """
    Insert virtual gene loss nodes between
    a parent species and a child species.

    :param gene: parent gene that is lost
    :param start_species: lower species in which the gene is conserved
    :param end_species: parent of the species from which the
        gene originated
    :returns: first virtual child created in the process
    """
    prev_gene = gene
    prev_species = start_species
    start_species = start_species.up

    while start_species != end_species:
        is_left = prev_species == start_species.children[0]
        is_right = prev_species == start_species.children[1]
        state = layout_state[start_species]
        cur_gene = PseudoGene()

        state["anchor_nodes"].add(cur_gene)
        state["branches"][cur_gene] = {
            "kind": EdgeEvent.FULL_LOSS,
            "name": "",
            "left": prev_gene if is_left else None,
            "right": prev_gene if is_right else None,
        }

        prev_gene = cur_gene
        prev_species = start_species
        start_species = start_species.up

    return prev_gene


def _format_synteny(synteny: Synteny) -> str:
    if all(len(family) == 1 for family in synteny):
        return "".join(map(tex.escape, synteny))
    else:
        return ",".join(map(tex.escape, synteny))


def _compute_branches(  # pylint:disable=too-many-locals
    layout_state: Dict[TreeNode, Dict],
    rec: ReconciliationOutput,
) -> None:
    """Create the branching nodes for each species."""
    gene_tree = rec.input.object_tree
    species_lca = rec.input.species_lca
    species_tree = species_lca.tree
    mapping = rec.object_species
    syntenies = (
        rec.syntenies if isinstance(rec, SuperReconciliationOutput) else {}
    )

    for root_species in species_tree.traverse("postorder"):
        state: Dict[str, Any] = {
            "branches": {},
            "anchor_nodes": set(),
        }
        layout_state[root_species] = state

        for root_gene in gene_tree.traverse("postorder"):
            if mapping[root_gene] != root_species:
                continue

            if root_gene.is_leaf():
                # Create branches even for leaf genes
                if root_gene in syntenies:
                    name = _format_synteny(syntenies[root_gene])
                else:
                    species_name, gene_name = root_gene.name.split("_")
                    name = (
                        rf"{tex.escape(species_name)}"
                        rf"\textsubscript{{{tex.escape(gene_name)}}}"
                    )

                state["anchor_nodes"].add(root_gene)
                state["branches"][root_gene] = {
                    "kind": NodeEvent.LEAF,
                    "name": name,
                }
            else:
                # Create branches for actual internal nodes
                left_gene, right_gene = root_gene.children
                event = rec.node_event(root_gene)

                name = (
                    rf"{_format_synteny(syntenies[root_gene])}"
                    if root_gene in syntenies
                    else ""
                )

                if event == NodeEvent.SPECIATION:
                    # Speciation nodes are located below the trunk
                    # and linked to child species’s gene anchors
                    left_species = root_species.children[0]

                    if species_lca.is_ancestor_of(
                        left_species, mapping[right_gene]
                    ):
                        # Left gene and right gene are swapped relative
                        # to the left and right species
                        left_gene, right_gene = right_gene, left_gene

                    left_gene = _add_losses(
                        layout_state,
                        left_gene,
                        mapping[left_gene],
                        root_species,
                    )
                    right_gene = _add_losses(
                        layout_state,
                        right_gene,
                        mapping[right_gene],
                        root_species,
                    )

                    state["anchor_nodes"].add(root_gene)
                    state["branches"][root_gene] = {
                        "kind": NodeEvent.SPECIATION,
                        "name": name,
                        "left": left_gene,
                        "right": right_gene,
                    }
                elif event == NodeEvent.DUPLICATION:
                    # Duplications are located in the trunk and linked
                    # to other nodes in the same species
                    left_gene = _add_losses(
                        layout_state,
                        left_gene,
                        mapping[left_gene],
                        root_species.up,
                    )
                    right_gene = _add_losses(
                        layout_state,
                        right_gene,
                        mapping[right_gene],
                        root_species.up,
                    )

                    state["anchor_nodes"].add(root_gene)
                    state["anchor_nodes"].remove(left_gene)
                    state["anchor_nodes"].remove(right_gene)
                    state["branches"][root_gene] = {
                        "kind": NodeEvent.DUPLICATION,
                        "name": name,
                        "left": left_gene,
                        "right": right_gene,
                    }
                elif event == NodeEvent.HORIZONTAL_TRANSFER:
                    # Transfers are located in the trunk, like duplications,
                    # but are linked to a node outside the current subtree
                    conserv_gene, foreign_gene = (
                        (left_gene, right_gene)
                        if species_lca.is_ancestor_of(
                            root_species, mapping[left_gene]
                        )
                        else (right_gene, left_gene)
                    )
                    conserv_gene = _add_losses(
                        layout_state,
                        conserv_gene,
                        mapping[conserv_gene],
                        root_species.up,
                    )

                    state["anchor_nodes"].add(root_gene)
                    state["anchor_nodes"].remove(conserv_gene)
                    state["branches"][root_gene] = {
                        "kind": NodeEvent.HORIZONTAL_TRANSFER,
                        "name": name,
                        "left": conserv_gene,
                        "right": foreign_gene,
                    }
                else:
                    raise ValueError("Invalid event")


def _layout_measure_nodes(
    layout_state: Dict[TreeNode, Dict],
    species_tree: Tree,
    params: DrawParams,
) -> Dict[TreeNode, tex.MeasureBox]:
    """Compute the rendered size of each branching node."""
    pending_measures = {}

    for root_species in species_tree.traverse():
        layout = layout_state[root_species]

        for root_gene, branch in layout["branches"].items():
            if branch["kind"] == NodeEvent.LEAF:
                pending_measures[
                    root_gene
                ] = rf"\tikz\node[extant gene={{{branch['name']}}}] {{}};"
            elif branch["kind"] == EdgeEvent.FULL_LOSS:
                continue
            else:
                if branch["kind"] == NodeEvent.SPECIATION:
                    node_type = "[speciation]"
                elif branch["kind"] == NodeEvent.DUPLICATION:
                    node_type = "[duplication]"
                elif branch["kind"] == NodeEvent.HORIZONTAL_TRANSFER:
                    node_type = "[horizontal gene transfer]"
                else:
                    raise ValueError("Invalid node type")

                pending_measures[
                    root_gene
                ] = rf"\tikz\node{node_type} {{{branch['name']}}};"

    return dict(
        zip(
            pending_measures.keys(),
            tex.measure(
                pending_measures.values(),
                preamble=(
                    r"\usepackage{tikz}"
                    r"\usetikzlibrary{arrows.meta}"
                    r"\usetikzlibrary{shapes}" + params.get_tikz_definitions()
                ),
            ),
        )
    )


def _layout_branches(  # pylint:disable=too-many-locals
    layout_state: Dict[TreeNode, Dict],
    species_tree: Tree,
    params: DrawParams,
):
    """Compute the size and relative position of each branch."""
    measures = _layout_measure_nodes(layout_state, species_tree, params)

    for root_species in species_tree.traverse():
        next_pos_across: float = 0
        next_pos_sequence = params.species_branch_padding
        layout = layout_state[root_species]
        layout["anchors"] = {}

        for root_gene, branch in layout["branches"].items():
            size = (
                measures[root_gene].overall_size()
                if root_gene in measures
                else Size(0, 0)
            )

            if branch["kind"] == NodeEvent.LEAF:
                if params.orientation == Orientation.VERTICAL:
                    next_pos_across -= size.w
                    pos = Position(next_pos_across, -size.h)
                else:
                    next_pos_across -= size.h
                    pos = Position(-size.w, next_pos_across)

                next_pos_across -= params.gene_branch_spacing
            elif branch["kind"] in (
                NodeEvent.SPECIATION,
                EdgeEvent.FULL_LOSS,
            ):
                if params.orientation == Orientation.VERTICAL:
                    next_pos_across -= size.w
                    pos = Position(next_pos_across, next_pos_sequence)
                    next_pos_sequence += size.h
                else:
                    next_pos_across -= size.h
                    pos = Position(next_pos_sequence, next_pos_across)
                    next_pos_sequence += size.w

                next_pos_across -= params.gene_branch_spacing
                next_pos_sequence += params.gene_branch_spacing
            elif branch["kind"] == NodeEvent.DUPLICATION:
                left_rect = layout["branches"][branch["left"]]["rect"]
                right_rect = layout["branches"][branch["right"]]["rect"]

                if params.orientation == Orientation.VERTICAL:
                    across = (
                        (left_rect.center() + right_rect.center()).x
                        - size.w
                    ) / 2
                    sequence = min(
                        params.species_branch_padding,
                        left_rect.y,
                        right_rect.y
                    ) - params.species_branch_padding - size.h
                    pos = Position(across, sequence)
                else:
                    across = (
                        (left_rect.center() + right_rect.center()).y
                        - size.h
                    ) / 2
                    sequence = min(
                        params.species_branch_padding,
                        left_rect.x,
                        right_rect.x
                    ) - params.species_branch_padding - size.w
                    pos = Position(sequence, across)
            elif branch["kind"] == NodeEvent.HORIZONTAL_TRANSFER:
                cons_rect = layout["branches"][branch["left"]]["rect"]

                if params.orientation == Orientation.VERTICAL:
                    across = cons_rect.center().x - size.w / 2
                    sequence = (
                        min(params.species_branch_padding, cons_rect.y)
                        - params.species_branch_padding
                        - size.h
                    )
                    pos = Position(across, sequence)
                else:
                    across = cons_rect.center().y - size.h / 2
                    sequence = (
                        min(params.species_branch_padding, cons_rect.x)
                        - params.species_branch_padding
                        - size.w
                    )
                    pos = Position(sequence, across)
            else:
                raise ValueError("Invalid node type")

            rect = Rect.make_from(pos, size)
            branch["rect"] = rect

            if root_gene in layout["anchor_nodes"]:
                if params.orientation == Orientation.VERTICAL:
                    layout["anchors"][root_gene] = Position(rect.center().x, 0)
                else:
                    layout["anchors"][root_gene] = Position(0, rect.center().y)

        del layout["anchor_nodes"]

        # Shift all nodes to the left or up to make room
        # for the initial padding
        if layout["branches"]:
            if params.orientation == Orientation.VERTICAL:
                padding_shift = Position(
                    x=(
                        min(
                            -branch["rect"].right().x
                            for branch in layout["branches"].values()
                        )
                        - params.species_branch_padding
                    ),
                    y=0,
                )
            else:
                padding_shift = Position(
                    x=0,
                    y=(
                        min(
                            -branch["rect"].bottom().y
                            for branch in layout["branches"].values()
                        )
                        - params.species_branch_padding
                    ),
                )

            for root_gene in layout["branches"]:
                branch = layout["branches"][root_gene]
                branch["rect"] += padding_shift

            for root_gene in layout["anchors"]:
                layout["anchors"][root_gene] += padding_shift


def _layout_subtrees(
    layout_state: Dict[TreeNode, Dict],
    species_tree: Tree,
    params: DrawParams,
):
    """Compute the size and absolute position of each subtree."""
    # Compute the size of each subtree
    for root_species in species_tree.traverse("postorder"):
        state = layout_state[root_species]

        if state["branches"]:
            if params.orientation == Orientation.VERTICAL:
                trunk_width = max(
                    -branch["rect"].top_left().x
                    for branch in state["branches"].values()
                ) + params.species_branch_padding
                trunk_height = max(0, max(
                    -branch["rect"].top_left().y
                    for branch in state["branches"].values()
                )) + params.trunk_overhead
                fork_thickness = max(0, max(
                    branch["rect"].bottom_right().y
                    for branch in state["branches"].values()
                )) + params.species_branch_padding
            else:
                trunk_width = max(0, max(
                    -branch["rect"].top_left().x
                    for branch in state["branches"].values()
                )) + params.trunk_overhead
                trunk_height = max(
                    -branch["rect"].top_left().y
                    for branch in state["branches"].values()
                ) + params.species_branch_padding
                fork_thickness = max(0, max(
                    branch["rect"].bottom_right().x
                    for branch in state["branches"].values()
                )) + params.species_branch_padding
        else:
            # Empty subtree
            fork_thickness = 0

            if params.orientation == Orientation.VERTICAL:
                trunk_width = 0
                trunk_height = params.trunk_overhead
            else:
                trunk_width = params.trunk_overhead
                trunk_height = 0

        trunk_size = Size(trunk_width, trunk_height)

        if root_species.is_leaf():
            # Extant species
            state["size"] = trunk_size
            state["trunk"] = Rect.make_from(Position(0, 0), trunk_size)
            state["fork_thickness"] = 0
        else:
            # Ancestral species
            left_species, right_species = root_species.children
            left_info = layout_state[left_species]
            right_info = layout_state[right_species]

            if params.orientation == Orientation.VERTICAL:
                subtree_span = (
                    max(left_info["size"].h, right_info["size"].h)
                    + trunk_height
                )
            else:
                subtree_span = (
                    max(left_info["size"].w, right_info["size"].w)
                    + trunk_width
                )

            subtree_span += params.level_spacing + fork_thickness

            if params.orientation == Orientation.VERTICAL:
                subtree_spacing = max(params.min_subtree_spacing, trunk_width)
                state["size"] = Size(
                    left_info["size"].w
                    + subtree_spacing
                    + right_info["size"].w,
                    subtree_span,
                )
                state["left_pos"] = Position(
                    0,
                    subtree_span - left_info["size"].h,
                )
                state["right_pos"] = Position(
                    left_info["size"].w + subtree_spacing,
                    subtree_span - right_info["size"].h,
                )
                trunk_pos = Position(
                    left_info["size"].w + (subtree_spacing - trunk_width) / 2,
                    0,
                )
            else:
                subtree_spacing = max(params.min_subtree_spacing, trunk_height)
                state["size"] = Size(
                    subtree_span,
                    left_info["size"].h
                    + subtree_spacing
                    + right_info["size"].h,
                )
                state["left_pos"] = Position(
                    subtree_span - left_info["size"].w,
                    0,
                )
                state["right_pos"] = Position(
                    subtree_span - right_info["size"].w,
                    left_info["size"].h + subtree_spacing,
                )
                trunk_pos = Position(
                    0,
                    left_info["size"].h + (subtree_spacing - trunk_height) / 2,
                )

            state["trunk"] = Rect.make_from(trunk_pos, trunk_size)
            state["fork_thickness"] = fork_thickness

    # Compute the absolute position of each subtree
    layout_state[species_tree]["rect"] = Rect.make_from(
        position=Position(0, 0),
        size=layout_state[species_tree]["size"],
    )
    del layout_state[species_tree]["size"]

    for root_species in species_tree.traverse("preorder"):
        this_layout = layout_state[root_species]
        this_rect = this_layout["rect"]

        # Position child subtrees
        if not root_species.is_leaf():
            left_species, right_species = root_species.children

            layout_state[left_species]["rect"] = Rect.make_from(
                position=this_rect.top_left() + this_layout["left_pos"],
                size=layout_state[left_species]["size"],
            )
            del this_layout["left_pos"]
            del layout_state[left_species]["size"]

            layout_state[right_species]["rect"] = Rect.make_from(
                position=this_rect.top_left() + this_layout["right_pos"],
                size=layout_state[right_species]["size"],
            )
            del this_layout["right_pos"]
            del layout_state[right_species]["size"]

        # Make trunk, anchor, and branch nodes positions absolute
        this_layout["trunk"] += this_rect.top_left()
        trunk_rect = this_layout["trunk"]

        for anchor in this_layout["anchors"]:
            if params.orientation == Orientation.VERTICAL:
                this_layout["anchors"][anchor] += trunk_rect.top_right()
            else:
                this_layout["anchors"][anchor] += trunk_rect.bottom_left()

        for branch in this_layout["branches"].values():
            branch["rect"] += trunk_rect.bottom_right()


def _finalize_layout(
    layout_state: Dict[TreeNode, Dict],
    species_tree: Tree,
) -> Layout:
    """Turn a computed layout into final immutable structures."""
    result: Dict[TreeNode, SubtreeLayout] = {}

    for root_species in species_tree.traverse("preorder"):
        this_layout = layout_state[root_species]

        for anchor, branch in this_layout["branches"].items():
            this_layout["branches"][anchor] = Branch(**branch)

        result[root_species] = SubtreeLayout(**layout_state[root_species])

    return result


def compute_layout(
    rec: ReconciliationOutput,
    params: DrawParams = DrawParams(),
) -> Layout:
    """
    Compute the layout of a gene tree embedded in a species tree.

    :param rec: reconciliation object defining the gene and species trees
        their embedding, and an optional synteny labelling
    :param params: layout parameters
    :returns: layout information for each species node
    """
    layout_state: Dict[TreeNode, Dict] = {}
    species_tree = rec.input.species_lca.tree
    _compute_branches(layout_state, rec)
    _layout_branches(layout_state, species_tree, params)
    _layout_subtrees(layout_state, species_tree, params)
    return _finalize_layout(layout_state, species_tree)


def _tikz_draw_fork(  # pylint:disable=too-many-arguments
    species_node: TreeNode,
    layout: SubtreeLayout,
    left_layout: Optional[SubtreeLayout],
    right_layout: Optional[SubtreeLayout],
    layers: Dict[str, List[str]],
    params: DrawParams,
) -> None:
    """Draw the exterior fork of a species subtree."""
    if not species_node.is_leaf():
        # Draw fork
        assert left_layout is not None
        assert right_layout is not None

        if params.orientation == Orientation.VERTICAL:
            left_fork = (
                left_layout.trunk.top_left(),
                "|-",
                layout.trunk.bottom_left(),
                "--",
                layout.trunk.top_left(),
            )
            right_fork = (
                right_layout.trunk.top_right(),
                "|-",
                layout.trunk.bottom_right(),
                "--",
                layout.trunk.top_right(),
            )
            inner_fork = (
                left_layout.trunk.top_right(),
                "|-",
                layout.trunk.bottom_left()
                    + Position(0, layout.fork_thickness),
                "-|",
                right_layout.trunk.top_left(),
            )
        else:
            left_fork = (
                left_layout.trunk.top_left(),
                "-|",
                layout.trunk.top_right(),
                "--",
                layout.trunk.top_left(),
            )
            right_fork = (
                right_layout.trunk.bottom_left(),
                "-|",
                layout.trunk.bottom_right(),
                "--",
                layout.trunk.bottom_left(),
            )
            inner_fork = (
                left_layout.trunk.bottom_left(),
                "-|",
                layout.trunk.bottom_right()
                    + Position(layout.fork_thickness, 0),
                "|-",
                right_layout.trunk.top_left(),
            )

        layers["species"].append(
            rf"""\draw[species border] ({
                left_fork[0]
            }) {left_fork[1]} ({
                left_fork[2]
            }) {left_fork[3]} ({
                left_fork[4]
            });"""
        )
        layers["species"].append(
            rf"""\draw[species border] ({
                right_fork[0]
            }) {right_fork[1]} ({
                right_fork[2]
            }) {right_fork[3]} ({
                right_fork[4]
            });"""
        )
        layers["species"].append(
            rf"""\draw[species border] ({
                inner_fork[0]
            }) {inner_fork[1]} ({
                inner_fork[2]
            }) {inner_fork[3]} ({
                inner_fork[4]
            });"""
        )
    else:
        # Draw leaf
        if params.orientation == Orientation.VERTICAL:
            leaf_shift = Position(0, params.species_leaf_spacing)
            path = (
                layout.trunk.top_left(),
                layout.trunk.bottom_left() + leaf_shift,
                layout.trunk.bottom_right() + leaf_shift,
                layout.trunk.top_right(),
            )
        else:
            leaf_shift = Position(params.species_leaf_spacing, 0)
            path = (
                layout.trunk.top_left(),
                layout.trunk.top_right() + leaf_shift,
                layout.trunk.bottom_right() + leaf_shift,
                layout.trunk.bottom_left(),
            )

        layers["species"].append(
            rf"""\draw[species border] ({
                path[0]
            }) -- ({
                path[1]
            }) -- node[species label] {{{
                tex.escape(species_node.name)
            }}} ({
                path[2]
            }) -- ({
                path[3]
            });"""
        )


def _tikz_draw_branches(  # pylint:disable=too-many-locals,disable=too-many-arguments
    layout: SubtreeLayout,
    left_layout: Optional[SubtreeLayout],
    right_layout: Optional[SubtreeLayout],
    all_layouts: Layout,
    mapping: TreeMapping,
    layers: Dict[str, List[str]],
    params: DrawParams,
) -> None:
    """Draw the interior branches of a species subtree."""
    if params.orientation == Orientation.VERTICAL:
        fork_links = ("|-", "-|")
    else:
        fork_links = ("-|", "|-")

    for root_gene, branch in layout.branches.items():
        branch_pos = branch.rect.center()
        left_gene = branch.left
        right_gene = branch.right

        if root_gene in layout.anchors:
            layers["gene branches"].append(
                rf"""\draw[branch] ({
                    branch_pos
                }) -- ({
                    layout.anchors[root_gene]
                });"""
            )

        if branch.kind == NodeEvent.LEAF:
            if params.orientation == Orientation.VERTICAL:
                leaf_pos = (
                    branch.rect.top()
                    + Position(0, params.extant_gene_diameter / 2)
                )
            else:
                leaf_pos = (
                    branch.rect.left()
                    + Position(params.extant_gene_diameter / 2, 0)
                )

            layers["events"].append(
                rf"\node[extant gene={{{branch.name}}}] "
                rf"at ({leaf_pos}) {{}};"
            )
        elif branch.kind == EdgeEvent.FULL_LOSS:
            if right_gene is None:
                assert left_layout is not None
                keep_pos = left_layout.anchors[left_gene]

                if params.orientation == Orientation.VERTICAL:
                    loss_pos = Position(layout.trunk.right().x, branch_pos.y)
                else:
                    loss_pos = Position(branch_pos.x, layout.trunk.bottom().y)
            else:
                assert right_layout is not None
                keep_pos = right_layout.anchors[right_gene]

                if params.orientation == Orientation.VERTICAL:
                    loss_pos = Position(layout.trunk.left().x, branch_pos.y)
                else:
                    loss_pos = Position(branch_pos.x, layout.trunk.top().y)

            layers["gene branches"].append(
                rf"""\draw[branch] ({branch_pos}) -- ({loss_pos});"""
            )
            layers["events"].append(rf"""\node [loss] at ({loss_pos}) {{}};""")
            layers["gene branches"].append(
                rf"""\draw[branch] ({branch_pos}) {
                    fork_links[1]
                } ({keep_pos});"""
            )
        elif branch.kind == NodeEvent.SPECIATION:
            assert left_layout is not None
            assert right_layout is not None
            layers["gene branches"].append(
                rf"""\draw[branch] ({
                    left_layout.anchors[left_gene]
                }) {fork_links[0]} ({branch_pos}) {fork_links[1]} ({
                    right_layout.anchors[right_gene]
                });"""
            )
            layers["events"].append(
                rf"\node[speciation] at ({branch_pos}) {{{branch.name}}};"
            )
        elif branch.kind == NodeEvent.DUPLICATION:
            layers["gene branches"].append(
                rf"""\draw[branch] ({
                    layout.branches[left_gene].rect.center()
                }) {fork_links[0]} ({branch_pos}) {fork_links[1]} ({
                    layout.branches[right_gene].rect.center()
                });"""
            )
            layers["events"].append(
                rf"\node[duplication] at ({branch_pos}) {{{branch.name}}};"
            )
        elif branch.kind == NodeEvent.HORIZONTAL_TRANSFER:
            foreign_layout = all_layouts[mapping[right_gene]]
            foreign_pos = foreign_layout.anchors[right_gene]

            if params.orientation == Orientation.VERTICAL:
                bend_direction = (
                    "bend left"
                    if branch_pos.x < foreign_pos.x
                    else "bend right"
                )
            else:
                bend_direction = (
                    "bend left"
                    if branch_pos.y > foreign_pos.y
                    else "bend right"
                )

            layers["gene branches"].append(
                rf"""\draw[branch] ({
                    layout.branches[left_gene].rect.center()
                }) |- ({branch_pos});"""
            )
            layers["gene transfers"].append(
                rf"""\draw[transfer branch] ({
                    branch_pos
                }) to[{bend_direction}=35] ({foreign_pos});"""
            )
            layers["events"].append(
                r"\node[horizontal gene transfer] at "
                rf"({branch_pos}) {{{branch.name}}};"
            )
        else:
            raise ValueError("Invalid node type")


def render_to_tikz(
    rec: ReconciliationOutput,
    layout: Layout,
    params: DrawParams = DrawParams(),
):
    r"""
    Generate TikZ code for drawing a laid out reconciliation.

    The `tikz` LaTeX package and the following TikZ libraries are required
    to be loaded for the generated code to compile:

    - `shapes`
    - `arrows.meta`

    Here’s a basic skeleton in which the generated code can be inserted:

    ```
    \documentclass[crop, tikz, border=20pt]{standalone}

    \usepackage{tikz}
    \usetikzlibrary{arrows.meta}
    \usetikzlibrary{shapes}

    \begin{document}
        <generated code>
    \end{document}
    ```

    :param rec: reconciliation object defining the gene and species trees
        their embedding, and an optional synteny labelling
    :param layout: reconciliation layout as computed by :meth:`layout`
    :param params: rendering parameters
    :returns: generated TikZ code
    """
    result = [
        params.get_tikz_definitions(),
        r"\begin{tikzpicture}",
    ]
    layers: Dict[str, List[str]] = {
        "species": [],
        "gene branches": [],
        "gene transfers": [],
        "events": [],
    }

    for species_node in rec.input.species_lca.tree.traverse("preorder"):
        node_layout = layout[species_node]

        if species_node.is_leaf():
            left_layout = None
            right_layout = None
        else:
            left, right = species_node.children
            left_layout = layout[left]
            right_layout = layout[right]

        _tikz_draw_fork(
            species_node, node_layout, left_layout, right_layout, layers, params
        )
        _tikz_draw_branches(
            node_layout,
            left_layout,
            right_layout,
            layout,
            rec.object_species,
            layers,
            params,
        )

    for name, layer in layers.items():
        result.append(f"% {name}")
        result.extend(layer)

    result.append("\\end{tikzpicture}\n")
    return "\n".join(result)
