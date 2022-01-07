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
from .model.synteny import SyntenyMapping


class Branch(NamedTuple):
    """Branch in a gene tree embedded in a species tree."""

    # Type of branch, i.e. the event that triggered the gene branching
    kind: Union[NodeEvent, EdgeEvent]

    # Size and position of the branch node relative to the species trunk’s
    # bottom left corner
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

    - If it’s not a leaf node, its left and right subtrees.
    - A fork that joins its left and right subtrees and contains
      the gene copies that are sent to the child species.
    - A trunk that protrudes above the center of the fork and
      contains the branching nodes. At the top of the trunk are
      anchor points that link the subtree to its parent subtree.

    Visual representation of those elements:

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
    # fork and the trunk but also the rects of the child subtrees
    rect: Rect

    # Size and position of the trunk relative to this subtree’s
    # top left corner
    trunk: Rect

    # Number of gene branches in the fork
    fork_thickness: int

    # Anchor points at the top of the trunk relative
    # to the left border
    anchors: Mapping[GeneAnchor, int]

    # Branching nodes in this subtree
    branches: Mapping[GeneAnchor, Branch]

    def get_anchor_pos(self, anchor: GeneAnchor) -> Position:
        """Get the absolute position of an anchor point of this subtree."""
        return (
            self.rect.top_left()
            + self.trunk.top_left()
            + Position(self.anchors[anchor], 0)
        )


Layout = Mapping[TreeNode, SubtreeLayout]


class DrawParams(NamedTuple):
    """Parameters for drawing a reconciliation."""

    # Horizontal unit scale (unit-less parameters are multiples of this unit)
    x_unit: str = "1pt"

    # Vertical unit scale (unit-less parameters are multiples of this unit)
    y_unit: str = "1pt"

    # Minimum space between the outline of the species tree
    # and one of the gene branches it contains
    species_branch_padding: float = 10

    # Minimum space between two gene branches in the tree
    gene_branch_spacing: float = 5

    # Space above trunks
    trunk_overhead: float = 10

    # Horizontal space between the trunk of a parent node
    # and its children subtrees
    subtree_spacing: float = 0

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

    # Length of the dashed line that represent lost genes
    full_loss_size: float = 20

    # Minimum size of the hollow circles that represent speciation events
    speciation_size: float = 10

    # Minimum size of the square that represent duplication events
    duplication_size: float = 10

    # Minimum size of the sideways square that represent transfer events
    transfer_size: float = 10

    def get_tikz_definitions(self):
        """Get TikZ definitions matching a set of drawing parameters."""
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
                    yshift=-{self.species_label_spacing},
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
                    branch, dashed,
                }},
                extant gene/.style={{
                    circle, fill,
                    outer sep=0pt, inner sep=0pt,
                    minimum size={{{self.extant_gene_diameter}}},
                    label={{
                        [font={{\strut}},
                         inner xsep=0pt,
                         inner ysep=2pt,
                         outer ysep=0pt,
                         fill=white]
                        below:#1
                    }},
                }},
                branch node/.style={{
                    draw, fill=white,
                    outer sep=0pt, inner sep=0pt,
                    line width={{{self.branch_thickness}}},
                    font={{\strut}},
                }},
                speciation/.style={{
                    branch node, rounded rectangle,
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
                    name = "".join(syntenies[root_gene])
                else:
                    species_name, gene_name = root_gene.name.split("_")
                    name = rf"{species_name}\textsubscript{{{gene_name}}}"

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
                    rf"{''.join(syntenies[root_gene])}"
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
    """Compute branching node sizes."""
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
        next_pos_x: float = 0
        next_pos_y = params.species_branch_padding
        layout = layout_state[root_species]
        layout["anchors"] = {}

        for root_gene, branch in layout["branches"].items():
            size = (
                measures[root_gene].overall_size()
                if root_gene in measures
                else Size(0, 0)
            )

            if branch["kind"] == NodeEvent.LEAF:
                pos = Position(next_pos_x, -size.h)
                next_pos_x += params.gene_branch_spacing + size.w
            elif branch["kind"] in (
                NodeEvent.SPECIATION,
                EdgeEvent.FULL_LOSS,
            ):
                pos = Position(next_pos_x, next_pos_y)
                next_pos_x += params.gene_branch_spacing + size.w
                next_pos_y += params.gene_branch_spacing + size.h
            elif branch["kind"] == NodeEvent.DUPLICATION:
                left_rect = layout["branches"][branch["left"]]["rect"]
                right_rect = layout["branches"][branch["right"]]["rect"]
                pos = Position(
                    x=((left_rect.center() + right_rect.center()).x) / 2
                    - size.w / 2,
                    y=min(
                        params.species_branch_padding, left_rect.y, right_rect.y
                    )
                    - params.species_branch_padding
                    - size.h,
                )
            elif branch["kind"] == NodeEvent.HORIZONTAL_TRANSFER:
                cons_rect = layout["branches"][branch["left"]]["rect"]
                pos = Position(
                    x=cons_rect.center().x - size.w / 2,
                    y=min(params.species_branch_padding, cons_rect.y)
                    - params.species_branch_padding
                    - size.h,
                )
            else:
                raise ValueError("Invalid node type")

            rect = Rect.make_from(pos, size)
            branch["rect"] = rect

            if root_gene in layout["anchor_nodes"]:
                layout["anchors"][root_gene] = rect.center().x

        del layout["anchor_nodes"]

        # Shift all nodes to the right to make room for the left padding
        if layout["branches"]:
            padding_shift = Position(
                x=(
                    -min(
                        branch["rect"].x
                        for branch in layout["branches"].values()
                    )
                    + params.species_branch_padding
                ),
                y=0,
            )

            for root_gene in layout["branches"]:
                branch = layout["branches"][root_gene]
                branch["rect"] += padding_shift
                layout["branches"][root_gene] = Branch(**branch)

            for root_gene in layout["anchors"]:
                layout["anchors"][root_gene] += padding_shift.x


def _layout_subtrees(
    layout_state: Dict[TreeNode, Dict],
    species_tree: Tree,
    params: DrawParams,
) -> Layout:
    """Compute the size and absolute position of each subtree."""
    # Compute the size of each subtree
    for root_species in species_tree.traverse("postorder"):
        state = layout_state[root_species]

        if state["branches"]:
            trunk_width = (
                max(
                    branch.rect.bottom_right().x
                    for branch in state["branches"].values()
                )
                + params.species_branch_padding
            )
            trunk_height = (
                max(
                    (0,)
                    + tuple(
                        -branch.rect.top_left().y
                        for branch in state["branches"].values()
                    )
                )
                + params.trunk_overhead
            )
            fork_thickness = (
                max(
                    (0,)
                    + tuple(
                        branch.rect.bottom_left().y
                        for branch in state["branches"].values()
                    )
                )
                + params.species_branch_padding
            )
        else:
            # Empty subtree
            trunk_width = 2 * params.species_branch_padding
            trunk_height = params.trunk_overhead
            fork_thickness = params.species_branch_padding

        if root_species.is_leaf():
            # Extant species
            state["size"] = Size(trunk_width, trunk_height)
            state["trunk"] = Rect(0, 0, trunk_width, trunk_height)
            state["fork_thickness"] = 0
        else:
            # Ancestral species
            left_species, right_species = root_species.children
            left_info = layout_state[left_species]
            right_info = layout_state[right_species]

            height = (
                max(left_info["size"].h, right_info["size"].h)
                + params.level_spacing
                + fork_thickness
                + trunk_height
            )

            state["size"] = Size(
                left_info["size"].w
                + params.subtree_spacing * 2
                + trunk_width
                + right_info["size"].w,
                height,
            )

            state["left_pos"] = Position(0, height - left_info["size"].h)
            state["right_pos"] = Position(
                left_info["size"].w + params.subtree_spacing * 2 + trunk_width,
                height - right_info["size"].h,
            )
            state["trunk"] = Rect(
                left_info["size"].w + params.subtree_spacing,
                0,
                trunk_width,
                trunk_height,
            )
            state["fork_thickness"] = fork_thickness

    # Compute the absolute position of each subtree
    result: Dict[TreeNode, SubtreeLayout] = {}
    layout_state[species_tree]["rect"] = Rect.make_from(
        position=Position(0, 0),
        size=layout_state[species_tree]["size"],
    )
    del layout_state[species_tree]["size"]

    for root_species in species_tree.traverse("preorder"):
        if not root_species.is_leaf():
            left_species, right_species = root_species.children
            this_rect = layout_state[root_species]["rect"]

            layout_state[left_species]["rect"] = Rect.make_from(
                position=this_rect.top_left()
                + layout_state[root_species]["left_pos"],
                size=layout_state[left_species]["size"],
            )
            del layout_state[root_species]["left_pos"]
            del layout_state[left_species]["size"]

            layout_state[right_species]["rect"] = Rect.make_from(
                position=this_rect.top_left()
                + layout_state[root_species]["right_pos"],
                size=layout_state[right_species]["size"],
            )
            del layout_state[root_species]["right_pos"]
            del layout_state[right_species]["size"]

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
    return _layout_subtrees(layout_state, species_tree, params)


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

        layers["species"].append(
            rf"""\draw[species border] ({
                left_layout.rect.top_left() + left_layout.trunk.top_left()
            }) |- ({
                layout.rect.top_left()
                    + layout.trunk.bottom_left()
            }) -- ({
                layout.rect.top_left() + layout.trunk.top_left()
            });"""
        )
        layers["species"].append(
            rf"""\draw[species border] ({
                right_layout.rect.top_left()
                    + right_layout.trunk.top_right()
            }) |- ({
                layout.rect.top_left()
                    + layout.trunk.bottom_right()
            }) -- ({
                layout.rect.top_left()
                    + layout.trunk.top_right()
            });"""
        )
        layers["species"].append(
            rf"""\draw[species border] ({
                left_layout.rect.top_left() + left_layout.trunk.top_right()
            }) |- ({
                layout.rect.top_left()
                    + layout.trunk.bottom_left()
                    + Position(0, layout.fork_thickness)
            }) -| ({
                right_layout.rect.top_left()
                    + right_layout.trunk.top_left()
            });"""
        )
    else:
        # Draw leaf
        leaf_shift = Position(0, params.species_leaf_spacing)
        layers["species"].append(
            rf"""\draw[species border] ({
                layout.rect.top_left() + layout.trunk.top_left()
            }) -- ({
                layout.rect.top_left()
                    + layout.trunk.bottom_left()
                    + leaf_shift
            }) -- node[species label] {{{species_node.name}}} ({
                layout.rect.top_left()
                    + layout.trunk.bottom_right()
                    + leaf_shift
            }) -- ({
                layout.rect.top_left() + layout.trunk.top_right()
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
    trunk_offset = layout.rect.top_left() + layout.trunk.bottom_left()

    for root_gene, branch in layout.branches.items():
        branch_pos = trunk_offset + branch.rect.center()

        left_gene = branch.left
        right_gene = branch.right

        if root_gene in layout.anchors:
            layers["gene branches"].append(
                rf"""\draw[branch] ({
                    branch_pos
                }) -- ({
                    layout.get_anchor_pos(root_gene)
                });"""
            )

        if branch.kind == NodeEvent.LEAF:
            leaf_pos = (
                trunk_offset
                + branch.rect.top()
                + Position(0, params.extant_gene_diameter / 2)
            )
            layers["events"].append(
                rf"\node[extant gene={{{branch.name}}}] "
                rf"at ({leaf_pos}) {{}};"
            )
        elif branch.kind == EdgeEvent.FULL_LOSS:
            if right_gene is None:
                assert left_layout is not None
                loss_pos = f"{params.full_loss_size}, 0"
                keep_pos = left_layout.get_anchor_pos(left_gene)
            else:
                assert right_layout is not None
                loss_pos = f"-{params.full_loss_size}, 0"
                keep_pos = right_layout.get_anchor_pos(right_gene)

            layers["gene branches"].append(
                rf"""\draw[loss] ({branch_pos}) -- ++({loss_pos});"""
            )
            layers["gene branches"].append(
                rf"""\draw[branch] ({branch_pos}) -| ({keep_pos});"""
            )
        elif branch.kind == NodeEvent.SPECIATION:
            assert left_layout is not None
            assert right_layout is not None
            layers["gene branches"].append(
                rf"""\draw[branch] ({
                    left_layout.get_anchor_pos(left_gene)
                }) |- ({branch_pos}) -| ({
                    right_layout.get_anchor_pos(right_gene)
                });"""
            )
            layers["events"].append(
                rf"\node[speciation] at ({branch_pos}) {{{branch.name}}};"
            )
        elif branch.kind == NodeEvent.DUPLICATION:
            layers["gene branches"].append(
                rf"""\draw[branch] ({
                    trunk_offset
                        + layout.branches[left_gene].rect.center()
                }) |- ({branch_pos}) -| ({
                    trunk_offset
                        + layout.branches[right_gene].rect.center()
                });"""
            )
            layers["events"].append(
                rf"\node[duplication] at ({branch_pos}) {{{branch.name}}};"
            )
        elif branch.kind == NodeEvent.HORIZONTAL_TRANSFER:
            foreign_layout = all_layouts[mapping[right_gene]]
            foreign_pos = foreign_layout.get_anchor_pos(right_gene)
            bend_direction = (
                "bend left" if branch_pos.x < foreign_pos.x else "bend right"
            )

            layers["gene branches"].append(
                rf"""\draw[branch] ({
                    trunk_offset
                        + layout.branches[left_gene].rect.center()
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
