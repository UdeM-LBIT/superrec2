from ete3 import PhyloTree, PhyloNode
from enum import Enum, auto
import textwrap
from typing import Dict, Mapping, NamedTuple
from .reconciliation import Reconciliation
from .utils import invert_mapping, sort_tree_nodes


class Position(NamedTuple):
    x: int
    y: int


class Size(NamedTuple):
    w: int
    h: int


class Rect(NamedTuple):
    x: int
    y: int
    w: int
    h: int


class BranchKind(Enum):
    Speciation = auto()
    Duplication = auto()


class Branch(NamedTuple):
    # Branch event type
    kind: BranchKind

    # Position of the branching node relative to
    # the upper left corner of the branching origin
    pos: Position

    # Left and right children nodes
    left: PhyloNode
    right: PhyloNode


class SubtreeLayout(NamedTuple):
    # Size of this subtree
    size: Size

    # Relative position of the left child’s subtree
    left_pos: Position

    # Relative position of the right child’s subtree
    right_pos: Position

    # Size and relative position of the upper trunk
    trunk_rect: Rect

    # Thickness of the branching node
    branch_thickness: int

    # Anchor points at the top of the trunk relative
    # to the left of the trunk
    anchors: Mapping[PhyloNode, int]

    # Branching nodes at the root of this subtree
    branches: Mapping[PhyloNode, Branch]


Layout = Mapping[PhyloNode, SubtreeLayout]


class LayoutParams(NamedTuple):
    branch_spacing: int = 20
    subtree_spacing: int = 20
    level_spacing: int = 20


def layout(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
    rec: Reconciliation,
    params: LayoutParams = LayoutParams()
) -> Layout:
    result: Dict[PhyloNode, SubtreeLayout] = {}
    rev_rec = invert_mapping(rec)

    for root_species in species_tree.traverse("postorder"):
        if root_species.is_leaf():
            # Extant species
            width = (
                params.branch_spacing
                * (len(rev_rec[root_species]) + 1)
            )

            result[root_species] = SubtreeLayout(
                size=Size(width, 0),
                left_pos=Position(0, 0),
                right_pos=Position(0, 0),
                trunk_rect=Rect(0, 0, width, 0),
                branch_thickness=0,
                anchors={
                    gene_node: (i + 1) * params.branch_spacing
                    for i, gene_node in enumerate(rev_rec[root_species])
                },
                branches=[],
            )
        else:
            # Ancestral species
            left_species, right_species = root_species.children
            left_info = result[left_species]
            right_info = result[right_species]

            # Find branching nodes that belong to this species’ ancestry
            trunk_width = params.branch_spacing
            trunk_height = 0

            anchors: Dict[PhyloNode, int] = {}
            branches: Dict[PhyloNode, Branch] = {}

            for root_gene in sort_tree_nodes(rev_rec[root_species]):
                left_gene, right_gene = root_gene.children

                # Speciation node
                if (
                    rec[left_gene] == left_species
                    and rec[right_gene] == right_species
                ):
                    offset = trunk_width
                    anchors[root_gene] = offset
                    trunk_width += params.branch_spacing

                    branches[root_gene] = Branch(
                        kind=BranchKind.Speciation,
                        pos=Position(offset, offset),
                        left=left_gene,
                        right=right_gene,
                    )

                # Duplication node
                if (
                    rec[left_gene] == root_species
                    and rec[right_gene] == root_species
                ):
                    y_pos = -trunk_height
                    trunk_height += params.branch_spacing

                    left_pos = branches[left_gene].pos
                    right_pos = branches[right_gene].pos

                    del anchors[left_gene]
                    del anchors[right_gene]
                    anchors[root_gene] = (left_pos.x + right_pos.x) / 2

                    branches[root_gene] = Branch(
                        kind=BranchKind.Duplication,
                        pos=Position(anchors[root_gene], y_pos),
                        left=left_gene,
                        right=right_gene,
                    )

            height = (
                max(left_info.size.h, right_info.size.h)
                + params.level_spacing
                + trunk_width
                + trunk_height
            )

            result[root_species] = SubtreeLayout(
                size=Size(
                    left_info.size.w + params.subtree_spacing * 2
                        + trunk_width + right_info.size.w,
                    height,
                ),
                left_pos=Position(0, height - left_info.size.h),
                right_pos=Position(
                    left_info.size.w + params.subtree_spacing * 2
                        + trunk_width,
                    height - right_info.size.h,
                ),
                trunk_rect=Rect(
                    left_info.size.w + params.subtree_spacing, 0,
                    trunk_width, trunk_height,
                ),
                branch_thickness=trunk_width,
                anchors=anchors,
                branches=branches,
            )

    return result


class RenderParams(NamedTuple):
    x_unit: str = "1pt"
    y_unit: str = "-1pt"
    species_border_thickness: str = "1pt"
    branch_thickness: str = "0.5pt"
    branch_outer_thickness: str = "4pt"
    extant_gene_diameter: str = "3pt"
    speciation_diameter: str = "10pt"
    duplication_size: str = "9pt"


def render(
    species_tree: PhyloTree,
    layout: Layout,
    params: RenderParams = RenderParams(),
):
    result = []
    result.append(textwrap.dedent(rf"""
        \begin{{tikzpicture}}[
            x={{{params.x_unit}}},
            y={{{params.y_unit}}},
            species border/.style={{
                line width={{{params.species_border_thickness}}},
            }},
            branch/.style={{
                line width={{{params.branch_thickness}}},
                preaction={{
                    draw, white,
                    line width={{{params.branch_outer_thickness}}},
                }},
            }},
            extant gene/.style={{
                circle, fill,
                outer sep=0pt, inner sep=0pt,
                minimum width={{{params.extant_gene_diameter}}},
                minimum height={{{params.extant_gene_diameter}}},
                label={{below:#1}},
            }},
            branch node/.style={{
                draw, fill=white,
                outer sep=0pt, inner sep=0pt,
                line width={{{params.branch_thickness}}},
            }},
            speciation/.style={{
                branch node, circle,
                minimum width={{{params.speciation_diameter}}},
                minimum height={{{params.speciation_diameter}}},
            }},
            duplication/.style={{
                branch node, rectangle,
                minimum width={{{params.duplication_size}}},
                minimum height={{{params.duplication_size}}},
            }},
        ]
    """))

    pos: Dict[PhyloNode, Position] = {}
    pos[species_tree] = Position(0, 0)

    for species_node in species_tree.traverse("preorder"):
        x, y = pos[species_node]
        node_layout = layout[species_node]

        if not species_node.is_leaf():
            left, right = species_node.children
            pos[left] = Position(
                x + node_layout.left_pos.x,
                y + node_layout.left_pos.y
            )
            pos[right] = Position(
                x + node_layout.right_pos.x,
                y + node_layout.right_pos.y
            )
            left_layout = layout[left]
            right_layout = layout[right]

            # Draw outer fork
            result.append(rf"""\draw[species border] ({
                pos[left].x + left_layout.trunk_rect.x}, {
                pos[left].y + left_layout.trunk_rect.y
            }) |- ({
                x + node_layout.trunk_rect.x}, {
                y + node_layout.trunk_rect.y + node_layout.trunk_rect.h
            }) -- ({
                x + node_layout.trunk_rect.x}, {
                y + node_layout.trunk_rect.y
            }) ({
                x + node_layout.trunk_rect.x + node_layout.trunk_rect.w}, {
                y + node_layout.trunk_rect.y
            }) -- ({
                x + node_layout.trunk_rect.x + node_layout.trunk_rect.w}, {
                y + node_layout.trunk_rect.y + node_layout.trunk_rect.h
            }) -| ({
                pos[right].x + right_layout.trunk_rect.x
                    + right_layout.trunk_rect.w}, {
                pos[right].y + right_layout.trunk_rect.y
            });""")

            # Draw inner fork
            result.append(rf"""\draw[species border] ({
                pos[left].x + left_layout.trunk_rect.x
                    + left_layout.trunk_rect.w}, {
                pos[left].y + left_layout.trunk_rect.y
            }) |- ({
                pos[right].x + right_layout.trunk_rect.x}, {
                y + node_layout.trunk_rect.y + node_layout.trunk_rect.h
                    + node_layout.branch_thickness
            }) -- ({
                pos[right].x + right_layout.trunk_rect.x}, {
                pos[right].y + right_layout.trunk_rect.y
            });""")

            # Draw branches
            branching_nodes = []
            offset_x = x + node_layout.trunk_rect.x
            offset_y = y + node_layout.trunk_rect.y \
                + node_layout.trunk_rect.h

            for root_gene, branch in node_layout.branches.items():
                abs_x = offset_x + branch.pos.x
                abs_y = offset_y + branch.pos.y

                left_gene = branch.left
                right_gene = branch.right

                if root_gene in node_layout.anchors:
                    result.append(rf"""\draw[branch] ({abs_x}, {abs_y}) -- ({
                        offset_x + node_layout.anchors[root_gene]}, {
                        offset_y - node_layout.trunk_rect.h
                    });""")

                if branch.kind == BranchKind.Speciation:
                    left_anchor = left_layout.anchors[left_gene]
                    right_anchor = right_layout.anchors[right_gene]

                    result.append(rf"""\draw[branch] ({
                        pos[left].x + left_layout.trunk_rect.x
                            + left_anchor}, {
                        pos[left].y + left_layout.trunk_rect.y
                    }) |- ({abs_x}, {abs_y}) -| ({
                        pos[right].x + right_layout.trunk_rect.x
                            + right_anchor}, {
                        pos[right].y + right_layout.trunk_rect.y
                    });""")
                    branching_nodes.append(
                        rf"\node[speciation] at ({abs_x}, {abs_y}) {{}};"
                    )
                elif branch.kind == BranchKind.Duplication:
                    left_pos = node_layout.branches[left_gene].pos
                    right_pos = node_layout.branches[right_gene].pos

                    result.append(rf"""\draw[branch] ({
                        offset_x + left_pos.x}, {
                        offset_y + left_pos.y
                    }) |- ({abs_x}, {abs_y}) -| ({
                        offset_x + right_pos.x}, {
                        offset_y + right_pos.y
                    });""")
                    branching_nodes.append(
                        rf"\node[duplication] at ({abs_x}, {abs_y}) {{}};"
                    )

            result.extend(branching_nodes)
        else:
            # Draw species
            result.append(
                r"\draw[species border] "
                rf"({x}, {y}) -- ({x + node_layout.size.w}, {y});"
            )

            for gene, anchor_pos in node_layout.anchors.items():
                result.append(
                    rf"\node[extant gene={{\({gene.name}\)}}] "
                    rf"at ({x + anchor_pos}, {y}) {{}};"
                )

    result.append("\end{tikzpicture}")
    return "\n".join(result)
