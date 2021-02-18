from ete3 import PhyloTree, PhyloNode
from enum import Enum, auto
import textwrap
from typing import Dict, Mapping, NamedTuple
from .geometry import Position, Rect, Size
from .reconciliation import Reconciliation
from .utils import invert_mapping, sort_tree_nodes


class BranchKind(Enum):
    Speciation = auto()
    Duplication = auto()


class Branch(NamedTuple):
    # Branch event type
    kind: BranchKind

    # Position of the branching node relative to
    # the bottom left corner of the trunk
    pos: Position

    # Left and right children nodes
    left: PhyloNode
    right: PhyloNode


class SubtreeLayout(NamedTuple):
    # Size and position of this subtree
    rect: Rect

    # Size and relative position of the trunk
    trunk: Rect

    # Thickness of the subtree branches
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
    first_pass: Dict[PhyloNode, Dict] = {}
    final_pass: Dict[PhyloNode, SubtreeLayout] = {}
    rev_rec = invert_mapping(rec)

    # Compute the individual size and layout of each subtree
    for root_species in species_tree.traverse("postorder"):
        if root_species.is_leaf():
            # Extant species
            width = (
                params.branch_spacing
                * (len(rev_rec[root_species]) + 1)
            )

            first_pass[root_species] = {
                "size": Size(width, 0),
                "trunk": Rect(0, 0, width, 0),
                "branch_thickness": 0,
                "anchors": {
                    gene_node: (i + 1) * params.branch_spacing
                    for i, gene_node in enumerate(rev_rec[root_species])
                },
                "branches": [],
            }
        else:
            # Ancestral species
            left_species, right_species = root_species.children
            left_info = first_pass[left_species]
            right_info = first_pass[right_species]

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
                max(left_info["size"].h, right_info["size"].h)
                + params.level_spacing
                + trunk_width
                + trunk_height
            )

            first_pass[root_species] = {
                "size": Size(
                    left_info["size"].w + params.subtree_spacing * 2
                        + trunk_width + right_info["size"].w,
                    height,
                ),
                "left_pos": Position(0, height - left_info["size"].h),
                "right_pos": Position(
                    left_info["size"].w + params.subtree_spacing * 2
                        + trunk_width,
                    height - right_info["size"].h,
                ),
                "trunk": Rect(
                    left_info["size"].w + params.subtree_spacing, 0,
                    trunk_width, trunk_height,
                ),
                "branch_thickness": trunk_width,
                "anchors": anchors,
                "branches": branches,
            }

    # Compute the absolute position of each subtree
    first_pass[species_tree]["rect"] = Rect.make_from(
        position=Position(0, 0),
        size=first_pass[species_tree]["size"],
    )
    del first_pass[species_tree]["size"]

    for root_species in species_tree.traverse("preorder"):
        if not root_species.is_leaf():
            left_species, right_species = root_species.children
            this_rect = first_pass[root_species]["rect"]

            first_pass[left_species]["rect"] = Rect.make_from(
                position=this_rect.top_left()
                    + first_pass[root_species]["left_pos"],
                size=first_pass[left_species]["size"],
            )
            del first_pass[root_species]["left_pos"]
            del first_pass[left_species]["size"]

            first_pass[right_species]["rect"] = Rect.make_from(
                position=this_rect.top_left()
                    + first_pass[root_species]["right_pos"],
                size=first_pass[right_species]["size"],
            )
            del first_pass[root_species]["right_pos"]
            del first_pass[right_species]["size"]

        final_pass[root_species] = SubtreeLayout(**first_pass[root_species])

    return final_pass


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

    for species_node in species_tree.traverse("preorder"):
        node_layout = layout[species_node]

        if not species_node.is_leaf():
            left, right = species_node.children
            left_layout = layout[left]
            right_layout = layout[right]

            # Draw outer fork
            result.append(rf"""\draw[species border] ({
                left_layout.rect.top_left() + left_layout.trunk.top_left()
            }) |- ({
                node_layout.rect.top_left() + node_layout.trunk.bottom_left()
            }) -- ({
                node_layout.rect.top_left() + node_layout.trunk.top_left()
            }) ({
                node_layout.rect.top_left() + node_layout.trunk.top_right()
            }) -- ({
                node_layout.rect.top_left() + node_layout.trunk.bottom_right()
            }) -| ({
                right_layout.rect.top_left() + right_layout.trunk.bottom_right()
            });""")

            # Draw inner fork
            result.append(rf"""\draw[species border] ({
                left_layout.rect.top_left() + left_layout.trunk.top_right()
            }) |- ({
                node_layout.rect.top_left() + node_layout.trunk.bottom_left()
                    + Position(0, node_layout.branch_thickness)
            }) -| ({
                right_layout.rect.top_left() + right_layout.trunk.top_left()
            });""")

            # Draw branches
            branching_nodes = []
            trunk_offset = (
                node_layout.rect.top_left()
                + node_layout.trunk.bottom_left()
            )

            for root_gene, branch in node_layout.branches.items():
                branch_pos = trunk_offset + branch.pos

                left_gene = branch.left
                right_gene = branch.right

                if root_gene in node_layout.anchors:
                    result.append(rf"""\draw[branch] ({branch_pos}) -- ({
                        node_layout.rect.top_left()
                            + node_layout.trunk.top_left()
                            + Position(node_layout.anchors[root_gene], 0)
                    });""")

                if branch.kind == BranchKind.Speciation:
                    result.append(rf"""\draw[branch] ({
                        left_layout.rect.top_left()
                            + left_layout.trunk.top_left()
                            + Position(left_layout.anchors[left_gene], 0)
                    }) |- ({branch_pos}) -| ({
                        right_layout.rect.top_left()
                            + right_layout.trunk.top_left()
                            + Position(right_layout.anchors[right_gene], 0)
                    });""")
                    branching_nodes.append(
                        rf"\node[speciation] at ({branch_pos}) {{}};"
                    )
                elif branch.kind == BranchKind.Duplication:
                    result.append(rf"""\draw[branch] ({
                        trunk_offset + node_layout.branches[left_gene].pos
                    }) |- ({branch_pos}) -| ({
                        trunk_offset + node_layout.branches[right_gene].pos
                    });""")
                    branching_nodes.append(
                        rf"\node[duplication] at ({branch_pos}) {{}};"
                    )

            result.extend(branching_nodes)
        else:
            # Draw species
            result.append(rf"""\draw[species border] ({
                node_layout.rect.top_left()
            }) -- ({
                node_layout.rect.top_right()
            });""")

            for gene, anchor_pos in node_layout.anchors.items():
                result.append(rf"""\node[extant gene={{\({gene.name}\)}}] at ({
                    node_layout.rect.top_left() + Position(anchor_pos, 0)
                }) {{}};""")

    result.append("\end{tikzpicture}")
    return "\n".join(result)
