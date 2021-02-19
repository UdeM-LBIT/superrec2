from ete3 import PhyloTree, PhyloNode
from enum import Enum, auto
import textwrap
from typing import Dict, Mapping, NamedTuple, Union
from .geometry import Position, Rect, Size
from .lowest_common_ancestor import LowestCommonAncestor
from .reconciliation import Event, get_event, Reconciliation
from .utils import invert_mapping, sort_tree_nodes


class BranchKind(Enum):
    FullLoss = auto()
    Speciation = auto()
    Duplication = auto()
    HorizontalGeneTransfer = auto()


class Branch(NamedTuple):
    # Branch event type
    kind: BranchKind

    # Position of the branching node relative to
    # the bottom left corner of the trunk
    pos: Position

    # Left and right children nodes
    left: PhyloNode
    right: PhyloNode


class PseudoNode():
    pass


AnchorNode = Union[PhyloNode, PseudoNode]


class SubtreeLayout(NamedTuple):
    # Size and position of this subtree
    rect: Rect

    # Size and relative position of the trunk
    trunk: Rect

    # Thickness of the subtree branches
    branch_thickness: int

    # Anchor points at the top of the trunk relative
    # to the left of the trunk
    anchors: Mapping[AnchorNode, int]

    # Branching nodes at the root of this subtree
    branches: Mapping[AnchorNode, Branch]

    def get_anchor_pos(self, anchor: AnchorNode):
        return (
            self.rect.top_left()
            + self.trunk.top_left()
            + Position(self.anchors[anchor], 0)
        )


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
    work_state: Dict[PhyloNode, Dict] = {}
    result: Dict[PhyloNode, SubtreeLayout] = {}

    species_lca = LowestCommonAncestor(species_tree)
    rev_rec = invert_mapping(rec)

    def add_losses(gene, start_species, end_species):
        prev_gene = gene
        prev_species = start_species
        start_species = start_species.up

        while start_species != end_species:
            is_left = prev_species == start_species.children[0]
            is_right = prev_species == start_species.children[1]
            state = work_state[start_species]
            cur_gene = PseudoNode()

            x_pos = state["next_x"]
            state["next_x"] += params.branch_spacing

            state["anchors"][cur_gene] = x_pos
            state["branches"][cur_gene] = Branch(
                kind=BranchKind.FullLoss,
                pos=Position(x_pos, x_pos),
                left=prev_gene if is_left else None,
                right=prev_gene if is_right else None,
            )

            prev_gene = cur_gene
            prev_species = start_species
            start_species = start_species.up

        return prev_gene

    # First pass: Compute the branching nodes for each species
    for root_species in species_tree.traverse("postorder"):
        if not root_species.is_leaf():
            left_species, right_species = root_species.children

            work_state[root_species] = state = {
                "branches": {},
                "anchors": {},
                "next_x": params.branch_spacing,
            }

            for root_gene in sort_tree_nodes(rev_rec[root_species]):
                left_gene, right_gene = root_gene.children
                event = get_event(root_gene, species_lca, rec)

                if event == Event.Speciation:
                    left_gene = add_losses(
                        left_gene, rec[left_gene], root_species)
                    right_gene = add_losses(
                        right_gene, rec[right_gene], root_species)

                    x_pos = state["next_x"]
                    state["next_x"] += params.branch_spacing
                    state["anchors"][root_gene] = x_pos

                    state["branches"][root_gene] = Branch(
                        kind=BranchKind.Speciation,
                        pos=Position(x_pos, x_pos),
                        left=left_gene,
                        right=right_gene,
                    )
                elif event == Event.Duplication:
                    left_gene = add_losses(
                        left_gene, rec[left_gene], root_species.up)
                    right_gene = add_losses(
                        right_gene, rec[right_gene], root_species.up)

                    left_pos = state["branches"][left_gene].pos
                    right_pos = state["branches"][right_gene].pos
                    x_pos = (left_pos.x + right_pos.x) / 2
                    y_pos = (
                        min(params.branch_spacing, left_pos.y, right_pos.y)
                        - params.branch_spacing
                    )
                    state["anchors"][root_gene] = x_pos

                    del state["anchors"][left_gene]
                    del state["anchors"][right_gene]

                    state["branches"][root_gene] = Branch(
                        kind=BranchKind.Duplication,
                        pos=Position(x_pos, y_pos),
                        left=left_gene,
                        right=right_gene,
                    )
                elif event == Event.HorizontalGeneTransfer:
                    conserv_gene, foreign_gene = (
                        (left_gene, right_gene)
                        if species_lca.is_ancestor_of(
                            root_species, rec[left_gene]
                        )
                        else (right_gene, left_gene)
                    )
                    conserv_gene = add_losses(
                        conserv_gene, rec[conserv_gene], root_species.up)
                    conserv_pos = state["branches"][conserv_gene].pos

                    x_pos = conserv_pos.x
                    y_pos = (
                        min(params.branch_spacing, conserv_pos.y)
                        - params.branch_spacing
                    )
                    state["anchors"][root_gene] = x_pos
                    del state["anchors"][conserv_gene]

                    state["branches"][root_gene] = Branch(
                        kind=BranchKind.HorizontalGeneTransfer,
                        pos=Position(x_pos, y_pos),
                        left=conserv_gene,
                        right=foreign_gene,
                    )
                else:
                    raise ValueError("Invalid event")

    # Second pass: Compute the individual size and layout of each subtree
    for root_species in species_tree.traverse("postorder"):
        if root_species.is_leaf():
            width = (
                params.branch_spacing
                * (len(rev_rec[root_species]) + 1)
            )

            work_state[root_species] = {
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
            state = work_state[root_species]

            trunk_width = state["next_x"]
            trunk_height = (
                max(-branch.pos.y for branch in state["branches"].values())
                + params.branch_spacing
            )

            del state["next_x"]

            left_info = work_state[left_species]
            right_info = work_state[right_species]

            height = (
                max(left_info["size"].h, right_info["size"].h)
                + params.level_spacing
                + trunk_width
                + trunk_height
            )

            state["size"] = Size(
                left_info["size"].w + params.subtree_spacing * 2
                    + trunk_width + right_info["size"].w,
                height,
            )

            state["left_pos"] = Position(0, height - left_info["size"].h)
            state["right_pos"] = Position(
                left_info["size"].w + params.subtree_spacing * 2
                    + trunk_width,
                height - right_info["size"].h,
            )
            state["trunk"] = Rect(
                left_info["size"].w + params.subtree_spacing, 0,
                trunk_width, trunk_height,
            )
            state["branch_thickness"] = trunk_width

    # Third pass: Compute the absolute position of each subtree
    work_state[species_tree]["rect"] = Rect.make_from(
        position=Position(0, 0),
        size=work_state[species_tree]["size"],
    )
    del work_state[species_tree]["size"]

    for root_species in species_tree.traverse("preorder"):
        if not root_species.is_leaf():
            left_species, right_species = root_species.children
            this_rect = work_state[root_species]["rect"]

            work_state[left_species]["rect"] = Rect.make_from(
                position=this_rect.top_left()
                    + work_state[root_species]["left_pos"],
                size=work_state[left_species]["size"],
            )
            del work_state[root_species]["left_pos"]
            del work_state[left_species]["size"]

            work_state[right_species]["rect"] = Rect.make_from(
                position=this_rect.top_left()
                    + work_state[root_species]["right_pos"],
                size=work_state[right_species]["size"],
            )
            del work_state[root_species]["right_pos"]
            del work_state[right_species]["size"]

        result[root_species] = SubtreeLayout(**work_state[root_species])

    return result


class RenderParams(NamedTuple):
    x_unit: str = "1pt"
    y_unit: str = "-1pt"
    species_border_thickness: str = "1pt"
    branch_thickness: str = "0.5pt"
    branch_outer_thickness: str = "4pt"
    extant_gene_diameter: str = "3pt"
    full_loss_size: str = "20pt"
    speciation_diameter: str = "10pt"
    duplication_size: str = "9pt"


def render_to_tikz(
    species_tree: PhyloTree,
    rec: Reconciliation,
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
                shorten <={{-{params.species_border_thickness} / 2 + 0.05pt}},
                shorten >={{-{params.species_border_thickness} / 2 + 0.05pt}},
            }},
            branch/.style={{
                line width={{{params.branch_thickness}}},
                preaction={{
                    draw, white, -{{}},
                    line width={{{params.branch_outer_thickness}}},
                    shorten <={{{params.branch_thickness}}},
                    shorten >={{{params.branch_thickness}}},
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
            horizontal gene transfer/.style={{
                branch node, diamond,
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
            });""")
            result.append(rf"""\draw[species border] ({
                right_layout.rect.top_left() + right_layout.trunk.top_right()
            }) |- ({
                node_layout.rect.top_left() + node_layout.trunk.bottom_right()
            }) -- ({
                node_layout.rect.top_left() + node_layout.trunk.top_right()
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
                        node_layout.get_anchor_pos(root_gene)
                    });""")

                if branch.kind == BranchKind.FullLoss:
                    if right_gene is None:
                        result.append(
                            rf"\draw[loss] ({branch_pos}) -- "
                            rf"++({params.full_loss_size}, 0);"
                        )
                        result.append(rf"""\draw[branch] ({branch_pos}) -| ({
                            left_layout.get_anchor_pos(left_gene)
                        });""")
                    elif left_gene is None:
                        result.append(
                            rf"\draw[loss] ({branch_pos}) -- "
                            rf"++(-{params.full_loss_size}, 0);"
                        )
                        result.append(rf"""\draw[branch] ({branch_pos}) -| ({
                            right_layout.get_anchor_pos(right_gene)
                        });""")
                elif branch.kind == BranchKind.Speciation:
                    result.append(rf"""\draw[branch] ({
                        left_layout.get_anchor_pos(left_gene)
                    }) |- ({branch_pos}) -| ({
                        right_layout.get_anchor_pos(right_gene)
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
                elif branch.kind == BranchKind.HorizontalGeneTransfer:
                    foreign_layout = layout[rec[right_gene]]
                    result.append(rf"""\draw[branch] ({
                        trunk_offset + node_layout.branches[left_gene].pos
                    }) |- ({branch_pos});""")
                    result.append(rf"""\draw[transfer branch] ({
                        branch_pos
                    }) to[bend left=45] ({
                        foreign_layout.get_anchor_pos(right_gene)
                    });""")
                    branching_nodes.append(
                        r"\node[horizontal gene transfer] at "
                        rf"({branch_pos}) {{}};"
                    )

            result.extend(branching_nodes)
        else:
            # Draw species
            result.append(rf"""\draw[species border] ({
                node_layout.rect.top_left()
            }) -- ({
                node_layout.rect.top_right()
            });""")

            for gene in node_layout.anchors:
                result.append(rf"""\node[extant gene={{\({gene.name}\)}}] at ({
                    node_layout.get_anchor_pos(gene)
                }) {{}};""")

    result.append("\end{tikzpicture}")
    return "\n".join(result)
