from ete3 import PhyloTree, PhyloNode
from enum import Enum, auto
import textwrap
from typing import Dict, Mapping, NamedTuple, Optional, Union
from .geometry import Position, Rect, Size
from .lowest_common_ancestor import LowestCommonAncestor
from .reconciliation import Event, get_event, Reconciliation
from .utils import invert_mapping, sort_tree_nodes


class BranchKind(Enum):
    """Type of branch event."""
    # Not actually a branch but a leaf node
    Leaf = auto()

    # Loss of the parent gene in one of the two children species
    FullLoss = auto()

    # Transmission of the parent gene to both children species
    Speciation = auto()

    # Duplication of the parent gene in the same genome
    Duplication = auto()

    # Transfer of the parent gene to a foreign genome
    HorizontalGeneTransfer = auto()


class Branch(NamedTuple):
    """Branch in a gene tree embedded in a species tree."""
    # Type of branch, i.e. the event that triggered the gene branching
    kind: BranchKind

    # Branch node position relative to the species trunk’s bottom left corner
    pos: Position

    # Children genes
    # * If this is a full loss branch, either `left` or `right` will be None,
    #   which corresponds to the lost gene copy.
    # * If this is a horizontal gene transfer branch, `left` is the conserved
    #   copy and `right` is the transfered copy.
    # * If this is a leaf branch, both children will be None
    left: Optional[PhyloNode]
    right: Optional[PhyloNode]


class PseudoGene():
    """Objects used as virtual nodes for lost genes."""
    pass


GeneAnchor = Union[PhyloNode, PseudoGene]


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


Layout = Mapping[PhyloNode, SubtreeLayout]


class LayoutParams(NamedTuple):
    """Parameters for laying out subtrees of a reconciliation."""
    # Space between the outline of the species tree and the gene branches
    species_branch_padding: int = 20

    # Space between two gene branches in the tree
    gene_branch_spacing: int = 15

    # Space above trunks
    trunk_overhead: int = 20

    # Horizontal space between the trunk of a parent node
    # and its children subtrees
    subtree_spacing: int = 0

    # Vertical space between the fork of the parent node
    # and its children subtrees
    level_spacing: int = 0


def _layout_compute_branches(
    layout_state: Dict[PhyloNode, Dict],
    species_tree: PhyloTree,
    rec: Reconciliation,
    rev_rec: Mapping,
    params: LayoutParams
) -> None:
    """Compute the branching nodes for each species."""
    species_lca = LowestCommonAncestor(species_tree)

    def add_losses(gene, start_species, end_species):
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

            x_pos = state["next_x"]
            state["next_x"] += params.gene_branch_spacing

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

    for root_species in species_tree.traverse("postorder"):
        layout_state[root_species] = state = {
            "branches": {},
            "anchors": {},
            "next_x": params.species_branch_padding,
        }

        # Create branches even for leaf genes so that duplications
        # above leaves can be properly represented
        for extant_gene in rev_rec[root_species]:
            if extant_gene.is_leaf():
                x_pos = state["next_x"]
                state["next_x"] += params.gene_branch_spacing
                state["anchors"][extant_gene] = x_pos
                state["branches"][extant_gene] = Branch(
                    kind=BranchKind.Leaf,
                    pos=Position(x_pos, 0),
                    left=None,
                    right=None,
                )

        # Create branches for actual internal nodes
        internal_genes = sort_tree_nodes([
            gene for gene in rev_rec[root_species]
            if not gene.is_leaf()
        ])

        for root_gene in internal_genes:
            left_gene, right_gene = root_gene.children
            event = get_event(root_gene, species_lca, rec)

            if event == Event.Speciation:
                # Speciation nodes are located below the trunk
                # and linked to child species’s gene anchors
                left_gene = add_losses(
                    left_gene, rec[left_gene], root_species)
                right_gene = add_losses(
                    right_gene, rec[right_gene], root_species)

                x_pos = state["next_x"]
                state["next_x"] += params.gene_branch_spacing
                state["anchors"][root_gene] = x_pos
                state["branches"][root_gene] = Branch(
                    kind=BranchKind.Speciation,
                    pos=Position(x_pos, x_pos),
                    left=left_gene,
                    right=right_gene,
                )
            elif event == Event.Duplication:
                # Duplications are located in the trunk and linked
                # to other nodes in the same species
                left_gene = add_losses(
                    left_gene, rec[left_gene], root_species.up)
                right_gene = add_losses(
                    right_gene, rec[right_gene], root_species.up)

                left_pos = state["branches"][left_gene].pos
                right_pos = state["branches"][right_gene].pos
                x_pos = (left_pos.x + right_pos.x) / 2
                y_pos = (
                    min(
                        params.gene_branch_spacing,
                        left_pos.y,
                        right_pos.y
                    ) - params.gene_branch_spacing
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
                # Transfers are located in the trunk, like duplications,
                # but are linked to a node outside the current subtree
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
                    min(params.gene_branch_spacing, conserv_pos.y)
                    - params.gene_branch_spacing
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


def _layout_compute_sizes(
    layout_state: Dict[PhyloNode, Dict],
    species_tree: PhyloTree,
    rev_rec: Mapping,
    params: LayoutParams
):
    """Compute the individual size and layout of each subtree."""
    for root_species in species_tree.traverse("postorder"):
        state = layout_state[root_species]

        if state["branches"]:
            trunk_width = (
                state["next_x"]
                - params.gene_branch_spacing
                + params.species_branch_padding
            )
        else:
            # Empty subtree
            trunk_width = (
                state["next_x"]
                + params.species_branch_padding
            )

        trunk_height = max((0,) + tuple(
            -branch.pos.y
            for branch in state["branches"].values()
        )) + params.trunk_overhead

        del state["next_x"]

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
            state["fork_thickness"] = trunk_width


def _layout_compute_positions(
    layout_state: Dict[PhyloNode, Dict],
    species_tree: PhyloTree,
) -> Layout:
    """Compute the absolute position of each subtree."""
    result: Layout = {}
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


def layout(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
    rec: Reconciliation,
    params: LayoutParams = LayoutParams()
) -> Layout:
    """
    Compute the layout of a gene tree embedded in a species tree.

    :param gene_tree: embedded gene tree
    :param species_tree: host species tree
    :param rec: mapping of the gene tree onto the species tree
    :param params: layout parameters
    :returns: layout information for each species node
    """
    layout_state: Dict[PhyloNode, Dict] = {}
    rev_rec = invert_mapping(rec)

    _layout_compute_branches(layout_state, species_tree, rec, rev_rec, params)
    _layout_compute_sizes(layout_state, species_tree, rev_rec, params)
    return _layout_compute_positions(layout_state, species_tree)


class RenderParams(NamedTuple):
    x_unit: str = "1pt"
    y_unit: str = "1pt"
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
            y={{-{params.y_unit}}},
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
            # Draw fork
            left, right = species_node.children
            left_layout = layout[left]
            right_layout = layout[right]

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
            result.append(rf"""\draw[species border] ({
                left_layout.rect.top_left() + left_layout.trunk.top_right()
            }) |- ({
                node_layout.rect.top_left() + node_layout.trunk.bottom_left()
                    + Position(0, node_layout.fork_thickness)
            }) -| ({
                right_layout.rect.top_left() + right_layout.trunk.top_left()
            });""")
        else:
            # Draw leaf
            result.append(rf"""\draw[species border] ({
                node_layout.rect.top_left() + node_layout.trunk.top_left()
            }) -- ({
                node_layout.rect.top_left() + node_layout.trunk.bottom_left()
            }) -- ({
                node_layout.rect.top_left() + node_layout.trunk.bottom_right()
            }) -- ({
                node_layout.rect.top_left() + node_layout.trunk.top_right()
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

            if branch.kind == BranchKind.Leaf:
                branching_nodes.append(
                    rf"\node[extant gene={{\({root_gene.name}\)}}] "
                    rf"at ({branch_pos}) {{}};"
                )
            elif branch.kind == BranchKind.FullLoss:
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
                foreign_pos = foreign_layout.get_anchor_pos(right_gene)
                bend_direction = (
                    "bend left" if branch_pos.x < foreign_pos.x
                    else "bend right"
                )

                result.append(rf"""\draw[branch] ({
                    trunk_offset + node_layout.branches[left_gene].pos
                }) |- ({branch_pos});""")
                result.append(rf"""\draw[transfer branch] ({
                    branch_pos
                }) to[{bend_direction}=35] ({foreign_pos});""")
                branching_nodes.append(
                    r"\node[horizontal gene transfer] at "
                    rf"({branch_pos}) {{}};"
                )

        result.extend(branching_nodes)

    result.append("\end{tikzpicture}")
    return "\n".join(result)
