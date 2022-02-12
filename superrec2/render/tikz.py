"""Generate a TikZ drawing from a reconciliation layout."""
from typing import Dict, List, Sequence, Optional, Tuple
import textwrap
from ete3 import TreeNode
from .model import DrawParams, Orientation, SubtreeLayout, Layout
from ..utils import tex
from ..model.reconciliation import (
    Event, NodeEvent, EdgeEvent, ReconciliationOutput
)
from ..model.tree_mapping import TreeMapping
from ..utils.geometry import Position


def get_tikz_definitions(params: DrawParams):
    """Get TikZ definitions matching a set of drawing parameters."""
    if params.orientation == Orientation.VERTICAL:
        leaf_label_style = textwrap.dedent(
            r"""
            [font={{\strut}},
                fill=white,
                inner xsep=0pt, inner ysep=2pt,
                outer xsep=0pt, outer ysep=0pt]
            below:#1
            """
        )
        species_label_style = textwrap.dedent(
            rf"""
            font=\bfseries,
            midway,
            anchor=north,
            yshift=-{params.species_label_spacing}
            """
        )
    else:
        leaf_label_style = textwrap.dedent(
            r"""
            [fill=white,
                inner xsep=4pt, inner ysep=0pt,
                outer xsep=0pt, outer ysep=0pt]
            right:#1
            """
        )
        species_label_style = textwrap.dedent(
            rf"""
            font=\bfseries,
            midway,
            anchor=west,
            xshift={params.species_label_spacing}
            """
        )

    return textwrap.dedent(
        rf"""
        \tikzset{{
            x={{{params.x_unit}}},
            y={{-{params.y_unit}}},
            species border/.style={{
                line width={{{params.species_border_thickness}}},
                shorten <={{-{params.species_border_thickness} / 2 + 0.05pt}},
                shorten >={{-{params.species_border_thickness} / 2 + 0.05pt}},
            }},
            species label/.style={{
                {textwrap.indent(species_label_style, " " * 16).strip()}
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
                draw, cross out, thick,
                line width={{{params.branch_thickness}}},
                inner sep=0pt,
                outer sep=0pt,
                minimum width={{{params.loss_size}}},
                minimum height={{{params.loss_size}}},
            }},
            extant gene/.style={{
                circle, fill,
                outer sep=0pt, inner sep=0pt,
                minimum size={{{params.extant_gene_diameter}}},
                label={{
                    {textwrap.indent(leaf_label_style, " " * 20).strip()}
                }},
            }},
            branch node/.style={{
                draw, fill=white,
                outer sep=0pt, inner xsep=0pt, inner ysep=2pt,
                line width={{{params.branch_thickness}}},
            }},
            speciation/.style={{
                branch node, rounded rectangle,
                inner xsep=4pt,
                minimum width={{{params.speciation_size}}},
                minimum height={{{params.speciation_size}}},
            }},
            duplication/.style={{
                branch node, rectangle,
                inner xsep=4pt,
                minimum width={{{params.duplication_size}}},
                minimum height={{{params.duplication_size}}},
            }},
            horizontal gene transfer/.style={{
                branch node, signal, signal to=east and west,
                minimum width={{{params.transfer_size}}},
                minimum height={{{params.transfer_size}}},
            }},
        }}"""
    ).lstrip()


def measure_nodes(
    nodes: Sequence[Tuple[Event, str]],
    params: DrawParams,
) -> Sequence[tex.MeasureBox]:
    boxes = []

    for kind, name in nodes:
        if kind == NodeEvent.LEAF:
            boxes.append(rf"\tikz\node[extant gene={{{name}}}] {{}};")
        else:
            if kind == NodeEvent.SPECIATION:
                node_type = "[speciation]"
            elif kind == NodeEvent.DUPLICATION:
                node_type = "[duplication]"
            elif kind == NodeEvent.HORIZONTAL_TRANSFER:
                node_type = "[horizontal gene transfer]"
            elif kind == EdgeEvent.FULL_LOSS:
                node_type = "[loss]"
            else:
                raise ValueError("Invalid node type")

            boxes.append(rf"\tikz\node{node_type} {{{name}}};")

    return tex.measure(
        boxes,
        preamble=(
            r"\usepackage{tikz}"
            r"\usetikzlibrary{arrows.meta}"
            r"\usetikzlibrary{shapes}" + get_tikz_definitions(params)
        ),
    )


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


def render(
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
        get_tikz_definitions(params),
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
