"""Generate a TikZ drawing from a reconciliation layout."""
from typing import Callable, Dict, List, Sequence, Optional, Tuple
import textwrap
from ete3 import TreeNode
from .model import DrawParams, Orientation, SubtreeLayout, Layout
from ..utils import tex
from ..utils.text import balanced_wrap
from ..model.reconciliation import (
    Event,
    NodeEvent,
    EdgeEvent,
    ReconciliationOutput,
)
from ..model.tree_mapping import TreeMapping
from ..utils.geometry import Position


# Round all coordinates to this number of places in generated TikZ code
MAX_DIGITS = 4


def get_tikz_definitions(params: DrawParams):
    """Get TikZ definitions matching a set of drawing parameters."""
    if params.orientation == Orientation.VERTICAL:
        leaf_label_style = textwrap.dedent(
            r"""
            [font={\strut\color{#1}},
                align=center,
                inner xsep=0pt, inner ysep=2pt,
                outer xsep=0pt, outer ysep=0pt]
            below:#2
            """
        )
        species_label_style = textwrap.dedent(
            rf"""
            font=\bfseries,
            midway,
            anchor=north,
            align=center,
            yshift=-{params.species_label_spacing},
            """
        )
    else:
        leaf_label_style = textwrap.dedent(
            r"""
            [font={\color{#1}},
                align=justify,
                inner xsep=4pt, inner ysep=0pt,
                outer xsep=0pt, outer ysep=0pt]
            right:#2
            """
        )
        species_label_style = textwrap.dedent(
            rf"""
            font=\bfseries,
            midway,
            anchor=west,
            align=left,
            xshift={params.species_label_spacing},
            """
        )

    return textwrap.dedent(
        rf"""
        \colorlet{{species background color}}{{black!15}}
        \tikzset{{
            x={{{params.x_unit}}},
            y={{-{params.y_unit}}},
            species background/.style={{
                fill=species background color,
                draw=species background color,
                line width={{{params.species_border_thickness}}},
            }},
            species label/.style={{
                {textwrap.indent(species_label_style, " " * 16).strip()}
            }},
            branch/.style={{
                draw={{#1}},
                line width={{{params.branch_thickness}}},
            }},
            transfer branch/.style={{
                branch={{#1}},
                -Stealth,
            }},
            loss/.style={{
                draw={{#1}}, cross out, thick,
                line width={{{params.branch_thickness}}},
                inner sep=0pt,
                outer sep=0pt,
                minimum width={{{params.loss_size}}},
                minimum height={{{params.loss_size}}},
            }},
            extant gene/.style 2 args={{
                circle, fill={{#1}},
                outer sep=0pt, inner sep=0pt,
                minimum size={{{params.extant_gene_diameter}}},
                label={{
                    {textwrap.indent(leaf_label_style, " " * 20).strip()}
                }},
            }},
            extant gene/.default={{black}}{{}},
            branch node/.style={{
                draw={{#1}}, fill={{species background color!50!white}},
                align=center,
                font={{\color{{#1}}}},
                outer sep=0pt, inner xsep=0pt, inner ysep=2pt,
                line width={{{params.branch_thickness}}},
            }},
            branch node/.default={{black}},
            speciation/.style={{
                branch node={{#1}}, rectangle, rounded corners,
                inner xsep=4pt,
                minimum width={{{params.speciation_size}}},
                minimum height={{{params.speciation_size}}},
            }},
            duplication/.style={{
                branch node={{#1}}, rectangle,
                inner xsep=4pt,
                minimum width={{{params.duplication_size}}},
                minimum height={{{params.duplication_size}}},
            }},
            horizontal gene transfer/.style={{
                branch node={{#1}}, chamfered rectangle,
                chamfered rectangle sep={{{params.transfer_size} / 2.4}},
                inner xsep=2pt,
                inner ysep=-1pt,
                minimum width={{{params.transfer_size}}},
                minimum height={{{params.transfer_size}}},
            }},
        }}"""
    ).lstrip()


def measure_nodes(
    nodes: Sequence[Tuple[Event, str]],
    params: DrawParams,
) -> Sequence[tex.MeasureBox]:
    """
    Measure the overall space occupied by each node in a set of nodes.

    :param nodes: event nodes to measure
    :param params: drawing settings
    :returns: list of measurements, in the same order as the node sequence
    """
    boxes = []

    for kind, name in nodes:
        if kind == NodeEvent.LEAF:
            boxes.append(rf"\tikz\node[extant gene={{black}}{{{name}}}] {{}};")
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
                left_layout.trunk.top_left().meet_vh(layout.trunk.bottom_left()),
                layout.trunk.bottom_left(),
                layout.trunk.top_left(),
            )
            right_fork = (
                right_layout.trunk.top_right(),
                right_layout.trunk.top_right().meet_vh(layout.trunk.bottom_right()),
                layout.trunk.bottom_right(),
                layout.trunk.top_right(),
            )
            fork_join = layout.trunk.bottom_left() + Position(0, layout.fork_thickness)
            inner_fork = (
                left_layout.trunk.top_right(),
                left_layout.trunk.top_right().meet_vh(fork_join),
                fork_join.meet_hv(right_layout.trunk.top_left()),
                right_layout.trunk.top_left(),
            )
        else:
            left_fork = (
                left_layout.trunk.top_left(),
                left_layout.trunk.top_left().meet_hv(layout.trunk.top_right()),
                layout.trunk.top_right(),
                layout.trunk.top_left(),
            )
            right_fork = (
                right_layout.trunk.bottom_left(),
                right_layout.trunk.bottom_left().meet_hv(layout.trunk.bottom_right()),
                layout.trunk.bottom_right(),
                layout.trunk.bottom_left(),
            )
            fork_join = layout.trunk.bottom_right() + Position(layout.fork_thickness, 0)
            inner_fork = (
                left_layout.trunk.bottom_left(),
                left_layout.trunk.bottom_left().meet_hv(fork_join),
                fork_join.meet_vh(right_layout.trunk.top_left()),
                right_layout.trunk.top_left(),
            )

        layers["background"].append(
            rf"""\path[species background] ({
                left_fork[0] : {MAX_DIGITS}
            }) [rounded corners={{{params.species_border_rounding}}}] -- ({
                left_fork[1] : {MAX_DIGITS}
            }) -- ({
                left_fork[2] : {MAX_DIGITS}
            }) [sharp corners] -- ({
                left_fork[3] : {MAX_DIGITS}
            }) -- ({
                right_fork[3] : {MAX_DIGITS}
            }) [rounded corners={{{params.species_border_rounding}}}] -- ({
                right_fork[2] : {MAX_DIGITS}
            }) -- ({
                right_fork[1] : {MAX_DIGITS}
            }) [sharp corners] -- ({
                right_fork[0] : {MAX_DIGITS}
            }) -- ({
                inner_fork[3] : {MAX_DIGITS}
            }) [rounded corners={{{params.species_border_rounding}}}] -- ({
                inner_fork[2] : {MAX_DIGITS}
            }) -- ({
                inner_fork[1] : {MAX_DIGITS}
            }) [sharp corners] -- ({
                inner_fork[0] : {MAX_DIGITS}
            }) -- cycle;"""
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

        species_name = tex.escape(species_node.name)

        if params.species_label_width is not None:
            species_name = balanced_wrap(
                species_name, params.species_label_width
            ).replace("\n", "\\\\")

        layers["background"].append(
            r"\path[species background, "
            rf"""rounded corners={{{params.species_border_rounding}}}] ({
                path[0] : {MAX_DIGITS}
            }) -- ({
                path[1] : {MAX_DIGITS}
            }) -- node[species label] {{{
                species_name
            }}} ({
                path[2] : {MAX_DIGITS}
            }) -- ({
                path[3] : {MAX_DIGITS}
            });"""
        )


def _tikz_draw_branches(  # pylint:disable=too-many-locals,disable=too-many-arguments
    layout: SubtreeLayout,
    left_layout: Optional[SubtreeLayout],
    right_layout: Optional[SubtreeLayout],
    all_layouts: Layout,
    mapping: TreeMapping,
    layers: Dict[str, List[str]],
    get_color: Callable[str, str],
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
                rf"""\path[branch={{{get_color(branch.color)}}}] ({
                    branch.anchor_parent : {MAX_DIGITS}
                }) -- ({
                    layout.anchors[root_gene] : {MAX_DIGITS}
                });"""
            )

        if branch.kind == NodeEvent.LEAF:
            if params.orientation == Orientation.VERTICAL:
                leaf_pos = branch.rect.top() + Position(
                    0, params.extant_gene_diameter / 2
                )
            else:
                leaf_pos = branch.rect.left() + Position(
                    params.extant_gene_diameter / 2, 0
                )

            layers["events"].append(
                rf"""\node[extant gene={{{
                    get_color(branch.color)
                }}}{{{branch.name}}}] at ({leaf_pos : {MAX_DIGITS}}) {{}};"""
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
                rf"""\path[branch={{{
                    get_color(branch.color)
                }}}] ({branch_pos : {MAX_DIGITS}}) -- ({
                    loss_pos : {MAX_DIGITS}
                });"""
            )
            layers["events"].append(
                rf"""\node[loss={{{
                    get_color(branch.color)
                }}}] at ({loss_pos : {MAX_DIGITS}}) {{}};"""
            )
            layers["gene branches"].append(
                rf"""\path[branch={{{
                    get_color(branch.color)
                }}}] ({branch_pos : {MAX_DIGITS}}) {fork_links[1]} ({
                    keep_pos : {MAX_DIGITS}
                });"""
            )
        elif branch.kind == NodeEvent.SPECIATION:
            assert left_layout is not None
            assert right_layout is not None
            layers["gene branches"].append(
                rf"""\path[branch={{{get_color(branch.color)}}}] ({
                    left_layout.anchors[left_gene] : {MAX_DIGITS}
                }) {fork_links[0]} ({
                    branch.anchor_left : {MAX_DIGITS}
                }) ({
                    branch.anchor_right : {MAX_DIGITS}
                }) {fork_links[1]} ({
                    right_layout.anchors[right_gene] : {MAX_DIGITS}
                });"""
            )
            layers["events"].append(
                rf"""\node[speciation={{{get_color(branch.color)}}}] at ({
                    branch_pos : {MAX_DIGITS}
                }) {{{branch.name}}};"""
            )
        elif branch.kind == NodeEvent.DUPLICATION:
            layers["gene branches"].append(
                rf"""\path[branch={{{get_color(branch.color)}}}] ({
                    layout.branches[left_gene].anchor_parent : {MAX_DIGITS}
                }) {fork_links[0]} ({
                    branch.anchor_left : {MAX_DIGITS}
                }) ({
                    branch.anchor_right : {MAX_DIGITS}
                }) {fork_links[1]} ({
                    layout.branches[right_gene].anchor_parent : {MAX_DIGITS}
                });"""
            )
            layers["events"].append(
                rf"""\node[duplication={{{get_color(branch.color)}}}] at ({
                    branch_pos : {MAX_DIGITS}
                }) {{{branch.name}}};"""
            )
        elif branch.kind == NodeEvent.HORIZONTAL_TRANSFER:
            foreign_layout = all_layouts[mapping[right_gene]]
            foreign_pos = foreign_layout.anchors[right_gene]

            if params.orientation == Orientation.VERTICAL:
                if branch_pos.x < foreign_pos.x:
                    bend_out = "out=0, in=180"
                    anchor_out = branch.anchor_right
                else:
                    bend_out = "out=180, in=0"
                    anchor_out = branch.anchor_left
            else:
                if branch_pos.y > foreign_pos.y:
                    anchor_out = branch.anchor_left
                    bend_out = "out=90, in=-90"
                else:
                    anchor_out = branch.anchor_right
                    bend_out = "out=-90, in=90"

            layers["gene branches"].append(
                rf"""\path[branch={{{get_color(branch.color)}}}] ({
                    layout.branches[left_gene].anchor_parent : {MAX_DIGITS}
                }) |- ({
                    branch.anchor_child : {MAX_DIGITS}
                });"""
            )
            layers["gene transfers"].append(
                rf"""\path[transfer branch={{{get_color(branch.color)}}}] ({
                    anchor_out : {MAX_DIGITS}
                }) to[{bend_out}] ({
                    foreign_pos : {MAX_DIGITS}
                });"""
            )
            # Force content for empty nodes to workaround
            # rendering bug with TikZ chamfered rectangles
            name = branch.name or r"\phantom{-}"
            layers["events"].append(
                rf"""\node[horizontal gene transfer={{{
                    get_color(branch.color)
                }}}] at ({branch_pos : {MAX_DIGITS}}) {{{name}}};"""
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
    layers: Dict[str, List[str]] = {
        "background": [],
        "gene branches": [],
        "gene transfers": [],
        "events": [],
    }
    colors: List[str] = []
    color_prefix = "reccolor"

    def get_color(html: str) -> str:
        if html in colors:
            return f"{color_prefix}{colors.index(html)}"

        colors.append(html)
        return f"{color_prefix}{len(colors) - 1}"

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
            get_color,
            params,
        )

    result = [get_tikz_definitions(params)]

    # Define colors used in the rendering
    for i, html in enumerate(colors):
        result.append(rf"\definecolor{{{color_prefix}{i}}}{{HTML}}{{{html}}}")

    # Append layers in order
    result.append(r"\begin{tikzpicture}")

    for name, layer in layers.items():
        result.append(f"% {name}")
        result.extend(layer)

    result.append(r"\end{tikzpicture}")
    result.append("")
    return "\n".join(result)
