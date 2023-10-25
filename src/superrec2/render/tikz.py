"""Generate a TikZ drawing from a reconciliation layout."""
from typing import Sequence
from dataclasses import dataclass
import textwrap
import math
from .model import DrawParams, Orientation, HostLayout, Layout
from ..model.synteny import format_synteny
from ..model.history import Event, Extant, Codiverge, Diverge, Gain, Loss
from ..utils import tex
from ..utils.tex import measure
from ..utils.text import balanced_wrap
from ..utils.geometry import Position, Size


# Round all coordinates to this number of decimal places in generated TikZ code
MAX_DIGITS = 4


def get_tikz_definitions(params: DrawParams):
    """Get TikZ definitions matching a set of drawing parameters."""
    if params.orientation == Orientation.Vertical:
        leaf_label_style = textwrap.dedent(
            r"""
            [font={\color{#1}},
                align=center,
                text depth=0pt,
                inner xsep=0pt, inner ysep=4pt,
                outer xsep=0pt, outer ysep=0pt]
            below:#2
            """
        )
        host_label_style = textwrap.dedent(
            rf"""
            font=\bfseries,
            midway,
            anchor=north,
            align=center,
            yshift=-{params.host_label_spacing},
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
        host_label_style = textwrap.dedent(
            rf"""
            font=\bfseries,
            midway,
            anchor=west,
            align=left,
            xshift={params.host_label_spacing},
            """
        )

    return textwrap.dedent(
        rf"""
        \colorlet{{host background color}}{{black!15}}
        \tikzset{{
            x={{{params.x_unit}}},
            y={{-{params.y_unit}}},
            host background/.style={{
                fill=host background color,
                draw=host background color,
                line width={{{params.host_border_thickness}}},
            }},
            unsampled host background/.style={{
                pattern={{Lines[angle=0, distance=2.5pt, line width=1pt]}},
                pattern color=host background color!65,
            }},
            host label/.style={{
                {textwrap.indent(host_label_style, " " * 16).strip()}
            }},
            branch/.style={{
                draw={{#1}},
                line width={{{params.branch_thickness}}},
            }},
            transfer branch/.style={{
                branch={{#1}},
                densely dashed,
                -{{Latex[length=0pt 8]}},
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
            event/.style={{
                draw={{#1}}, fill={{host background color!50!white}},
                align=center,
                font={{\color{{#1}}}},
                outer sep=0pt, inner xsep=0pt, inner ysep=2pt,
                line width={{{params.branch_thickness}}},
            }},
            event/.default={{black}},
            speciation/.style={{
                event={{#1}}, rectangle, rounded corners,
                inner xsep=4pt,
                minimum width={{{params.speciation_size}}},
                minimum height={{{params.speciation_size}}},
            }},
            duplication/.style={{
                event={{#1}}, rectangle,
                inner xsep=4pt,
                minimum width={{{params.duplication_size}}},
                minimum height={{{params.duplication_size}}},
            }},
            transfer/.style={{
                event={{#1}}, chamfered rectangle,
                chamfered rectangle sep={{{params.transfer_size} / 2.4}},
                inner xsep=2pt,
                inner ysep=-1pt,
                minimum width={{{params.transfer_size}}},
                minimum height={{{params.transfer_size}}},
            }},
            gain/.style 2 args={{
                draw={{#1}}, fill=white,
                circle, inner sep=1.5pt,
                label={{[font={{\strut\scriptsize}},xshift=-2pt]right:#2}},
            }},
            loss/.style={{
                draw={{#1}}, cross out, thick,
                line width={{{params.branch_thickness}}},
                inner sep=0pt,
                outer sep=0pt,
                minimum width={{{params.loss_size}}},
                minimum height={{{params.loss_size}}},
            }},
        }}"""
    ).lstrip()


def render_event(event: Event, position: Position, params: DrawParams) -> str:
    """Generate the TikZ code for drawing an event node."""
    if event.contents is not None:
        label = format_synteny(
            map(tex.escape, event.contents),
            params.event_label_width,
        ).replace("\n", "\\\\")
    elif event.name is not None and isinstance(event, Extant):
        if "_" in event.name:
            base, index = event.name.rsplit("_", 1)
            label = tex.escape(base.replace("_", " "))
            label += rf"\textsubscript{{{tex.escape(index)}}}"
        else:
            label = tex.escape(event.name)
    else:
        label = ""

    match event:
        case Extant():
            return rf"\node[extant gene={{black}}{{{label}}}] at " f"({position}) {{}};"

        case Codiverge():
            return rf"\node[speciation] at ({position}) {{{label}}};"

        case Diverge(segment, _, transfer):
            kind = "transfer" if transfer else "duplication"
            return rf"\node[{kind}] at ({position}) {{{label}}};"

        case Gain(gained):
            return rf"\node[gain={{{gained or ''}}}{{{label}}}] at ({position}) {{}};"

        case Loss(segment):
            return rf"\node[loss={{{segment or ''}}}] at ({position}) {{{label}}};"

    return rf"\node at ({position}) {{{label}{segment}}};"


def measure_events(events: Sequence[Event], params: DrawParams) -> dict[Event, Size]:
    """
    Measure the overall space occupied by each node in a set of nodes.

    :param nodes: event nodes to measure
    :param params: drawing settings
    :returns: list of measurements, in the same order as the node sequence
    """
    events = list(events)
    results = measure(
        texts=(
            r"\tikz" + render_event(event, Position(0, 0), params) for event in events
        ),
        preamble=(
            r"\usepackage{tikz}"
            r"\usetikzlibrary{arrows.meta}"
            r"\usetikzlibrary{shapes}" + get_tikz_definitions(params)
        ),
    )
    return {event: result.overall_size() for event, result in zip(events, results)}


@dataclass(frozen=True)
class PathNode:
    coords: Position
    style: str | None = None
    node: str | None = None


def _tikz_path(path: Sequence[PathNode], rounding: int, close: bool) -> str:
    """Generate TikZ code for drawing a sequence of nodes."""
    result = ""

    for node in path:
        coords_str = f"({node.coords:{rounding}})"

        if not result:
            result += coords_str
        else:
            style_str = f"[{node.style}] " if node.style is not None else ""
            node_str = f"{node.node} " if node.node is not None else ""
            result += f" {style_str}-- {node_str}{coords_str}"

    if close:
        result += " -- cycle"

    return result


def _render_host(
    layout: HostLayout,
    children: list[HostLayout],
    sampled: bool,
    params: DrawParams,
) -> str:
    """Draw the fork connecting a host node to its children."""
    rounded = f"rounded corners={{{params.host_border_radius}}}"
    sharp = "sharp corners"

    if children:
        # Connect root with leftmost and rightmost children
        # TODO: Connect to intermediate children
        fork_area = layout.fork_events_area
        left_area = children[0].events_area
        right_area = children[-1].events_area
        path = (
            PathNode(left_area.top_left()),
            PathNode(left_area.left().meet_vh(fork_area.top()), rounded),
            PathNode(fork_area.top_left()),
            PathNode(layout.events_area.top_left(), sharp),
            PathNode(layout.events_area.top_right()),
            PathNode(fork_area.top_right(), rounded),
            PathNode(right_area.right().meet_vh(fork_area.top())),
            PathNode(right_area.top_right(), sharp),
            PathNode(right_area.top_left()),
            PathNode(right_area.left().meet_vh(fork_area.bottom()), rounded),
            PathNode(left_area.right().meet_vh(fork_area.bottom())),
            PathNode(left_area.top_right(), sharp),
        )
    else:
        leaf_shift = Position(0, params.host_leaf_spacing)

        if sampled:
            host_name = tex.escape(layout.host.name)

            if params.host_label_width is not None:
                host_name = balanced_wrap(host_name, params.host_label_width).replace(
                    "\n", "\\\\"
                )

            host_label = f"node[host label] {{{host_name}}}"
        else:
            host_label = ""

        path = (
            PathNode(layout.area.top_left()),
            PathNode(layout.area.bottom_left() + leaf_shift, rounded),
            PathNode(layout.area.bottom_right() + leaf_shift, node=host_label),
            PathNode(layout.area.top_right(), sharp),
        )

    style = "host background" if sampled else "unsampled host background"
    return rf"\path[{style}] {_tikz_path(path, MAX_DIGITS, close=True)};"


def _render_branch(
    start: Position, end: Position, transfer: bool, params: DrawParams
) -> str:
    """Draw a branch connecting two locations of the associate tree."""
    rounded = f"rounded corners={{{params.branch_border_radius}}}"
    sharp = "sharp corners"

    if math.isclose(start.x, end.x):
        path = (PathNode(start), PathNode(end))
        return rf"\path[branch] {_tikz_path(path, MAX_DIGITS, close=True)};"

    midpoint = start.meet_hv(end)

    if transfer:
        tr_path = (PathNode(start), PathNode(midpoint))
        end_path = (PathNode(midpoint), PathNode(end))
        return (
            rf"\path[transfer branch] {_tikz_path(tr_path, MAX_DIGITS, close=False)};"
            rf"\path[branch] {_tikz_path(end_path, MAX_DIGITS, close=False)};"
        )

    path = (PathNode(start), PathNode(midpoint, rounded), PathNode(end, sharp))
    return rf"\path[branch] {_tikz_path(path, MAX_DIGITS, close=False)};"


def render(
    layout: Layout,
    params: DrawParams = DrawParams(),
) -> str:
    r"""
    Generate TikZ code for drawing a history using layout information.

    The `tikz` LaTeX package and the following TikZ libraries are required
    to be loaded for the generated code to compile:

    - `shapes`
    - `arrows.meta`

    Here’s a basic skeleton in which the generated code can be inserted:

    ```
    \documentclass[crop, tikz, border=20pt]{standalone}

    \usepackage{tikz}
    \usetikzlibrary{patterns.meta}
    \usetikzlibrary{arrows.meta}
    \usetikzlibrary{shapes}

    \begin{document}
        <generated code>
    \end{document}
    ```

    :param layout: layout to render
    :param params: rendering parameters
    :returns: generated TikZ code
    """
    layers: dict[str, list[str]] = {
        "hosts": [],
        "branches": [],
        "gene transfers": [],
        "events": [],
    }
    colors: list[str] = []
    color_prefix = "reccolor"

    def get_color(html: str) -> str:
        if html in colors:
            return f"{color_prefix}{colors.index(html)}"

        colors.append(html)
        return f"{color_prefix}{len(colors) - 1}"

    for host, host_layout in layout.items():
        layers["hosts"].append(
            _render_host(
                host_layout,
                [layout[child] for child in host_layout.children],
                sampled=host_layout.host.sampled,
                params=params,
            )
        )

        for event_node, event_layout in host_layout.events.items():
            event = event_node.data

            for child in event_layout.in_children + event_layout.out_children:
                layers["branches"].append(
                    _render_branch(
                        event_layout.area.center(),
                        layout[child.data.host].events[child].area.center(),
                        transfer=(isinstance(event, Diverge) and event.transfer),
                        params=params,
                    )
                )

            layers["events"].append(
                render_event(event_node.data, event_layout.area.center(), params)
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
