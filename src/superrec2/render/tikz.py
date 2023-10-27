"""Generate a TikZ drawing from a reconciliation layout."""
from typing import Callable, Sequence
from dataclasses import dataclass
from textwrap import indent, dedent
import math
import re
from .model import DrawParams, Orientation, HostLayout, Layout
from ..model.history import Event, Extant, Codiverge, Diverge, Gain, Loss
from ..utils import tex
from ..utils.tex import measure
from ..utils.geometry import Position, Size


# Regex to match digit groups in a string
DIGITS = re.compile(r"([0-9]+)")


# Round all coordinates to this number of decimal places in generated TikZ code
MAX_DIGITS = 4


def get_tikz_definitions(params: DrawParams):
    """Get TikZ definitions matching a set of drawing parameters."""
    style = []

    if params.orientation == Orientation.Vertical:
        style.extend(
            (
                f"x={{({params.x_unit},0)}}",
                f"y={{(0,-{params.y_unit})}}",
            )
        )
    else:
        style.extend(
            (
                f"y={{({params.x_unit},0)}}",
                f"x={{(0,-{params.y_unit})}}",
            )
        )

    style.extend(
        (
            """\
            limit width/.style 2 args={
                execute at begin node={\\begin{varwidth}{#1}#2},
                execute at end node={\\end{varwidth}},
            }""",
            f"""\
            host background/.style={{
                fill=host background color,
                draw=host background color,
                line width={{{params.host_border_thickness}}},
            }}""",
            """\
            unsampled host background/.style={
                pattern={Lines[angle=45, distance=3pt, line width=2pt]},
                pattern color=host background color!50,
            }""",
        )
    )

    if params.orientation == Orientation.Vertical:
        host_label_anchor = "north"
        host_label_align = "\\centering"
    else:
        host_label_anchor = "west"
        host_label_align = "\\narrowragged"

    style.extend(
        (
            f"""\
            host label/.style={{
                font=\\bfseries,
                midway,
                anchor={host_label_anchor},
                limit width={{{params.host_label_width}em}}{{{host_label_align}}},
            }}""",
            f"""\
            branch/.style={{
                draw={{#1}},
                line width={{{params.branch_thickness}}},
            }}""",
            """\
            transfer branch/.style={
                branch={#1},
                densely dashed,
                -{Latex[length=0pt 8]},
            }""",
            f"""\
            event/.style={{
                draw={{#1}}, fill={{host background color!50!white}},
                font={{\\color{{#1}}}},
                outer sep=0pt, inner xsep=0pt, inner ysep=2pt,
                limit width={{{params.event_label_width}em}}{{\\centering}},
                line width={{{params.branch_thickness}}},
            }}""",
            "event/.default={black}",
        )
    )

    if params.orientation == Orientation.Vertical:
        extant_label_position = "below"
        extant_label_align = "\\centering"
    else:
        extant_label_position = "right"
        extant_label_align = "\\narrowragged"

    style.extend(
        (
            f"""\
            extant/.style 2 args={{
                circle, fill={{#1}},
                outer sep=0pt, inner sep=0pt,
                minimum size={{{params.extant_gene_diameter}}},
                label={{
                    [font={{\\color{{#1}}}},
                        text depth=0pt,
                        inner xsep=0pt, inner ysep=4pt,
                        outer xsep=0pt, outer ysep=0pt,
                        limit width={{{params.event_label_width}em}}{{{
                            extant_label_align}}}]
                    {extant_label_position}:#2
                }},
            }}""",
            "extant/.default={black}{}",
            f"""\
            speciation/.style={{
                event={{#1}}, rectangle, rounded corners,
                inner xsep=4pt,
                minimum width={{{params.speciation_size}}},
                minimum height={{{params.speciation_size}}},
            }}""",
            f"""\
            duplication/.style={{
                event={{#1}}, rectangle,
                inner xsep=4pt,
                minimum width={{{params.duplication_size}}},
                minimum height={{{params.duplication_size}}},
            }}""",
            f"""\
            transfer/.style={{
                event={{#1}}, chamfered rectangle,
                chamfered rectangle sep={{{params.transfer_size} / 2.4}},
                inner xsep=2pt,
                inner ysep=-1pt,
                minimum width={{{params.transfer_size}}},
                minimum height={{{params.transfer_size}}},
            }}""",
        )
    )

    if params.orientation == Orientation.Vertical:
        segment_position = "right"
        segment_align = "\\narrowragged"
    else:
        segment_position = "above"
        segment_align = "\\centering"

    style.extend(
        (
            f"""\
            gain/.style 2 args={{
                draw={{#1}}, fill=white,
                circle, inner sep=1.5pt,
                label={{[%
                    font={{\\strut\\scriptsize}},
                    limit width={{{params.event_label_width}em}}{{{segment_align}}},
                ]{segment_position}:#2}},
            }}""",
            f"""\
            loss/.style 2 args={{
                draw={{#1}}, cross out, thick,
                line width={{{params.branch_thickness}}},
                inner sep=0pt,
                outer sep=0pt,
                minimum width={{{params.loss_size}}},
                minimum height={{{params.loss_size}}},
                label={{[%
                    font={{\\strut\\scriptsize}},
                    limit width={{{params.event_label_width}em}}{{{segment_align}}},
                ]{segment_position}:#2}},
            }}""",
        )
    )

    result = "\\colorlet{host background color}{black!15}\n"
    complete_style = ",\n".join(indent(dedent(item), " " * 4) for item in style)
    result += f"\\tikzset{{\n{complete_style},\n}}"
    return result


def contents_key(item: str):
    parts = DIGITS.split(item)
    return [int(part) if part.isdigit() else part for part in parts]


def format_contents(
    contents: frozenset[str] | None,
    wrapper: Callable[[str], str] = lambda x: x,
) -> str:
    return ", ".join(
        wrapper(tex.escape(item)) for item in sorted(contents, key=contents_key)
    )


def render_event(event: Event, position: Position, params: DrawParams) -> str:
    """Generate the TikZ code for drawing an event node."""
    if event.contents is not None:
        label = format_contents(event.contents)
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
            return rf"\node[extant={{black}}{{{label}}}] at ({position}) {{}};"

        case Codiverge():
            return rf"\node[speciation] at ({position}) {{{label}}};"

        case Diverge(contents=contents, segment=segment, transfer=transfer, cut=cut):
            kind = "transfer" if transfer else "duplication"

            if contents is not None:
                if cut:
                    remainder = format_contents(contents - segment)
                    segment = format_contents(segment)
                    label = f"{segment} / {remainder or '---'}"
                else:
                    label = ", ".join(
                        [
                            rf"\underline{{{tex.escape(item)}}}"
                            for item in sorted(segment, key=contents_key)
                        ]
                        + [
                            tex.escape(item)
                            for item in sorted(contents - segment, key=contents_key)
                        ]
                    )

            return rf"\node[{kind}] at ({position}) {{{label}}};"

        case Gain(gained=gained):
            gained = format_contents(gained)
            return rf"\node[gain={{black}}{{{gained}}}] at ({position}) {{}};"

        case Loss(contents=contents, segment=segment):
            if segment is not None and segment != contents:
                segment = format_contents(segment)
            else:
                segment = ""

            return rf"\node[loss={{black}}{{{segment}}}] at ({position}) {{}};"

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
            r"\usepackage{varwidth}"
            r"\usepackage{tikz}"
            r"\usetikzlibrary{patterns.meta}"
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
        # TODO: Connect to intermediate children for non-binary trees
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
            host_label = f"node[host label] {{{tex.escape(layout.host.name)}}}"
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

    \usepackage{varwidth}
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
