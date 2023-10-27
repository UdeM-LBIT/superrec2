"""Compute layouts for reconciliations."""
from typing import Callable
from itertools import chain
from sowing import traversal
from .model import EventLayout, DrawParams, HostLayout, Layout, Orientation
from ..model.history import Event, History
from ..utils.geometry import Position, Rect, Size


def _process_events(
    history: History,
    measure_events: Callable[[list[Event], DrawParams], dict[Event, Size]],
    params: DrawParams,
) -> Layout:
    """
    Measure event nodes and group them by their host.

    :param history: history containing the event tree to process
    :param measure_events: callback to compute the sizes of a set of events
    :param params: drawing settings
    :returns: resulting layout
    """
    layout = {
        host_name: HostLayout(
            host=host_cursor.node.data,
            children=[edge.node.data.name for edge in host_cursor.node.edges],
        )
        for host_name, host_cursor in history.host_index.items()
    }

    sizes = measure_events(
        events=(cursor.node.data for cursor in traversal.depth(history.event_tree)),
        params=params,
    )

    for cursor in traversal.depth(history.event_tree):
        event = cursor.node.data
        in_children = []
        out_children = []

        for edge in cursor.node.edges:
            child = edge.node

            if child.data.host == event.host:
                in_children.append(child)
            else:
                out_children.append(child)

        size = sizes[event]

        if params.orientation == Orientation.Horizontal:
            size = Size(size.h, size.w)

        host_layout = layout[event.host]
        host_layout.events[cursor.node] = EventLayout(
            forking=bool(out_children) and not in_children,
            in_children=in_children,
            out_children=out_children,
            area=Rect(Position(0, 0), size),
        )

    return layout


def _layout_events(layout: Layout, params: DrawParams) -> Layout:
    """
    Compute the relative position of each event inside their host.

    :param layout: layout including the size of each event
    :param params: drawing settings
    :returns: layout updated with event positions
    """
    for host, host_layout in layout.items():
        next_pos_main = 0
        next_pos_cross = 0

        for node, event_layout in host_layout.events.items():
            size = event_layout.area.size

            if event_layout.forking:
                # Align forking nodes to the same diagonal
                next_pos_main -= size.w
                pos = Position(next_pos_main, next_pos_cross)
                next_pos_main -= params.gene_branch_spacing
                next_pos_cross += size.h + params.gene_branch_spacing
            elif event_layout.in_children:
                # Align internal nodes above the center of their children
                children_area = Rect.fit(
                    host_layout.events[child].area for child in event_layout.in_children
                )
                pos = children_area.top() + Position(
                    -size.w / 2,
                    -params.species_branch_padding - size.h,
                )
            else:
                # Align leaves across the same axis
                next_pos_main -= size.w
                pos = Position(next_pos_main, -size.h)
                next_pos_main -= params.events_spacing

            event_layout.area = Rect(pos, size)

        # Compute overall events area
        host_layout.events_area = Rect.fit(
            chain(
                (Rect.zero(),),
                (event.area for event in host_layout.events.values()),
            )
        ).grow(params.events_host_padding)

        # Make the upper-left corner be the origin
        delta = host_layout.events_area.position
        host_layout.events_area -= delta

        for event_layout in host_layout.events.values():
            event_layout.area -= delta

    return layout


def _layout_hosts(
    layout: Layout,
    epochs: dict[str, int],
    params: DrawParams,
):
    """
    Compute the overall size and position of each host subtree relative to its parent.

    :param layout: layout including the position and size of each event
    :param epochs: hosts listed by their epoch; all hosts of the same epoch will
        be rendered at the same level in time
    :param params: drawing settings
    :returns: layout updated with host subtree sizes and positions
    """
    epochs_inv = {}
    epochs_heights = {}

    for host, epoch in epochs.items():
        epochs_inv.setdefault(epoch, []).append(host)

    for epoch, hosts in sorted(epochs_inv.items(), reverse=True):
        # Compute minimum height to fit events across all hosts of the epoch
        epoch_height = max(layout[host].events_area.size.h for host in hosts)
        epoch_height += params.epoch_spacing
        epochs_heights[epoch] = epoch_height

        for host in hosts:
            host_layout = layout[host]

            # Grow events area to reach epoch height, if needed
            grow = epoch_height - params.epoch_spacing - host_layout.events_area.size.h
            host_layout.events_area = Rect(
                position=host_layout.events_area.position,
                size=Size(
                    host_layout.events_area.size.w,
                    host_layout.events_area.size.h + grow,
                ),
            )

            for event_layout in host_layout.events.values():
                event_layout.area += Position(0, grow)

            # Balance children hosts around the left and right of the fork
            children_widths = []
            total_children_width = 0

            for child_host in host_layout.children:
                child_layout = layout[child_host]
                children_widths.append(child_layout.area.size.w)
                total_children_width += child_layout.area.size.w

            running_children_width = 0
            split_at = 0

            for split_at in range(len(host_layout.children)):
                if running_children_width * 2 >= total_children_width:
                    break

                running_children_width += children_widths[split_at]

            # Position left-side children
            child_position = 0

            for child_host in host_layout.children[:split_at]:
                child_layout = layout[child_host]
                child_size = child_layout.area.size
                child_layout.area += Position(child_position - child_size.w, 0)

                for child_epoch in range(epoch, epochs[child_host]):
                    child_layout.area += Position(0, epochs_heights[child_epoch])

                child_position -= child_size.w + params.min_subtree_spacing

            # Position right-side children
            child_position = host_layout.events_area.size.w

            for child_host in host_layout.children[split_at:]:
                child_layout = layout[child_host]
                child_size = child_layout.area.size
                child_layout.area += Position(child_position, 0)

                for child_epoch in range(epoch, epochs[child_host]):
                    child_layout.area += Position(0, epochs_heights[child_epoch])

                child_position += child_size.w + params.min_subtree_spacing

            # Compute overall area
            host_layout.fork_events_area = Rect.fit(
                chain(
                    (
                        Rect(
                            position=Position(
                                params.events_host_padding,
                                params.events_host_padding,
                            ),
                            size=Size.zero(),
                        ),
                    ),
                    (
                        event.area
                        for event in host_layout.events.values()
                        if event.forking
                    ),
                )
            ).grow(params.events_host_padding)

            host_layout.area = Rect.fit(
                chain(
                    (host_layout.events_area,),
                    (layout[child_host].area for child_host in host_layout.children),
                )
            )

            # Recompute coordinates to be relative to the overall upper-left corner
            delta_vec = host_layout.area.position
            host_layout.area -= delta_vec
            host_layout.events_area -= delta_vec
            host_layout.fork_events_area -= delta_vec

            for event in host_layout.events.values():
                event.area -= delta_vec

            for child_host in host_layout.children:
                layout[child_host].area -= delta_vec

    return layout


def _layout_absolute(layout: Layout, epochs: dict[int, list[str]]) -> Layout:
    """
    Make all positions in the layout absolute.

    :param layout: layout with relative positions
    :param epochs: hosts listed by epoch; all hosts in the same epoch are
        to be rendered at the same level in time
    :returns: layout updated with absolute positions
    """
    for epoch, hosts in sorted(epochs.items()):
        for host in hosts:
            host_layout = layout[host]
            delta = host_layout.area.position

            # Make event positions absolute
            host_layout.events_area += delta
            host_layout.fork_events_area += delta

            for event_layout in host_layout.events.values():
                event_layout.area += delta

            # Make child host positions absolute
            for child_host in host_layout.children:
                layout[child_host].area += delta

    return layout


def compute(
    history: History,
    measure_events: Callable[[list[Event], DrawParams], dict[Event, Size]],
    params: DrawParams = DrawParams(),
) -> Layout:
    """
    Compute a layout for an evolutionary history.

    :param history: history to layout
    :param measure_events: callback to compute the size of a set of events
    :param params: layout parameters
    :returns: layout information for each host node
    """
    layout = _process_events(history, measure_events, params)
    layout = _layout_events(layout, params)

    epochs = history.epochs()
    epochs_inv = {}

    for host, epoch in epochs.items():
        epochs_inv.setdefault(epoch, []).append(host)

    layout = _layout_hosts(layout, epochs, params)
    layout = _layout_absolute(layout, epochs_inv)
    return layout
