"""Compute layouts for reconciliations."""

from sowing import traversal
from sowing.node import Node
from math import inf
from .model import EventLayout, DrawParams, HostLayout, Layout, Orientation
from ..model.history import Event, History
from ..utils.geometry import Position, Rect, Size


def _init_layout(
    history: History,
    rects: dict[Event, Size],
    params: DrawParams,
) -> Layout:
    """
    Initialize the structures used for computing the layout of each host and event.

    :param history: history to layout
    :param rects: sizes for each event node of the history
    :param params: layout parameters
    :returns: initial layout information for each host node
    """
    # Initialize host layout information
    layout = {}

    for cursor in traversal.depth(history.host_tree, preorder=False):
        node = cursor.node
        host = node.data
        layout[host.name] = HostLayout(
            params=params,
            host=host,
            children={
                edge.node.data.name: layout[edge.node.data.name] for edge in node.edges
            },
        )

    # Initialize event layout information
    for cursor in traversal.depth(history.event_tree):
        event = cursor.node.data
        descendants = [
            host.node.data.name for host in history.host_index[event.host].node.edges
        ]

        in_children = []
        desc_children = []
        side_children = []

        for edge in cursor.node.edges:
            child = edge.node
            host = child.data.host

            if host == event.host:
                in_children.append(child)
            elif host in descendants:
                desc_children.append(child)
            else:
                side_children.append(child)

        rect = rects[event]

        if params.orientation == Orientation.Horizontal:
            rect = Rect(
                Position(rect.position.y, rect.position.x),
                Size(rect.size.h, rect.size.w),
            )

        host_layout = layout[event.host]
        host_layout.events[cursor.node] = EventLayout(
            in_children=in_children,
            desc_children=desc_children,
            side_children=side_children,
            area=rect - rect.top_left(),
            anchor=-rect.top_left(),
        )

    # Add dummy events inside empty hosts
    for host_layout in layout.values():
        if not host_layout.events:
            host_layout.events[Node(object())] = EventLayout(
                in_children=[],
                desc_children=[],
                side_children=[],
            )

    return layout


def _layout_fork(host_layout: HostLayout, epoch_height: int) -> None:
    """
    Position the events inside the fork of an host.

    :param host_layout: host layout information to update
    :param epoch_height: total height on the main axis of the epochs below
    """
    params = host_layout.params

    # Position forking and leaf events inside the fork
    next_pos_main = -epoch_height - params.events_host_padding
    next_pos_cross = params.events_host_padding

    for node, event_layout in host_layout.events.items():
        size = event_layout.area.size
        event_layout += Position(next_pos_cross, next_pos_main - size.h)

        if event_layout.forking:
            next_pos_cross += size.w / 2 + params.events_spacing
            next_pos_main -= size.h + params.events_spacing
        elif event_layout.leaf:
            next_pos_cross += params.events_spacing + size.w


def _layout_inner(layout: Layout, event: Event) -> None:
    """
    Position an event inside the trunk of an host.

    :param layout: layout of all hosts
    :param event_layout: event layout information to update
    """
    host_layout = layout[event.data.host]
    event_layout = host_layout.events[event]
    params = host_layout.params

    if event_layout.forking or event_layout.leaf:
        return

    event_layout -= event_layout.area.position
    anchor = event_layout.anchor
    size = event_layout.area.size

    # Center event above the anchors of its inner children
    cross_area = Rect.fit(
        layout[child.data.host].events[child].anchor
        for child in event_layout.in_children
    )
    cross_offset = cross_area.center().x

    # Position event above the anchors of its inner and outside children
    # and above the current forking region
    main_area = Rect.fit(
        layout[child.data.host].events[child].area.top()
        for child in event_layout.in_children + event_layout.side_children
    )
    main_offset = min(
        main_area.top().y - params.events_spacing, host_layout.fork_area.top().y
    )

    event_layout += Position(cross_offset - anchor.x, main_offset - size.h)

    # Add dummy events at the transfer location of outside children
    for child in event_layout.side_children:
        child_anchor = layout[child.data.host].events[child].anchor
        position = event_layout.anchor.meet_hv(child_anchor)
        layout[child.data.host].events[Node(object())] = EventLayout(
            in_children=[child],
            desc_children=[],
            side_children=[],
            area=Rect.fit((position,)),
            anchor=position,
        )


def _layout_children(host_layout: HostLayout) -> None:
    """
    Position the children of an host along the cross axis.

    :param host_layout: host layout information to update
    """
    params = host_layout.params
    children = list(host_layout.children.values())

    if not children:
        return

    # Reset current host cross axis position to zero
    initial_offset = host_layout.trunk_area.left().x

    for node, event_layout in host_layout.events.items():
        event_layout += Position(-initial_offset, 0)

    # Position children hosts along the cross axis
    next_pos_cross = 0

    for child_layout in children:
        child_size = child_layout.area.size
        child_layout += Position(-child_layout.area.left().x + next_pos_cross, 0)
        next_pos_cross += child_size.w + params.subtree_spacing

    # Shift last child to leave enough room for the current fork
    leftmost_subfork = max(event.anchor.x for event in children[0].events.values())
    rightmost_subfork = min(event.anchor.x for event in children[-1].events.values())
    children_span = rightmost_subfork - leftmost_subfork
    missing_space = host_layout.trunk_area.size.w - children_span

    if missing_space > 0:
        children[-1] += Position(missing_space, 0)
        rightmost_subfork += missing_space

    # Move fork events at the center above the children hosts
    children_center = (leftmost_subfork + rightmost_subfork) / 2
    offset = children_center - host_layout.trunk_area.size.w / 2

    for node, event_layout in host_layout.events.items():
        event_layout += Position(offset, 0)


def compute(
    history: History,
    rects: dict[Event, Size],
    params: DrawParams = DrawParams(),
) -> Layout:
    """
    Compute a layout for an evolutionary history.

    :param history: history to layout
    :param rects: sizes for each event node of the history
    :param params: layout parameters
    :returns: layout information for each host node
    """
    epochs = history.epochs()
    layout = _init_layout(history, rects, params)
    epoch_height = 0

    for epoch in reversed(epochs.range()):
        # Position the forks of each ending host along the main axis
        for cursor in epochs.hosts_at(end=epoch):
            _layout_fork(layout[cursor.node.data.name], epoch_height)

        # Layout internal events of the current epoch in post-order
        for cursor in traversal.depth(history.event_tree, preorder=False):
            if epochs.events[cursor] == epoch:
                _layout_inner(layout, cursor.node)

        # Position the children of each starting host along the cross axis
        for cursor in epochs.hosts_at(start=epoch):
            _layout_children(layout[cursor.node.data.name])

        epoch_height = max(
            abs(layout[cursor.node.data.name].area.top().y)
            for cursor in epochs.hosts_at(start=epoch)
        )

    return layout
