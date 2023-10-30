from sowing.node import Node
from sowing.zipper import Zipper
from sowing import traversal
from dataclasses import replace
from ...model.history import Associate, Event, Gain, Loss, Diverge, History


Contents = frozenset[str]
AssociateNode = Node[Associate, None]
EXTRA_CONTENTS = "__extra__"


def compute_min_contents(
    associate_tree: AssociateNode,
) -> dict[AssociateNode, Contents]:
    """
    Compute the minimum contents that must be assigned to each associate of a tree.

    :param associate_tree: input associate tree
    :returns: sets minimum contents for each node
    """
    # Compute the contents appearing in the subtree of each node in an associate tree
    min_contents: dict[AssociateNode, Contents] = {}

    for cursor in traversal.depth(associate_tree, preorder=False):
        if cursor.is_leaf():
            min_contents[cursor] = frozenset()
        else:
            left_contents = min_contents[cursor.down(0)]
            right_contents = min_contents[cursor.down(1)]
            min_contents[cursor] = left_contents | right_contents

        if cursor.node.data.contents:
            min_contents[cursor] |= cursor.node.data.contents

    # Compute the contents gained above each node of the associate tree,
    # so that each item is gained right above its lowest common ancestor
    gains: dict[AssociateNode, Contents] = {}
    gains[associate_tree.unzip()] = min_contents[associate_tree.unzip()]

    for cursor in traversal.depth(associate_tree, preorder=True):
        if not cursor.is_leaf():
            # Contents to be gained somewhere in the subtree of this node
            gained_below = gains[cursor]

            # Move gains downwards unless they are shared by both children
            left_child = cursor.down(0)
            left = min_contents[left_child]

            right_child = cursor.down(1)
            right = min_contents[right_child]

            gains[cursor] = gained_below & left & right
            gains[left_child] = gained_below & (left - right)
            gains[right_child] = gained_below & (right - left)

            # Take away contents gained in the children
            min_contents[cursor] -= gains[left_child] | gains[right_child]

    return min_contents


def propagate_contents(history: History) -> History:
    """Replace all extra-contents placeholders with actual contents."""

    def propagate_at_node(cursor: Zipper[Event, None]) -> Zipper[Event, None]:
        parent_event = cursor.node.data
        edges = cursor.node.edges
        new_edges = []

        for edge in edges:
            event = edge.node.data
            extra = parent_event.contents - event.contents

            if EXTRA_CONTENTS in event.contents:
                event = replace(
                    event, contents=event.contents - {EXTRA_CONTENTS} | extra
                )

            if isinstance(event, (Loss, Diverge)) and EXTRA_CONTENTS in event.segment:
                event = replace(event, segment=event.segment - {EXTRA_CONTENTS} | extra)

            if isinstance(event, Gain) and EXTRA_CONTENTS in event.gained:
                event = replace(event, segment=event.segment - {EXTRA_CONTENTS} | extra)

            new_edges.append(edge.replace(node=edge.node.replace(data=event)))

        return cursor.replace(node=cursor.node.replace(edges=tuple(new_edges)))

    return History(
        host_tree=history.host_tree,
        event_tree=traversal.fold(
            propagate_at_node,
            traversal.depth(history.event_tree, preorder=True),
        ),
    )
