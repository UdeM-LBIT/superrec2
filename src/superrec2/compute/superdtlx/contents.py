from sowing.node import Node
from sowing import traversal
from ...model.history import Associate


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
