"""Handle and compare sequences and subsequences."""
from typing import List, Sequence, TypeVar


Element = TypeVar("Element")


def subseq_complete(sequence: Sequence[Element]) -> int:
    """Get a bitmask representing a complete subsequence."""
    return (1 << len(sequence)) - 1


def mask_from_subseq(
    child: Sequence[Element],
    parent: Sequence[Element],
) -> int:
    """
    Create a bitmask representing a subsequence.

    Note that :param:`child` must be a subsequence of :param:`parent`, i.e.,
    it can only contain elements from :param:`parent`, appearing in the same
    relative order as in :param:`parent`.

    :param child: subsequence
    :param parent: parent sequence
    :returns: bitmask representing the :param:`child` subsequence
    """
    child_i = 0
    mask = 0

    for parent_i, parent_v in enumerate(parent):
        if child_i == len(child):
            break

        if child[child_i] == parent_v:
            mask |= 1 << parent_i
            child_i += 1

    return mask


def subseq_from_mask(
    child: int,
    parent: Sequence[Element],
) -> List[Element]:
    """
    Reconstruct a subsequence from a bitmask.

    Note that :param:`child` must not have more set bits than the number of
    elements in :param:`parent`.

    :param child: subsequence bitmask
    :param parent: parent sequence
    :returns: reconstructed subsequence
    """
    result = []
    parent_i = 0

    while child:
        if child & 1:
            result.append(parent[parent_i])

        child >>= 1
        parent_i += 1

    return result


def subseq_segment_dist(
    child: int,
    parent: int,
    edges: bool,
) -> int:
    """
    Count the number of lost segments between two subsequences.

    :param child: subsequence bitmask
    :param parent: parent subsequence bitmask
    :param edges: if True, count lost sequences on the edges of
        :param:`parent`, if False, ignore them
    :returns: number of lost segments from :param:`parent` to :param:`child`,
        or -1 if :param:`child` is not a subsequence of :param:`parent`
    """
    in_segm = not edges
    dist = 0

    if parent.bit_length() < child.bit_length():
        return -1

    for _ in range(parent.bit_length()):
        bit_child = child & 1
        bit_parent = parent & 1

        if bit_child and not bit_parent:
            return -1

        if bit_parent:
            if not bit_child:
                if not in_segm:
                    dist += 1
                    in_segm = True
            elif in_segm:
                in_segm = False

        child >>= 1
        parent >>= 1

    if in_segm and not edges:
        dist -= 1

    return dist
