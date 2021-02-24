from collections import defaultdict
from typing import DefaultDict, List, Mapping, TypeVar


A = TypeVar('A')
B = TypeVar('B')


def invert_mapping(mapping: Mapping[A, B]) -> DefaultDict[B, List[A]]:
    """
    Compute the reciprocal of a mapping.

    :param mapping: initial mapping to invert
    :returns: for each value of the input mapping, a list of
        keys for which that value is used
    """
    result = defaultdict(list)

    for key, value in mapping.items():
        result[value].append(key)

    return result
