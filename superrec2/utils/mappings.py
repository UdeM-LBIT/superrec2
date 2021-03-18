"""Helpers for transforming discrete mappings."""
from collections import defaultdict
from typing import DefaultDict, List, Mapping, TypeVar


From = TypeVar("From")
To = TypeVar("To")


def invert_mapping(mapping: Mapping[From, To]) -> DefaultDict[To, List[From]]:
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
