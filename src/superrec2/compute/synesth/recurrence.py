from sowing import traversal
from sowing.indexed import IndexedTree
from itertools import product
from enum import Enum, auto
from collections import defaultdict
from typing import TypeVar
from ...utils.algebras import SemiRing, Structure
from ...model.history import (
    Host,
    Reconciliation,
    Event,
    Extant,
    Codiverge,
    Diverge,
)
from .paths import make_path
from .contents import AssociateNode, Contents, EXTRA_CONTENTS, compute_min_contents


T = TypeVar("T")


class HostChoice(Enum):
    # Starts at the specified parent host, ends at any non-ancestor
    Incoming = auto()

    # Starts at the left child of the parent host, ends at any non-ancestor
    Left = auto()

    # Starts at the right child of the parent host, ends at any non-ancestor
    Right = auto()

    # Starts and ends at any separate host
    Separate = auto()


class ContentsChoice(Enum):
    # Starts with the equivalent contents as the parent contents
    Incoming = auto()

    # Starts with the minimal contents subset
    Minimal = auto()


def compute_choices_at(
    incoming_host: str,
    incoming_contents: Contents,
    min_contents: Contents,
    host_index: IndexedTree[Host, None],
    structure: Structure[T, [Event]],
    table: dict[tuple[str, Contents], SemiRing[T]],
) -> dict[tuple[HostChoice, ContentsChoice], SemiRing[T]]:
    choices = defaultdict(lambda: structure.zero)
    try_start_hosts = []

    host_cursor = host_index[incoming_host]
    left_host = None if host_cursor.is_leaf() else host_cursor.down(0).node.data.name
    right_host = None if host_cursor.is_leaf() else host_cursor.down(1).node.data.name

    for item in host_index.keys():
        if item == incoming_host:
            try_start_hosts.append((HostChoice.Incoming, item))
        elif item == left_host:
            try_start_hosts.append((HostChoice.Left, item))
        elif item == right_host:
            try_start_hosts.append((HostChoice.Right, item))
        elif not host_index.is_comparable(item, incoming_host):
            try_start_hosts.append((HostChoice.Separate, item))

    try_start_contents = [(ContentsChoice.Minimal, min_contents & incoming_contents)]

    if EXTRA_CONTENTS in incoming_contents or not (incoming_contents <= min_contents):
        try_start_contents.append(
            (
                ContentsChoice.Incoming,
                (min_contents & incoming_contents) | {EXTRA_CONTENTS},
            )
        )

    if EXTRA_CONTENTS not in incoming_contents:
        try_start_contents.append((ContentsChoice.Incoming, incoming_contents))

    for (
        (host_choice, start_host),
        (contents_choice, start_contents),
    ) in product(try_start_hosts, try_start_contents):
        if host_choice == HostChoice.Separate:
            try_end_hosts = (start_host,)
        else:
            try_end_hosts = (
                item
                for item in host_index.keys()
                if not host_index.is_strict_ancestor_of(item, start_host)
            )

        try_end_contents = (min_contents, min_contents | {EXTRA_CONTENTS})

        for end_host, end_contents in product(try_end_hosts, try_end_contents):
            choices[(host_choice, contents_choice)] += make_path(
                start_host=start_host,
                start_contents=start_contents,
                end_host=end_host,
                end_contents=end_contents,
                host_index=host_index,
                structure=structure,
                path=table[(end_host, end_contents)],
            )

    return choices


def join_binary_event(
    host: str,
    contents: Contents,
    left_contents: Contents,
    right_contents: Contents,
    structure: Structure[T, [Event]],
    left_choices: dict[tuple[HostChoice, ContentsChoice], SemiRing[T]],
    right_choices: dict[tuple[HostChoice, ContentsChoice], SemiRing[T]],
) -> SemiRing[T]:
    # Speciation with matching host-children order
    results = (
        structure(Codiverge(apparent=True, host=host, contents=contents))
        * left_choices[(HostChoice.Left, ContentsChoice.Incoming)]
        * right_choices[(HostChoice.Right, ContentsChoice.Incoming)]
    )

    # Speciation with reverse host-children order
    results += (
        structure(Codiverge(apparent=True, host=host, contents=contents))
        * left_choices[(HostChoice.Right, ContentsChoice.Incoming)]
        * right_choices[(HostChoice.Left, ContentsChoice.Incoming)]
    )

    # Duplication
    if right_contents == contents:
        results += (
            structure(
                Diverge(
                    apparent=True,
                    host=host,
                    contents=contents,
                    segment=left_contents,
                    cut=False,
                    transfer=False,
                    result=0,
                )
            )
            * left_choices[(HostChoice.Incoming, ContentsChoice.Minimal)]
            * right_choices[(HostChoice.Incoming, ContentsChoice.Incoming)]
        )
    else:
        results += (
            structure(
                Diverge(
                    apparent=True,
                    host=host,
                    contents=contents,
                    segment=right_contents,
                    cut=False,
                    transfer=False,
                    result=1,
                )
            )
            * left_choices[(HostChoice.Incoming, ContentsChoice.Incoming)]
            * right_choices[(HostChoice.Incoming, ContentsChoice.Minimal)]
        )

    # Duplication-transfer to the left
    results += (
        structure(
            Diverge(
                apparent=True,
                host=host,
                contents=contents,
                segment=left_contents,
                cut=False,
                transfer=True,
                result=0,
            )
        )
        * left_choices[(HostChoice.Separate, ContentsChoice.Minimal)]
        * right_choices[(HostChoice.Incoming, ContentsChoice.Incoming)]
    )

    # Duplication-transfer to the right
    results += (
        structure(
            Diverge(
                apparent=True,
                host=host,
                contents=contents,
                segment=right_contents,
                cut=False,
                transfer=True,
                result=1,
            )
        )
        * left_choices[(HostChoice.Incoming, ContentsChoice.Incoming)]
        * right_choices[(HostChoice.Separate, ContentsChoice.Minimal)]
    )

    if (
        left_contents | right_contents == contents
        and not left_contents & right_contents
    ):
        # Cut (symmetric)
        results += (
            structure(
                Diverge(
                    apparent=True,
                    host=host,
                    contents=contents,
                    segment=left_contents,
                    cut=True,
                    transfer=False,
                    result=0,
                )
            )
            * left_choices[(HostChoice.Incoming, ContentsChoice.Minimal)]
            * right_choices[(HostChoice.Incoming, ContentsChoice.Minimal)]
        )

        # Cut-transfer to the left
        results += (
            structure(
                Diverge(
                    apparent=True,
                    host=host,
                    contents=contents,
                    segment=left_contents,
                    cut=True,
                    transfer=True,
                    result=0,
                )
            )
            * left_choices[(HostChoice.Separate, ContentsChoice.Minimal)]
            * right_choices[(HostChoice.Incoming, ContentsChoice.Minimal)]
        )

        # Cut-transfer to the right
        results += (
            structure(
                Diverge(
                    apparent=True,
                    host=host,
                    contents=contents,
                    segment=right_contents,
                    cut=True,
                    transfer=True,
                    result=1,
                )
            )
            * left_choices[(HostChoice.Incoming, ContentsChoice.Minimal)]
            * right_choices[(HostChoice.Separate, ContentsChoice.Minimal)]
        )

    return results


def solve_binary(
    setting: Reconciliation, structure: Structure[T, [Event]]
) -> SemiRing[T]:
    """
    Solve a reconciliation problem over a given algebraic structure.

    :param setting: setting containing a host tree and an associate tree to reconcile;
        both trees must be binary trees
    :param structure: semiring structure with a morphism from events to the semiring
    :return: resulting semiring value
    """
    results: dict[AssociateNode, dict[tuple[str, Contents], SemiRing[T]]] = {}
    root = setting.associate_tree
    min_contents = compute_min_contents(root)

    for cursor in traversal.depth(root, preorder=False):
        node = cursor.node
        table = defaultdict(lambda: structure.zero)

        if cursor.is_leaf():
            name = node.data.name
            host = node.data.host
            contents = node.data.contents
            value = structure(
                Extant(name=name, host=host, contents=contents, apparent=True)
            )
            table[(host, contents)] += value
        else:
            left = cursor.down(0)
            right = cursor.down(1)

            for host, contents in product(
                setting.host_index.keys(),
                (min_contents[cursor], min_contents[cursor] | {EXTRA_CONTENTS}),
            ):
                choices_args = {
                    "incoming_host": host,
                    "incoming_contents": contents,
                    "host_index": setting.host_index,
                    "structure": structure,
                }
                left_choices = compute_choices_at(
                    min_contents=min_contents[left],
                    table=results[left.node],
                    **choices_args,
                )
                right_choices = compute_choices_at(
                    min_contents=min_contents[right],
                    table=results[right.node],
                    **choices_args,
                )

                table[(host, contents)] += join_binary_event(
                    host=host,
                    contents=contents,
                    left_contents=min_contents[left] & contents,
                    right_contents=min_contents[right] & contents,
                    structure=structure,
                    left_choices=left_choices,
                    right_choices=right_choices,
                )

            del results[left.node]
            del results[right.node]

        results[node] = table

    root_contents = min_contents[root.unzip()]
    return sum(
        (
            make_path(
                start_host=host,
                end_host=host,
                start_contents=frozenset(),
                end_contents=root_contents,
                host_index=setting.host_index,
                structure=structure,
                path=results[root][(host, root_contents)],
            )
            for host in setting.host_index.keys()
        ),
        start=structure.zero,
    )
