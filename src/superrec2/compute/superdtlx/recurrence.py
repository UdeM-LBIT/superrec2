from sowing import traversal
from sowing.indexed import IndexedTree
from itertools import product
from enum import Enum, auto
from collections import defaultdict
from typing import TypeVar
from ...utils.algebras import Semiring
from ...model.history import (
    Host,
    Reconciliation,
    Extant,
    Codiverge,
    Diverge,
)
from ..util import reconciliation_algorithm
from .paths import make_path
from .contents import AssociateNode, Contents, EXTRA_CONTENTS, compute_min_contents


T = TypeVar("T")


class HostChoice(Enum):
    # Starts at the specified parent host, ends at any descendant
    Incoming = auto()

    # Starts at the left child of the parent host, ends at any of its descendants
    Left = auto()

    # Starts at the right child of the parent host, ends at any of its descendants
    Right = auto()

    # Starts and ends at any separate host
    Separate = auto()


class ContentsChoice(Enum):
    # Starts with the equivalent contents as the parent contents
    Incoming = auto()

    # Starts with the minimal contents subset
    Minimal = auto()


def compute_choices_at(
    node: AssociateNode,
    incoming_host: str,
    incoming_contents: Contents,
    min_contents: Contents,
    host_index: IndexedTree[Host, None],
    structure: type[Semiring[T]],
    table: dict[tuple[AssociateNode, str, Contents], Semiring[T]],
) -> dict[tuple[HostChoice, ContentsChoice], Semiring[T]]:
    choices = defaultdict(structure.null)
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
                path=table[(node, end_host, end_contents)],
            )

    return choices


def join_binary_event(
    host: str,
    contents: Contents,
    left_contents: Contents,
    right_contents: Contents,
    structure: type[Semiring[T]],
    left_choices: dict[tuple[HostChoice, ContentsChoice], Semiring[T]],
    right_choices: dict[tuple[HostChoice, ContentsChoice], Semiring[T]],
) -> Semiring[T]:
    # Speciation with matching host-children order
    results = (
        structure.make(Codiverge(host=host, contents=contents))
        * left_choices[(HostChoice.Left, ContentsChoice.Incoming)]
        * right_choices[(HostChoice.Right, ContentsChoice.Incoming)]
    )

    # Speciation with reverse host-children order
    results += (
        structure.make(Codiverge(host=host, contents=contents))
        * left_choices[(HostChoice.Right, ContentsChoice.Incoming)]
        * right_choices[(HostChoice.Left, ContentsChoice.Incoming)]
    )

    # Duplication
    if right_contents == contents:
        results += (
            structure.make(
                Diverge(
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
            structure.make(
                Diverge(
                    host=host,
                    contents=contents,
                    segment=left_contents,
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
        structure.make(
            Diverge(
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
        structure.make(
            Diverge(
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
            structure.make(
                Diverge(
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
            structure.make(
                Diverge(
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
            structure.make(
                Diverge(
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


@reconciliation_algorithm
def reconcile(setting: Reconciliation, structure: type[Semiring[T]]) -> Semiring[T]:
    results = defaultdict(structure.null)
    root = setting.associate_tree
    min_contents = compute_min_contents(root)

    for cursor in traversal.depth(root, preorder=False):
        node = cursor.node

        if cursor.is_leaf():
            name = node.data.name
            host = node.data.host
            contents = node.data.contents
            value = structure.make(Extant(name=name, host=host, contents=contents))
            results[(node, host, contents)] += value
        else:
            for host, contents in product(
                setting.host_index.keys(),
                (min_contents[cursor], min_contents[cursor] | {EXTRA_CONTENTS}),
            ):
                left = cursor.down(0)
                right = cursor.down(1)

                choices_args = {
                    "incoming_host": host,
                    "incoming_contents": contents,
                    "host_index": setting.host_index,
                    "structure": structure,
                    "table": results,
                }
                left_choices = compute_choices_at(
                    node=left.node, min_contents=min_contents[left], **choices_args
                )
                right_choices = compute_choices_at(
                    node=right.node, min_contents=min_contents[right], **choices_args
                )

                results[(node, host, contents)] += join_binary_event(
                    host=host,
                    contents=contents,
                    left_contents=min_contents[left] & contents,
                    right_contents=min_contents[right] & contents,
                    structure=structure,
                    left_choices=left_choices,
                    right_choices=right_choices,
                )

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
                path=results[(root, host, root_contents)],
            )
            for host in setting.host_index.keys()
        ),
        start=structure.null(),
    )
