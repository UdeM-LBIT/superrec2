from sowing.node import Node
from collections.abc import Mapping
from superrec2.utils.algebras import (
    make_unit_magma,
    make_unit_generator,
    make_generator,
    min_plus,
    count,
)
from superrec2.model.history import Event, Codiverge, Diverge, Gain, Loss, Extant


def reconciliation_algorithm(algo):
    """
    Wrap a constitent interface around a reconciliation algorithm.

    The wrapped algorithm can accept any input, even non-binary ones, and
    automatically searches through all possible binarizations.

    :param algo: original algorithm
    :param progress: callback to report progress, if needed (default: no-op)
    :returns: wrapped algorithm
    """

    def reconcile(setting, structure, progress=lambda iterator, **kwargs: iterator):
        count_bins = sum(1 for _ in setting.binarize())
        return sum(
            (
                algo(bin_setting, structure)
                for bin_setting in progress(setting.binarize(), total=count_bins)
            ),
            start=structure.null(),
        ).value

    return reconcile


def make_cost_algebra(typename: str, costs: Mapping[str, float]):
    """Create a structure that computes reconciliation costs."""

    def make(event: Event):
        match event:
            case Extant():
                return 0

            case Gain():
                return 0

            case Codiverge():
                return costs.get("speciation", 0)

            case Diverge():
                if not event.cut and not event.transfer:
                    return costs.get("dup", 1)

                if not event.cut and event.transfer:
                    return costs.get("transfer-dup", 1)

                if event.cut and not event.transfer:
                    return costs.get("cut", 1)

                if event.cut and event.transfer:
                    return costs.get("transfer-cut", 1)

            case Loss():
                return costs.get("loss", 1)

        raise ValueError

    return min_plus(typename, make)


def join_event_nodes(node1: Node, node2: Node) -> Node:
    """Append an event node to another."""
    if node1.data is None:
        return node2

    if node2.data is None:
        return node1

    return node1.add(node2)


history_builder = make_unit_magma(
    "history_builder",
    unit=Node(),
    mul=join_event_nodes,
    make=Node,
)

history_counter = count("history_counter", lambda _: 1)
history_generator = make_generator("history_generator", history_builder)
history_unit_generator = make_unit_generator("history_unit_generator", history_builder)
