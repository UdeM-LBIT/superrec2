from sowing.node import Node
from typing import NamedTuple, Self
from superrec2.utils.algebras import (
    tuple_ordered_magma,
    pareto,
    make_unit_magma,
    make_unit_generator,
    make_generator,
    min_plus,
    count,
)
from superrec2.model.history import Event, Codiverge, Diverge, Gain, Loss, Extant


class DummyProgress:
    def __init__(self, *args, **kwargs):
        pass

    def update(self, count=1):
        pass


def reconciliation_algorithm(algo):
    """
    Wrap a consistent interface around a reconciliation algorithm.

    The wrapped algorithm can accept any input, even non-binary ones, and
    automatically searches through all possible binarizations.

    :param algo: original algorithm
    :param progress: callback to report progress, if needed (default: no-op)
    :returns: wrapped algorithm
    """

    def reconcile(
        setting,
        structure,
        progress=DummyProgress,
    ):
        bar = progress(total=sum(1 for _ in setting.binarize()))
        result = structure.null()

        for item in map(
            lambda binary_setting: algo(binary_setting, structure),
            setting.binarize(),
        ):
            result += item
            bar.update()

        return result.value

    return reconcile


class EventCosts(NamedTuple):
    speciation: float = 0
    duplication: float = 1
    transfer_duplication: float = 1
    cut: float = 1
    transfer_cut: float = 1
    loss: float = 1


def make_cost_algebra(typename: str, costs: EventCosts):
    """Create a structure that computes reconciliation costs."""

    def make(event: Event):
        match event:
            case Extant():
                return 0

            case Gain():
                return 0

            case Codiverge():
                return costs.speciation

            case Diverge():
                if not event.cut and not event.transfer:
                    return costs.duplication

                if not event.cut and event.transfer:
                    return costs.transfer_duplication

                if event.cut and not event.transfer:
                    return costs.cut

                if event.cut and event.transfer:
                    return costs.transfer_cut

            case Loss():
                return costs.loss

        raise ValueError

    return min_plus(typename, make)


@tuple_ordered_magma
class EventVector(NamedTuple):
    speciation: int = 0
    duplication: int = 0
    transfer_duplication: int = 0
    cut: int = 0
    transfer_cut: int = 0
    loss: int = 0

    @classmethod
    def make(cls, event: Event) -> Self:
        match event:
            case Extant() | Gain():
                return EventVector()

            case Codiverge():
                return EventVector(speciation=1)

            case Diverge():
                if not event.cut and not event.transfer:
                    return EventVector(duplication=1)

                if not event.cut and event.transfer:
                    return EventVector(transfer_duplication=1)

                if event.cut and not event.transfer:
                    return EventVector(cut=1)

                if event.cut and event.transfer:
                    return EventVector(transfer_cut=1)

            case Loss():
                return EventVector(loss=1)

            case _:
                raise ValueError(f"unknown event type {type(event)}")


event_vector_pareto = pareto("event_vector_pareto", EventVector)


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
