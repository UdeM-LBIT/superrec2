from sowing.node import Node
from typing import NamedTuple, Self
from dataclasses import dataclass
from superrec2.utils.algebras import (
    vector,
    Structure,
    Counter,
    pareto_of,
    generator_of,
    projector_of,
)
from superrec2.model.history import Event, Codiverge, Diverge, Gain, Loss, Extant


class DummyProgress:
    def __init__(self, *args, **kwargs):
        pass

    def update(self, count=1):
        pass


class DummyPool:
    def __init__(self, *args, **kwargs):
        pass

    def map(self, func, iterable):
        return map(func, iterable)


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
        pool=DummyPool(),
    ):
        bar = progress(total=sum(1 for _ in setting.binarize()))
        result = structure.zero

        for item in pool.map(
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

    def event_cost_morphism(self, event: Event):
        match event:
            case Extant():
                return 0

            case Gain():
                return 0

            case Codiverge():
                return self.speciation

            case Diverge():
                if not event.cut and not event.transfer:
                    return self.duplication

                if not event.cut and event.transfer:
                    return self.transfer_duplication

                if event.cut and not event.transfer:
                    return self.cut

                if event.cut and event.transfer:
                    return self.transfer_cut

            case Loss():
                return self.loss

            case _:
                raise ValueError(f"unknown event type {type(event)}")


@vector
class EventVector(NamedTuple):
    speciation: int = 0
    duplication: int = 0
    transfer_duplication: int = 0
    cut: int = 0
    transfer_cut: int = 0
    loss: int = 0


def event_vector_morphism(event: Event):
    match event:
        case Extant() | Gain():
            return frozenset({EventVector()})

        case Codiverge():
            return frozenset({EventVector(speciation=1)})

        case Diverge():
            if not event.cut and not event.transfer:
                return frozenset({EventVector(duplication=1)})

            if not event.cut and event.transfer:
                return frozenset({EventVector(transfer_duplication=1)})

            if event.cut and not event.transfer:
                return frozenset({EventVector(cut=1)})

            if event.cut and event.transfer:
                return frozenset({EventVector(transfer_cut=1)})

        case Loss():
            return frozenset({EventVector(loss=1)})

        case _:
            raise ValueError(f"unknown event type {type(event)}")


EventVectorPareto = pareto_of(EventVector)
event_vector_pareto = Structure(EventVectorPareto, event_vector_morphism)


@dataclass(frozen=True, slots=True)
class HistoryBuilder:
    value: Node

    def __mul__(node1, node2):
        if node1.value.data is None:
            return node2

        if node2.value.data is None:
            return node1

        return HistoryBuilder(node1.value.add(node2.value))


HistoryGenerator = generator_of(HistoryBuilder(Node()))
history_generator = Structure(
    HistoryGenerator,
    lambda event: frozenset({HistoryBuilder(Node(event))}),
)

HistoryProjector = projector_of(HistoryBuilder(Node()))
history_projector = Structure(
    HistoryProjector,
    lambda event: HistoryBuilder(Node(event)),
)


class PartialHistoryBuilder:
    def __init__(self, value):
        self.value = value

    def __eq__(self, other):
        return self.value == other.value

    def __hash__(self):
        return hash(self.value)

    def __mul__(node1, node2):
        if node1.value.data is None or not node1.value.data.apparent:
            return node2

        if node2.value.data is None:
            return node1

        return PartialHistoryBuilder(node1.value.add(node2.value))


PartialHistoryGenerator = generator_of(PartialHistoryBuilder(Node()))
partial_history_generator = Structure(
    PartialHistoryGenerator,
    lambda event: frozenset({PartialHistoryBuilder(Node(event))}),
)

PartialHistoryProjector = projector_of(PartialHistoryBuilder(Node()))
partial_history_projector = Structure(
    PartialHistoryProjector,
    lambda event: PartialHistoryBuilder(Node(event)),
)

history_counter = Structure(Counter, lambda event: 1)
