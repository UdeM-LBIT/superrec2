from sowing.node import Node
from typing import NamedTuple
from superrec2.utils.algebras import (
    vector,
    UnitMagma,
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
    duplication: float = 1
    transfer_duplication: float = 1
    cut: float = 1
    transfer_cut: float = 1
    loss: float = 1

    def event_cost_morphism(self, event: Event) -> float:
        match event:
            case Extant() | Gain() | Codiverge():
                return 0

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
    duplication: int = 0
    transfer_duplication: int = 0
    cut: int = 0
    transfer_cut: int = 0
    loss: int = 0


def event_vector_morphism(event: Event):
    match event:
        case Extant() | Gain() | Codiverge():
            return frozenset({EventVector()})

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


EventVectorPareto = pareto_of(EventVector())
event_vector_pareto = Structure(EventVectorPareto, event_vector_morphism)


class HistoryBuilder(UnitMagma[Node[Event, None]]):
    _one = Node()

    def _mul(node1, node2):
        if node1.data is None:
            return node2

        if node2.data is None:
            return node1

        return node1.add(node2)


HistoryGenerator = generator_of(HistoryBuilder.one)
history_generator = Structure(
    HistoryGenerator,
    lambda event: frozenset({HistoryBuilder(Node(event))}),
)

HistoryProjector = projector_of(HistoryBuilder.one)
history_projector = Structure(
    HistoryProjector,
    lambda event: HistoryBuilder(Node(event)),
)


class PartialHistoryBuilder(UnitMagma[Node[Event, None]]):
    _one = Node()

    def _mul(node1, node2):
        if node1.data is None or not node1.data.apparent:
            return node2

        if node2.data is None:
            return node1

        return node1.add(node2)


PartialHistoryGenerator = generator_of(PartialHistoryBuilder.one)
partial_history_generator = Structure(
    PartialHistoryGenerator,
    lambda event: frozenset({PartialHistoryBuilder(Node(event))}),
)

PartialHistoryProjector = projector_of(PartialHistoryBuilder.one)
partial_history_projector = Structure(
    PartialHistoryProjector,
    lambda event: PartialHistoryBuilder(Node(event)),
)

history_counter = Structure(Counter, lambda event: 1)
