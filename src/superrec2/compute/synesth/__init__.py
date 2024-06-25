from sowing.node import Node
from typing import NamedTuple, TypeVar, Self, Mapping
from .recurrence import solve_binary
from .contents import propagate_contents
from ...utils.algebras import (
    vector,
    UnitMagma,
    SemiRing,
    Structure,
    Counter,
    MinPlus,
    pareto_of,
    generator_of,
    projector_of,
)
from ...model.history import (
    Reconciliation,
    History,
    Event,
    Codiverge,
    Diverge,
    Gain,
    Loss,
    Extant,
)


T = TypeVar("T")


def solve(setting: Reconciliation, structure: Structure[T, [Event]]) -> SemiRing[T]:
    """
    Solve a reconciliation problem over a given algebraic structure.

    :param setting: setting containing a host tree and an associate tree to reconcile
    :param structure: semiring structure with a morphism from events to the semiring
    :return: resulting semiring value
    """
    return sum(
        (
            solve_binary(binary_setting, structure)
            for binary_setting in setting.binarize()
        ),
        start=structure.zero,
    )


def make_history(event_tree: Node[Event, None], setting: Reconciliation) -> History:
    """Construct and validate a history object from a reconciliation result."""
    history = History(setting.host_tree, event_tree)
    history = propagate_contents(history.prune_unsampled())
    history.validate()
    return history


class EventCosts(NamedTuple):
    """Set of costs for the possible events in a reconciled history."""

    duplication: float = 1
    transfer_duplication: float = 1
    cut: float = 1
    transfer_cut: float = 1
    loss: float = 1

    def morphism(self, event: Event) -> float:
        """Map event instances to their cost according to the current set of costs."""
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


# Count the number of possible histories in the reconciliation recurrence
history_counter = Structure(Counter, lambda event: 1)


class HistoryBuilder(UnitMagma[Node[Event, None]]):
    """Assemble an event tree along the reconciliation recurrence."""

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
    """
    Assemble an event tree along the reconciliation recurrence,
    ignoring non-apparent nodes.
    """

    _one = Node()

    def _mul(node1, node2):
        def apparent(node: Node[Event, None]) -> bool:
            return node.data is not None and node.data.apparent

        if not apparent(node1) and not apparent(node2):
            return PartialHistoryBuilder._one

        if not apparent(node1):
            return node2

        if not apparent(node2):
            return node1

        return node1.add(node2)


def hom_partial_history(event: Event) -> PartialHistoryBuilder:
    if event.apparent:
        node = Node(event)
    else:
        node = Node()

    return PartialHistoryBuilder(node)


PartialHistoryGenerator = generator_of(PartialHistoryBuilder.one)
partial_history_generator = Structure(
    PartialHistoryGenerator,
    lambda event: frozenset({hom_partial_history(event)}),
)

PartialHistoryProjector = projector_of(PartialHistoryBuilder.one)
partial_history_projector = Structure(
    PartialHistoryProjector,
    lambda event: hom_partial_history(event),
)


def min_cost_single(
    setting: Reconciliation, costs: EventCosts
) -> tuple[int, int, History]:
    """
    Find a minimum-cost history for a reconciliation problem
    and count the number of co-optimal solutions.

    :param setting: setting containing a host tree and an associate tree to reconcile
    :param costs: set of costs to use
    :returns: minimum total cost, number of solution histories with that cost, and
        one of the histories with that cost
    """
    min_cost = Structure(MinPlus, costs.morphism)
    structure = min_cost * (history_counter + history_projector)
    result = solve(setting, structure).value

    cost, (count, solution) = result
    history = make_history(solution.value, setting)

    return cost, count, history


def min_cost_all(
    setting: Reconciliation, costs: EventCosts
) -> tuple[int, frozenset[History]]:
    """
    Find all minimum-cost histories for a reconciliation problem
    (may be exponentially slow!).

    :param setting: setting containing a host tree and an associate tree to reconcile
    :param costs: set of costs to use
    :returns: minimum total cost and set of histories with that cost
    """
    min_cost = Structure(MinPlus, costs.morphism)
    structure = min_cost * history_generator
    result = solve(setting, structure).value

    cost, solutions = result
    histories = frozenset(
        make_history(solution.value, setting) for solution in solutions
    )

    return cost, histories


@vector
class EventCounts(NamedTuple):
    """Vector counting the number of occurrences of each event type in a history."""

    duplication: int = 0
    transfer_duplication: int = 0
    cut: int = 0
    transfer_cut: int = 0
    loss: int = 0

    @staticmethod
    def morphism(event: Event) -> frozenset[Self]:
        match event:
            case Extant() | Gain() | Codiverge():
                return frozenset({EventCounts()})

            case Diverge():
                if not event.cut and not event.transfer:
                    return frozenset({EventCounts(duplication=1)})

                if not event.cut and event.transfer:
                    return frozenset({EventCounts(transfer_duplication=1)})

                if event.cut and not event.transfer:
                    return frozenset({EventCounts(cut=1)})

                if event.cut and event.transfer:
                    return frozenset({EventCounts(transfer_cut=1)})

            case Loss():
                return frozenset({EventCounts(loss=1)})

            case _:
                raise ValueError(f"unknown event type {type(event)}")


EventCountsPareto = pareto_of(EventCounts())
event_counts_pareto = Structure(EventCountsPareto, EventCounts.morphism)


def pareto_single(setting: Reconciliation) -> Mapping[EventCounts, tuple[int, History]]:
    """
    Find a history for each event count vector that is Pareto-optimal for a
    given reconciliation problem.

    :param setting: setting containing a host tree and an associate tree to reconcile
    :returns: mapping of each Pareto-optimal event count vector to a tuple with
        the number of histories with that event count vector and one of those histories
    """
    structure = event_counts_pareto @ (history_counter + history_projector)
    result = solve(setting, structure).value

    return {
        vector: (count, make_history(solution.value, setting))
        for vector, (count, solution) in result.items()
    }


def pareto_all(setting: Reconciliation) -> Mapping[EventCounts, frozenset[History]]:
    """
    Find all histories for each event count vector that is Pareto-optimal for a
    given reconciliation problem (may be exponentially slow!).

    :param setting: setting containing a host tree and an associate tree to reconcile
    :returns: mapping of each Pareto-optimal event count vector to the set
        of histories with that event count vector

    """
    structure = event_counts_pareto @ history_generator
    result = solve(setting, structure).value

    return {
        vector: frozenset(
            make_history(solution.value, setting) for solution in solutions
        )
        for vector, solutions in result.items()
    }
