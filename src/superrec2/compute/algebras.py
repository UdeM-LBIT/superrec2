from sowing.node import Node
from collections.abc import Mapping
from superrec2.utils.algebras import make_unit_magma, make_generator, min_plus
from superrec2.model.history import Event, Codiverge, Diverge, Gain, Loss, Extant


def make_cost_algebra(typename: str, costs: Mapping[str, float]):
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
                    return costs.get("transfer_dup", 1)

                if event.cut and not event.transfer:
                    return costs.get("cut", 1)

                if event.cut and event.transfer:
                    return costs.get("transfer_cut", 1)

            case Loss():
                return costs.get("loss", 1)

        raise ValueError

    return min_plus(typename, make)


def join_event_nodes(node1: Node, node2: Node) -> Node:
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

history_generator = make_generator("history_generator", history_builder)
