"""Compute a minimum-cost (super-)reconciliation of two trees."""

import argparse
import json
import inspect
from ast import literal_eval
from functools import partial
from .util import add_arg_input, add_arg_output
from ..model.history import Reconciliation, History, graft_unsampled_hosts
from ..utils.algebras import Structure, MinPlus
from ..compute.util import (
    DummyProgress,
    DummyPool,
    EventCosts,
    event_vector_pareto,
    history_counter,
    history_projector,
    history_generator,
)
from ..compute import superdtlx


methods = {}


def register_method(method):
    pretty_name = method.__name__.replace("_", "-")
    methods[pretty_name] = method


def display_history(event_tree, setting, output):
    history = History(setting.host_tree, event_tree)
    history = superdtlx.finalize_history(history)
    history.validate()
    json.dump(history.to_mapping(), output)


@register_method
def single_solution(setting, costs, run, output):
    """
    Report a single arbitrary minimum-cost solution
    along with the total number of optimal solutions.
    """
    min_cost = Structure(MinPlus, costs.event_cost_morphism)
    structure = min_cost * (history_counter + history_projector)
    result = run(setting=setting, structure=structure)
    cost, (count, solution) = result

    print(f"cost={cost}", file=output)
    print(f"count={count}", file=output)
    display_history(solution.value, setting, output)


@register_method
def all_solutions(setting, costs, run, output):
    """Report all minimum-cost solutions."""
    min_cost = Structure(MinPlus, costs.event_cost_morphism)
    structure = min_cost * history_generator
    result = run(setting=setting, structure=structure)
    cost, solutions = result

    print(f"cost={cost}", file=output)
    print(f"count={len(solutions)}", file=output)

    for solution in solutions:
        display_history(solution.value, setting, output)
        print(file=output)


@register_method
def pareto(setting, _, run, output):
    """
    Compute all Pareto-optimal event count vectors and
    the number of corresponding solutions for each vector.
    """
    structure = event_vector_pareto @ history_counter
    result = run(setting=setting, structure=structure)

    for key in sorted(result.keys(), key=tuple):
        print(f"{key}: {result[key]}")


def reconcile(args):
    """Run the reconcile subcommand with the given arguments."""
    costs_dict: dict[str, float] = {}

    if args.cost:
        for cost_list in args.cost:
            for cost_entry in cost_list:
                kind, value = cost_entry.split("=", maxsplit=1)
                costs_dict[kind.replace("-", "_")] = literal_eval(value)

    costs = EventCosts(**costs_dict)
    setting = Reconciliation.from_mapping(json.load(args.input))

    if args.allow_unsampled:
        setting = Reconciliation(
            host_tree=graft_unsampled_hosts(setting.host_tree),
            associate_tree=setting.associate_tree,
        )

    setting.validate()

    methods[args.method](
        setting,
        costs,
        partial(
            superdtlx.reconcile,
            progress=DummyProgress,
            pool=DummyPool(),
            # FIXME: Restore multiprocess support
            # progress=tqdm,
            # pool=Pool(args.processes),
        ),
        args.output,
    )


def add_args(parser):
    """Add the reconcile subcommand to a command-line argument parser."""
    desc = __doc__

    subparser = parser.add_parser(
        "reconcile",
        description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    add_arg_input(subparser, "file defining the reconciliation problem")
    add_arg_output(subparser, "file where the results will be stored")

    subparser.add_argument(
        "--allow-unsampled",
        "-u",
        action="store_true",
        help="augment the host tree with candidate unsampled species",
    )

    event_types = [field.replace("_", "-") for field in EventCosts._fields]
    subparser.add_argument(
        "--cost",
        "-c",
        nargs="*",
        action="append",
        metavar="TYPE=VALUE",
        help=(
            "set the cost of an event type. repeat to set multiple costs. available "
            f"event types: {', '.join(event_types)} "
            "(default: use unit costs)"
        ),
    )

    methods_help = []

    for name, method in methods.items():
        doc = inspect.getdoc(method).replace("\n", " ")[:-1].lower()
        methods_help.append(name + " (" + doc + ")")

    subparser.add_argument(
        "--method",
        "-m",
        metavar="METHOD",
        default="single-solution",
        help=(
            "select what to compute (default: %(default)s). "
            "available methods: " + ", ".join(methods_help)
        ),
    )

    subparser.add_argument(
        "--processes",
        "-p",
        metavar="NPROC",
        default=1,
        type=int,
        help="number of processes to spawn for computing (default: %(default)s)",
    )

    subparser.set_defaults(func=reconcile)
