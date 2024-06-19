"""Compute a minimum-cost (super-)reconciliation of two trees."""

import argparse
import json
import inspect
from ast import literal_eval
from .util import add_arg_input, add_arg_output
from ..model.history import Reconciliation, graft_unsampled_hosts
from ..compute import synesth


methods = {}


def register_method(method):
    pretty_name = method.__name__.replace("_", "-")
    methods[pretty_name] = method


@register_method
def min_cost_single(setting, costs, output):
    """
    Report a single arbitrary minimum-cost history along with the total number
    of co-optimal solutions.
    """
    cost, count, history = synesth.min_cost_single(setting, costs)
    print(f"cost={cost}", file=output)
    print(f"count={count}", file=output)
    json.dump(history.to_mapping(), output)


@register_method
def min_cost_all(setting, costs, output):
    """Report all minimum-cost histories (may be exponentially slow!)."""
    cost, histories = synesth.min_cost_all(setting, costs)
    print(f"cost={cost}", file=output)
    print(f"count={len(histories)}", file=output)

    for history in histories:
        json.dump(history.to_mapping(), output)
        print(file=output)


@register_method
def pareto_single(setting, _, output):
    """
    Report Pareto-optimal event count vectors along with the number of
    histories having that event count and an arbitrarily-selected history
    for each event count vector.
    """
    result = synesth.pareto_single(setting)

    for key in sorted(result.keys(), key=tuple):
        count, history = result[key]
        print(f"events={key}", file=output)
        print(f"count={count}", file=output)
        json.dump(history.to_mapping(), output)
        print("\n", file=output)


@register_method
def pareto_all(setting, _, output):
    """
    Report Pareto-optimal event count vectors along with all histories
    having that event count vector (may be exponentially slow!).
    """
    result = synesth.pareto_all(setting)

    for key in sorted(result.keys(), key=tuple):
        histories = result[key]
        print(f"events={key}", file=output)
        print(f"count={len(histories)}", file=output)

        for history in histories:
            json.dump(history.to_mapping(), output)
            print(file=output)

        print(file=output)


def reconcile(args):
    """Run the reconcile subcommand with the given arguments."""
    costs_dict: dict[str, float] = {}

    if args.cost:
        for cost_list in args.cost:
            for cost_entry in cost_list:
                kind, value = cost_entry.split("=", maxsplit=1)
                costs_dict[kind.replace("-", "_")] = literal_eval(value)

    costs = synesth.EventCosts(**costs_dict)
    setting = Reconciliation.from_mapping(json.load(args.input))

    if args.allow_unsampled:
        setting = Reconciliation(
            host_tree=graft_unsampled_hosts(setting.host_tree),
            associate_tree=setting.associate_tree,
        )

    setting.validate()
    methods[args.method](setting, costs, args.output)


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

    event_types = [field.replace("_", "-") for field in synesth.EventCosts._fields]
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
        raw_doc = inspect.getdoc(method)

        if raw_doc is None:
            doc = "no documentation"
        else:
            doc = inspect.getdoc(method).replace("\n", " ")[:-1].lower()

        methods_help.append(name + " (" + doc + ")")

    subparser.add_argument(
        "--method",
        "-m",
        metavar="METHOD",
        default="min-cost-single",
        choices=list(methods.keys()),
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
