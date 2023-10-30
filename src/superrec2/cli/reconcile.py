"""Compute a minimum-cost (super-)reconciliation of two trees."""
import argparse
import json
from ast import literal_eval
from tqdm import tqdm
from .util import add_arg_input, add_arg_output
from ..model.history import Reconciliation, History, graft_unsampled_hosts
from ..utils.algebras import make_selector, make_product
from ..compute.util import (
    make_cost_algebra,
    history_counter,
    history_unit_generator,
    history_generator,
)
from ..compute import superdtlx


methods = {}


def register_method(method):
    pretty_name = method.__name__.replace("_", "-")
    methods[pretty_name] = method


@register_method
def single_solution(algo, cost_algebra, setting, output):
    """Report a single arbitrary minimum-cost solution."""
    single_solution_algebra = make_selector(
        "single_solution_algebra",
        cost_algebra,
        make_product("history_count_unit_gen", history_counter, history_unit_generator),
    )

    result = algo.reconcile(setting, single_solution_algebra, progress=tqdm)
    print(f"cost={result.cost.value}", file=output)
    print(f"count={result.selected.value[0].value}", file=output)

    history = History(setting.host_tree, result.selected.value[1].value.value)
    json.dump(algo.finalize_history(history).to_mapping(), output)
    print(file=output)


@register_method
def all_solutions(algo, cost_algebra, setting, output):
    """Report all minimum-cost solutions."""
    all_solutions_algebra = make_selector(
        "all_solutions_algebra",
        cost_algebra,
        history_generator,
    )

    result = algo.reconcile(setting, all_solutions_algebra, progress=tqdm)
    print(f"cost={result.cost.value}", file=output)
    print(f"count={len(result.selected.value)}", file=output)

    for solution in result.selected.value:
        history = History(setting.host_tree, solution.value)
        json.dump(algo.finalize_history(history).to_mapping(), output)
        print(file=output)


def reconcile(args):
    """Run the reconcile subcommand with the given arguments."""
    costs = {}

    if args.cost:
        for cost_list in args.cost:
            for cost_entry in cost_list:
                kind, value = cost_entry.split("=", maxsplit=1)
                costs[kind] = literal_eval(value)

    cost_algebra = make_cost_algebra("cost_algebra", costs)
    setting = Reconciliation.from_mapping(json.load(args.input))

    if args.allow_unsampled:
        setting = Reconciliation(
            host_tree=graft_unsampled_hosts(setting.host_tree),
            associate_tree=setting.associate_tree,
        )

    setting.validate()
    methods[args.method](superdtlx, cost_algebra, setting, args.output)


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

    subparser.add_argument(
        "--cost",
        "-c",
        nargs="*",
        action="append",
        metavar="TYPE=VALUE",
        help=(
            "set the cost of an event type. repeat to set multiple costs. available "
            "event types: speciation, loss, dup, cut, transfer-dup, transfer-cut "
            "(default: use unit costs)"
        ),
    )

    methods_help = "select what to compute. available methods:"

    for name, method in methods.items():
        methods_help += " " + name + " (" + method.__doc__ + ")"

    subparser.add_argument(
        "--method",
        "-m",
        metavar="METHOD",
        default="single-solution",
        help=(
            "select what to compute (default: %(default)s). available methods: "
            + ", ".join(
                f"{name} ({method.__doc__.lower()[:-1]})"
                for name, method in methods.items()
            )
        ),
    )

    subparser.set_defaults(func=reconcile)
