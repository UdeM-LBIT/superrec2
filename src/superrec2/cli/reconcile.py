"""Compute a minimum-cost (super-)reconciliation of two trees."""
import argparse
import json
import inspect
from ast import literal_eval
from tqdm import tqdm
from functools import partial
from pathos.multiprocessing import Pool
from .util import add_arg_input, add_arg_output
from ..model.history import Reconciliation, History, graft_unsampled_hosts
from ..utils.algebras import make_single_selector, make_multiple_selector, make_product
from ..compute.util import (
    EventCosts,
    make_cost_algebra,
    event_vector_pareto,
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
def single_solution(run, cost_algebra, output):
    """
    Report a single arbitrary minimum-cost solution
    along with the total number of optimal solutions.
    """
    single_solution_algebra = make_single_selector(
        "single_solution_algebra",
        cost_algebra,
        make_product("history_count_unit_gen", history_counter, history_unit_generator),
    )

    result = run(structure=single_solution_algebra)
    print(f"cost={result.key.value}", file=output)
    print(f"count={result.value.value[0].value}", file=output)
    yield result.value.value[1].value.value


@register_method
def all_solutions(run, cost_algebra, output):
    """Report all minimum-cost solutions."""
    all_solutions_algebra = make_single_selector(
        "all_solutions_algebra",
        cost_algebra,
        history_generator,
    )

    result = run(structure=all_solutions_algebra)
    print(f"cost={result.key.value}", file=output)
    print(f"count={len(result.value.value)}", file=output)

    for solution in result.value.value:
        yield solution.value


@register_method
def pareto(run, cost_algebra, output):
    """
    Compute all Pareto-optimal event count vectors and
    the number of corresponding solutions for each vector.
    """
    event_vector_selector = make_multiple_selector(
        "event_vector_selector",
        event_vector_pareto,
        history_counter,
    )
    result = run(structure=event_vector_selector)

    for key in sorted(result.keys(), key=tuple):
        print(f"{key}: {result[key].value}")

    return
    yield


def reconcile(args):
    """Run the reconcile subcommand with the given arguments."""
    costs_dict: dict[str, float] = {}

    if args.cost:
        for cost_list in args.cost:
            for cost_entry in cost_list:
                kind, value = cost_entry.split("=", maxsplit=1)
                costs_dict[kind.replace("-", "_")] = literal_eval(value)

    costs = EventCosts(**costs_dict)
    cost_algebra = make_cost_algebra("cost_algebra", costs)
    setting = Reconciliation.from_mapping(json.load(args.input))

    if args.allow_unsampled:
        setting = Reconciliation(
            host_tree=graft_unsampled_hosts(setting.host_tree),
            associate_tree=setting.associate_tree,
        )

    setting.validate()

    for event_tree in methods[args.method](
        partial(
            superdtlx.reconcile,
            setting=setting,
            progress=tqdm,
            pool=Pool(args.processes),
        ),
        cost_algebra,
        args.output,
    ):
        history = History(setting.host_tree, event_tree)
        json.dump(superdtlx.finalize_history(history).to_mapping(), args.output)
        print(file=args.output)


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
