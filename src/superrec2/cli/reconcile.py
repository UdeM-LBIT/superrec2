"""Compute a minimum-cost (super-)reconciliation of two trees."""
import argparse
import inspect
import json
import textwrap
import sys
from .util import add_arg_input, add_arg_output
from superrec2.model.reconciliation import (
    NodeEvent,
    EdgeEvent,
    get_default_cost,
    ReconciliationInput,
    SuperReconciliationInput,
    ReconciliationOutput,
)
from superrec2.compute.exhaustive import reconcile_exhaustive
from superrec2.compute.reconciliation import (
    reconcile_lca,
    reconcile_thl,
)
from superrec2.compute.super_reconciliation import (
    sreconcile_base_spfs,
    sreconcile_extended_spfs,
)
from superrec2.compute.unordered_super_reconciliation import (
    usreconcile_base_uspfs,
    usreconcile_extended_uspfs,
)
from superrec2.utils.dynamic_programming import RetentionPolicy


algorithms = {
    "exh": reconcile_exhaustive,
    "lca": reconcile_lca,
    "thl": reconcile_thl,
    "base_spfs": sreconcile_base_spfs,
    "ext_spfs": sreconcile_extended_spfs,
    "base_uspfs": usreconcile_base_uspfs,
    "superdtl": usreconcile_extended_uspfs,
}


cost_events = {
    NodeEvent.SPECIATION: ("spe", "a speciation"),
    NodeEvent.DUPLICATION: ("dup", "a duplication"),
    NodeEvent.HORIZONTAL_TRANSFER: ("hgt", "an horizontal transfer"),
    EdgeEvent.FULL_LOSS: ("floss", "a full loss"),
    EdgeEvent.SEGMENTAL_LOSS: ("sloss", "a segmental loss"),
}


def eval_cost(cost):
    """Evaluate a cost expression in the context of this module."""
    return eval(cost)  # pylint: disable=eval-used


def read_input(args):
    """Read and parse the input data file."""
    data = json.load(args.input)
    data["costs"] = dict(
        (kind, getattr(args, f"cost_{argname}"))
        for kind, (argname, _) in cost_events.items()
    )

    if "leaf_syntenies" in data:
        rec_input = SuperReconciliationInput.from_dict(data)
    else:
        rec_input = ReconciliationInput.from_dict(data)

    return rec_input


def call_algorithm(args, rec_input):
    """Execute the specified algorithm on the input data."""
    algo = algorithms[args.algorithm]
    algo_signature = inspect.signature(algo)
    params = list(algo_signature.parameters.values())

    if params[0].annotation != type(rec_input):
        if params[0].annotation == ReconciliationInput:
            print(
                f"Warning: '{args.algorithm}' is not a super-reconciliation \
algorithm: declared leaf syntenies will be ignored",
                file=sys.stderr,
            )
        else:
            print(
                f"Error: '{args.algorithm}' is a super-reconciliation \
algorithm: you need to provide leaf syntenies",
                file=sys.stderr,
            )
            return None

    if len(params) == 1:
        output = algo(rec_input)
    elif len(params) == 2 and params[1].annotation == RetentionPolicy:
        output = algo(
            rec_input,
            getattr(RetentionPolicy, args.solutions.upper()),
        )
    else:
        return None

    if output is None:
        results = []
    elif isinstance(output, ReconciliationOutput):
        results = [output]
    else:
        results = list(output)

    if not results:
        return None

    print("Minimum cost:", results[0].cost(), file=sys.stderr)
    return results


def dump_results(args, results):
    """Write reconciliation results to output file."""
    for result in results:
        json.dump(result.to_dict(), args.output)
        print(file=args.output)


def reconcile(args):
    """Run the reconcile subcommand with the given arguments."""
    rec_input = read_input(args)
    results = call_algorithm(args, rec_input)

    if results is None:
        return 1

    return dump_results(args, results)


def add_args(parser):
    """Add the reconcile subcommand to a command-line argument parser."""
    desc = __doc__
    desc += "\n\navailable algorithms:\n\n"
    indent = 22

    for key, impl in algorithms.items():
        paragraphs = textwrap.dedent(impl.__doc__).split("\n\n")
        unwrapped = " ".join(paragraphs[0].split("\n")).strip()
        rewrapped = textwrap.indent(textwrap.fill(unwrapped, 70), " " * indent)
        desc += "  " + key + " " * (indent - 2 - len(key))
        desc += rewrapped.strip() + "\n\n"

    subparser = parser.add_parser(
        "reconcile",
        description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_arg_input(subparser, "file defining the reconciliation problem")
    add_arg_output(subparser, "file where the reconciliation results will be stored")
    subparser.add_argument(
        "algorithm",
        metavar="ALGO",
        help="reconciliation algorithm to use (see list above)",
        choices=set(algorithms.keys()),
    )
    subparser.add_argument(
        "--solutions",
        metavar="POLICY",
        default="any",
        help="whether to print all minimum-cost solutions or any such \
solution (choices: all, any; default: any)",
        choices=("any", "all"),
    )

    for kind, (argname, fullname) in cost_events.items():
        subparser.add_argument(
            f"--cost-{argname}",
            metavar="COST",
            default=get_default_cost()[kind],
            type=eval_cost,
            help=f"cost of {fullname} event (default: %(default)s)",
        )

    subparser.set_defaults(func=reconcile)
