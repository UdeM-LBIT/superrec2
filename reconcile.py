#!/usr/bin/env python3
"""Compute a minimum-cost (super-)reconciliation of two trees."""
import argparse
import inspect
import json
import textwrap
import sys
from infinity import inf
from ete3 import Tree
from superrec2.utils.trees import LowestCommonAncestor
from superrec2.utils.dynamic_programming import RetentionPolicy
from superrec2.model.synteny import parse_synteny_mapping
from superrec2.model.tree_mapping import get_species_mapping, parse_tree_mapping
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


algorithms = {
    "exh": reconcile_exhaustive,
    "lca": reconcile_lca,
    "thl": reconcile_thl,
    "base_spfs": sreconcile_base_spfs,
    "ext_spfs": sreconcile_extended_spfs,
}


cost_events = {
    NodeEvent.SPECIATION: ("spe", "a speciation"),
    NodeEvent.DUPLICATION: ("dup", "a duplication"),
    NodeEvent.HORIZONTAL_TRANSFER: ("hgt", "an horizontal transfer"),
    EdgeEvent.FULL_LOSS: ("floss", "a full loss"),
    EdgeEvent.SEGMENTAL_LOSS: ("sloss", "a segmental loss"),
}


# Append each algorithm’s description to this module help string
__doc__ += "\nAvailable algorithms:\n\n"

for key, impl in algorithms.items():
    paragraphs = textwrap.dedent(impl.__doc__).split("\n\n")
    unwrapped = " ".join(paragraphs[0].split("\n")).strip()
    rewrapped = textwrap.indent(textwrap.fill(unwrapped, 70), " " * 4)
    __doc__ += f"{key}\n{rewrapped}\n\n"


def eval_cost(cost):
    """Evaluate a cost expression in the context of this module."""
    return eval(cost)


def parse_arguments():
    """Retrieve arguments or show help."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "algorithm",
        metavar="ALGO",
        help="reconciliation algorithm to use",
        choices=set(algorithms.keys()),
    )
    parser.add_argument(
        "policy",
        metavar="POLICY",
        help="whether to generate any solution or all possible solutions",
        choices=("any", "all"),
    )
    parser.add_argument(
        "--input",
        metavar="PATH",
        default="-",
        help="path to a file defining the reconciliation input \
(default: read from stdin)",
    )
    parser.add_argument(
        "--output",
        metavar="PATH",
        default="-",
        help="path where the resulting reconciliations will be stored \
(default: output to stdout)",
    )
    for kind, (argname, fullname) in cost_events.items():
        parser.add_argument(
            f"--cost-{argname}",
            metavar="COST",
            default=get_default_cost()[kind],
            type=eval_cost,
            help=f"cost of {fullname} event (default: %(default)s)",
        )
    return parser.parse_args()


def call_algorithm(args):
    infile = open(args.input, "r") if args.input != "-" else sys.stdin
    outfile = open(args.output, "w") if args.output != "-" else sys.stdout
    data = json.load(infile)
    data["costs"] = dict(
        (kind, getattr(args, f"cost_{argname}"))
        for kind, (argname, _) in cost_events.items()
    )

    if "leaf_syntenies" in data:
        rec_input = SuperReconciliationInput.from_dict(data)
    else:
        rec_input = ReconciliationInput.from_dict(data)

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
            return 1

    if len(params) == 1:
        output = algo(rec_input)
    elif len(params) == 2 and params[1].annotation == RetentionPolicy:
        output = algo(
            rec_input,
            RetentionPolicy.ALL
            if args.policy == "all"
            else RetentionPolicy.ANY,
        )
    else:
        return 1

    if output is None:
        results = []
    elif isinstance(output, ReconciliationOutput):
        results = [output]
    else:
        results = list(output)

    if not results:
        return 1

    print("Minimum cost:", results[0].cost(), file=sys.stderr)

    for result in results:
        json.dump(result.to_dict(), outfile)
        print(file=outfile)

    return 0


def main():
    args = parse_arguments()
    return call_algorithm(args)


if __name__ == "__main__":
    sys.exit(main())
