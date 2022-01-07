#!/usr/bin/env python3
"""
Compute a minimum-cost (super-)reconciliation of two trees.

Input arguments OBJECT_TREE, SPECIES_TREE, and LEAF_SYNTENIES can be passed
either as plain strings or as a path to a file.

Available algorithms are:

"""
import argparse
import inspect
import textwrap
import sys
from infinity import inf
from ete3 import Tree
from superrec2.utils.trees import LowestCommonAncestor
from superrec2.utils.dynamic_programming import RetentionPolicy
from superrec2.model.synteny import parse_synteny_mapping
from superrec2.model.tree_mapping import (
    get_species_mapping, parse_tree_mapping
)
from superrec2.model.reconciliation import (
    NodeEvent, EdgeEvent, get_default_cost,
    ReconciliationInput, SuperReconciliationInput,
    ReconciliationOutput,
)
from superrec2.compute.exhaustive import reconcile_exhaustive
from superrec2.compute.reconciliation import (
    reconcile_lca, reconcile_thl,
)
from superrec2.compute.super_reconciliation import (
    sreconcile_base_spfs, sreconcile_extended_spfs,
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
        help="algorithm to use to reconcile",
        choices=set(algorithms.keys()),
    )
    parser.add_argument(
        "policy",
        metavar="POLICY",
        help="whether to generate any solution or all possible solutions",
        choices=("any", "all"),
    )
    parser.add_argument(
        "object_tree",
        metavar="OBJECT_TREE",
        help="object tree to embed (Newick format)",
    )
    parser.add_argument(
        "species_tree",
        metavar="SPECIES_TREE",
        help="host species tree (Newick format)",
    )
    parser.add_argument(
        "leaf_syntenies",
        metavar="LEAF_SYNTENIES",
        nargs="?",
        help="labelling of leaf objects with syntenies (if applicable)",
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
    object_tree = Tree(args.object_tree, format=1)
    species_tree = Tree(args.species_tree, format=1)
    costs = dict((
        (kind, getattr(args, f"cost_{argname}"))
        for kind, (argname, _) in cost_events.items()
    ))

    if args.leaf_syntenies is None:
        rec_input = ReconciliationInput(
            object_tree=object_tree,
            species_lca=LowestCommonAncestor(species_tree),
            leaf_object_species=get_species_mapping(object_tree, species_tree),
            costs=costs,
        )
    else:
        rec_input = SuperReconciliationInput(
            object_tree=object_tree,
            species_lca=LowestCommonAncestor(species_tree),
            leaf_object_species=get_species_mapping(
                object_tree, species_tree
            ),
            costs=costs,
            leaf_syntenies=parse_synteny_mapping(
                object_tree, args.leaf_syntenies
            ),
        )

    algo = algorithms[args.algorithm]
    algo_signature = inspect.signature(algo)
    params = list(algo_signature.parameters.values())

    if params[0].annotation != type(rec_input):
        if params[0].annotation == ReconciliationInput:
            print(f"Warning: '{args.algorithm}' is not a super-reconciliation \
algorithm: the LEAF_SYNTENIES argument will be ignored", file=sys.stderr)
        else:
            print(f"Error: '{args.algorithm}' is a super-reconciliation \
algorithm: you need to provide the LEAF_SYNTENIES argument", file=sys.stderr)
            return None

    if len(params) == 1:
        return algo(rec_input)

    if len(params) == 2 and params[1].annotation == RetentionPolicy:
        return algo(
            rec_input,
            RetentionPolicy.ALL if args.policy == "all"
            else RetentionPolicy.ANY
        )

    return None


def main():
    args = parse_arguments()
    result = call_algorithm(args)

    if result is None:
        return 1
    elif isinstance(result, ReconciliationOutput):
        print("Minimum cost:", result.cost(), file=sys.stderr)
        print(result)
    else:
        results = list(result)

        if not results:
            print("No solution", file=sys.stderr)
            return 1

        print("Minimum cost:", results[0].cost(), file=sys.stderr)
        for next_result in results:
            print(next_result)


if __name__ == "__main__":
    sys.exit(main())
