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
from superrec2.model.synteny import parse_synteny_mapping
from superrec2.model.tree_mapping import (
    get_species_mapping, parse_tree_mapping
)
from superrec2.model.reconciliation import (
    NodeEvent, EdgeEvent, get_default_cost,
    ReconciliationInput, SuperReconciliationInput,
    ReconciliationOutput,
)
from superrec2.compute.exhaustive import (
    exhaustive_any, exhaustive_all,
)
from superrec2.compute.reconciliation import (
    reconcile_lca, reconcile_thl_any, reconcile_thl_all,
)
from superrec2.compute.super_reconciliation import (
    label_syntenies_any, label_syntenies_all,
)


algorithms = {
    "exh_any": exhaustive_any,
    "exh_all": exhaustive_all,
    "lca": reconcile_lca,
    "thl_any": reconcile_thl_any,
    "thl_all": reconcile_thl_all,
    "spfs_any": label_syntenies_any,
    "spfs_all": label_syntenies_all,
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
    input_param = list(algo_signature.parameters.values())[0]

    if input_param.annotation != type(rec_input):
        if input_param.annotation == ReconciliationInput:
            print(f"Warning: '{args.algorithm}' is not a super-reconciliation \
algorithm: the LEAF_SYNTENIES argument will be ignored")
        else:
            print(f"Error: '{args.algorithm}' is a super-reconciliation \
algorithm: you need to provide the LEAF_SYNTENIES argument")
            return None

    return algo(rec_input)


def main():
    args = parse_arguments()
    result = call_algorithm(args)

    if result is None:
        return 1
    elif isinstance(result, ReconciliationOutput):
        print(result)
    else:
        for next_result in result:
            print(next_result)


if __name__ == "__main__":
    sys.exit(main())
