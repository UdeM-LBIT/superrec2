#!/usr/bin/env python3
"""Test reconciliation algorithms on exhaustively-generated inputs."""
import argparse
from itertools import product
import math
from time import time
from datetime import timedelta
from string import ascii_lowercase
from typing import Any, Generator, List, Mapping, Sequence, Tuple, TypeVar
from ete3 import PhyloTree, PhyloNode
from superrec2.reconciliation.tools import (
    SuperReconciliation,
    SuperReconciliationInput,
    CostType,
    CostVector,
    get_species_name,
    get_reconciliation_cost,
    get_labeling_cost,
    reconcile_all,
)
from superrec2.reconciliation.compute import reconcile_thl
from superrec2.reconciliation.syntenies import label_ancestral_syntenies
from superrec2.utils.trees import LowestCommonAncestor
from superrec2.utils.min_sequence import MinSequence


def trees_all(labels: Sequence[Any]) -> Generator[PhyloTree, None, None]:
    """
    Enumerate all binary trees that display the given list of labels.

    :param labels: list of labels of the leaves
    :yields: generated binary trees
    """
    if len(labels) == 1:
        node = PhyloNode()
        node.name = labels[0]
        yield node
        return

    for i in range(1, len(labels)):
        for left, right in product(
            trees_all(labels[:i]),
            trees_all(labels[i:]),
        ):
            root = PhyloNode()
            root.add_child(left.copy())
            root.add_child(right.copy())
            root.name = left.name + right.name
            yield root


def trees_all_count(labels_count: int) -> int:
    """
    Count the number of binary trees that display a given number of labels.

    This is the (n-1)th Catalan number.

    :param labels_count: number of labels
    :returns: number of distinct binary trees
    """
    return math.comb(2 * (labels_count - 1), labels_count - 1) // labels_count


Element = TypeVar("Element")


def map_all(
    domain: List[Element], values: int, next_value: int = 0
) -> Generator[Mapping[Element, int], None, None]:
    """
    Enumerate all equivalent mappings of a set to another.

    Two mappings are said to be equivalent if the equivalence classes formed by
    their images are the same. For example, the two following mappings are
    equivalent, and only one of them will be listed by this function:

    * {0: 0, 1: 0, 2: 1, 3: 0}
    * {0: 1, 1: 1, 2: 0, 3: 1}

    :param domain: set of values to map
    :param values: number of possible values
    :param next_value: next available value (used for recursion)
    :yields: generated mappings
    """
    if not domain:
        yield {}
        return

    for value in range(next_value):
        for rest in map_all(domain[1:], values, next_value):
            yield {domain[0]: value, **rest}

    if next_value < values:
        for rest in map_all(domain[1:], values, next_value + 1):
            yield {domain[0]: next_value, **rest}


def map_all_count(domain_size: int, values: int, next_value: int = 0) -> int:
    """
    Count the number of equivalent mappings of a set to another.

    :param domain_size: number of items to map
    :param values: number of possible values
    :param next_value: next available value (used for recursion)
    :returns: number of equivalent mappings
    """
    if not domain_size:
        return 1

    result = 0

    for _ in range(next_value):
        result += map_all_count(domain_size - 1, values, next_value)

    if next_value < values:
        result += map_all_count(domain_size - 1, values, next_value + 1)

    return result


def subseq_all(
    sequence: Sequence[Element],
) -> Generator[List[Element], None, None]:
    """
    Enumerate all subsequences of a sequence.

    :param sequence: sequence of elements
    :yields: generated subsequences
    """
    if len(sequence) == 0:
        yield []
        return

    for subsequence in subseq_all(sequence[1:]):
        yield subsequence
        yield [sequence[0]] + subsequence


def subseq_all_count(subsequence_length: int) -> int:
    """
    Count the number of subsequences in a sequence.

    :param sequence_length: length of the sequence
    :returns: number of subsequences
    """
    return 2 ** subsequence_length


def dtl_inputs_all(
    species_count: int, synteny_count: int, synteny_size: int
) -> Generator[SuperReconciliationInput, None, None]:
    """
    Generate valid inputs to the DTL-Super-Reconciliation problem.

    :param species_count: number of extant species in the species trees
    :param synteny_count: number of extant syntenies in the synteny trees
    :param synteny_size: size of the root synteny
    :yields: generated inputs
    """
    species_names = [str(i) for i in range(species_count)]
    synteny_indices = [str(i) for i in range(synteny_count)]
    root_synteny = ascii_lowercase[:synteny_size]

    for species_tree in trees_all(species_names):
        species_lca = LowestCommonAncestor(species_tree)
        sub_synts = (subseq for subseq in subseq_all(root_synteny) if subseq)

        for leaf_labeling_tuple in product(sub_synts, repeat=synteny_count):
            for leaf_mapping in map_all(
                synteny_indices,
                species_count,
            ):
                synteny_names = [
                    f"{species}_{synteny}"
                    for synteny, species in leaf_mapping.items()
                ]

                for synteny_tree in trees_all(synteny_names):
                    synteny_tree.set_species_naming_function(get_species_name)
                    leaf_labeling = {
                        synteny_tree & synteny_names[i]: synteny_value
                        for i, synteny_value in enumerate(leaf_labeling_tuple)
                    }

                    yield SuperReconciliationInput(
                        species_lca=species_lca,
                        synteny_tree=synteny_tree,
                        leaf_labeling=leaf_labeling,
                    )


def dtl_inputs_all_count(
    species_count: int, synteny_count: int, synteny_size: int
) -> int:
    """
    Count the number of valid inputs to the DTL-Super-Reconciliation problem.

    :param species_count: number of extant species in the species trees
    :param synteny_count: number of extant syntenies in the synteny trees
    :param synteny_size: size of the root synteny
    :returns: number of valid inputs
    """
    return (
        trees_all_count(species_count)
        * trees_all_count(synteny_count)
        * ((subseq_all_count(synteny_size) - 1) ** synteny_count)
        * map_all_count(synteny_count, species_count)
    )


def solve_thl_spfs(
    dtl_input: SuperReconciliationInput, costs: CostVector
) -> Tuple[int, List[SuperReconciliation]]:
    """
    Solve a DTL-Super-Reconciliation problem by first computing an optimal
    THL-reconciliation and then an optimal SPFS labeling for this
    reconciliation.

    :param dtl_input: input problem
    :param costs: cost of each event
    :returns: cost of optimal result and list of super-reconciliations
    """
    species_lca = dtl_input.species_lca
    synteny_tree = dtl_input.synteny_tree
    leaf_labeling = dtl_input.leaf_labeling

    rec_cost, optimal_recs = reconcile_thl(synteny_tree, species_lca, costs)
    results: MinSequence[SuperReconciliation] = MinSequence()

    for optimal_rec in optimal_recs:
        local_rec_cost = get_reconciliation_cost(
            synteny_tree, species_lca, optimal_rec, costs
        )
        assert local_rec_cost == rec_cost

        labeling_cost, labelings = label_ancestral_syntenies(
            synteny_tree, species_lca, optimal_rec, leaf_labeling
        )

        for labeling in labelings:
            local_labeling_cost = get_labeling_cost(
                synteny_tree, species_lca, optimal_rec, labeling
            )
            assert local_labeling_cost == labeling_cost

        results.update(
            (
                rec_cost + costs[CostType.SEGMENTAL_LOSS] * labeling_cost,
                SuperReconciliation(
                    synteny_tree=synteny_tree,
                    reconciliation=optimal_rec,
                    labeling=labelings[0],
                ),
            )
        )

    return results.min, list(results)


def solve_bruteforce(
    dtl_input: SuperReconciliationInput, costs: CostVector
) -> Tuple[int, List[SuperReconciliation]]:
    """
    Solve a DTL-Super-Reconciliation problem by enumerating all possible
    DTL-Reconciliations, computing an optimal SPFS for each of those and
    returning the optimal results.

    :param dtl_input: input problem
    :param costs: cost of each event
    :returns: cost of optimal result and list of super-reconciliations
    """
    species_lca = dtl_input.species_lca
    synteny_tree = dtl_input.synteny_tree
    leaf_labeling = dtl_input.leaf_labeling

    results: MinSequence[SuperReconciliation] = MinSequence()

    for valid_rec in reconcile_all(synteny_tree, species_lca):
        rec_cost = get_reconciliation_cost(
            synteny_tree, species_lca, valid_rec, costs
        )

        labeling_cost, labelings = label_ancestral_syntenies(
            synteny_tree, species_lca, valid_rec, leaf_labeling
        )

        for labeling in labelings:
            local_labeling_cost = get_labeling_cost(
                synteny_tree, species_lca, valid_rec, labeling
            )
            assert local_labeling_cost == labeling_cost

        results.update(
            (
                rec_cost + costs[CostType.SEGMENTAL_LOSS] * labeling_cost,
                SuperReconciliation(
                    synteny_tree=synteny_tree,
                    reconciliation=valid_rec,
                    labeling=labelings[0],
                ),
            )
        )

    return results.min, list(results)


def parse_arguments():
    """Retrieve arguments or show help."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--duplication-cost",
        "-D",
        metavar="COST",
        type=float,
        default=1,
        help="cost of a segmental duplication event, or 'inf' to disallow "
        "duplications (default: %(default)s)",
    )

    parser.add_argument(
        "--transfer-cost",
        "-T",
        metavar="COST",
        type=float,
        default=1,
        help="cost of a horizontal gene transfer, or 'inf' to disallow "
        "transfers (default: %(default)s)",
    )

    parser.add_argument(
        "--full-loss-cost",
        "-L",
        metavar="COST",
        type=float,
        default=1,
        help="cost of a full synteny loss, or 'inf' to disallow "
        "full losses (default: %(default)s)",
    )

    parser.add_argument(
        "--segmental-loss-cost",
        "-S",
        metavar="COST",
        type=float,
        default=1,
        help="cost of a synteny segment loss, or 'inf' to disallow "
        "segmental losses (default: %(default)s)",
    )

    parser.add_argument(
        "--species-count",
        "-s",
        metavar="COUNT",
        type=int,
        default=4,
        help="number of extant species in the generated inputs "
        "(default: %(default)s)",
    )

    parser.add_argument(
        "--synteny-count",
        "-y",
        metavar="COUNT",
        type=int,
        default=4,
        help="number of extant syntenies in the generated inputs "
        "(default: %(default)s)",
    )

    parser.add_argument(
        "--synteny-size",
        "-z",
        metavar="COUNT",
        type=int,
        default=3,
        help="number of genes in the root synteny of each generated input "
        "(default: %(default)s)",
    )

    parser.add_argument(
        "--map-to-root",
        "-r",
        action="store_true",
        help="only consider reconciliations where the root synteny is "
        "mapped to the root species (default: %(default)s)",
    )

    return parser.parse_args()


def do_search(
    args, costs: CostVector, params
):  # pylint:disable=too-many-locals
    """Perform exhaustive search for counter-examples."""
    total_inputs = dtl_inputs_all_count(**params)
    done_inputs = 0
    start_time = time()

    for dtl_input in dtl_inputs_all(**params):
        thl_spfs_cost, thl_spfs_results = solve_thl_spfs(dtl_input, costs)
        bf_cost, bf_results = solve_bruteforce(dtl_input, costs)

        if bf_cost != thl_spfs_cost:
            thl_spfs_result = thl_spfs_results[0]

            if args.map_to_root:
                # Only consider this case if there is a reconciliation
                # that maps the root synteny to the root species
                for bf_result in bf_results:
                    if (
                        bf_result.reconciliation[dtl_input.synteny_tree]
                        == dtl_input.species_lca.tree
                    ):
                        break
                else:
                    continue
            else:
                bf_result = bf_results[0]

            thl_rec_cost = get_reconciliation_cost(
                dtl_input.synteny_tree,
                dtl_input.species_lca,
                thl_spfs_result.reconciliation,
                costs,
            )

            bf_rec_cost = get_reconciliation_cost(
                dtl_input.synteny_tree,
                dtl_input.species_lca,
                bf_result.reconciliation,
                costs,
            )

            print("\n---\n")
            print(
                f"THL-SPFS = {thl_spfs_cost} ({thl_rec_cost} + "
                f"{thl_spfs_cost - thl_rec_cost})\n"
            )
            print(thl_spfs_result)
            print(
                f"\nBF = {bf_cost} ({bf_rec_cost} + "
                f"{bf_cost - bf_rec_cost})"
            )
            print(bf_result)

        done_inputs += 1

        if done_inputs % 1000 == 0:
            elapsed = time() - start_time
            eta = elapsed * (total_inputs / done_inputs - 1)

            print(
                f"Processed {done_inputs} of {total_inputs} "
                f"in {timedelta(seconds=round(elapsed))}s "
                f"(ETA {timedelta(seconds=round(eta))}s)"
            )


def main():  # pylint:disable=missing-function-docstring
    args = parse_arguments()

    costs = {
        CostType.DUPLICATION: args.duplication_cost,
        CostType.HORIZONTAL_GENE_TRANSFER: args.transfer_cost,
        CostType.FULL_LOSS: args.full_loss_cost,
        CostType.SEGMENTAL_LOSS: args.segmental_loss_cost,
    }

    params = {
        "species_count": args.species_count,
        "synteny_count": args.synteny_count,
        "synteny_size": args.synteny_size,
    }

    do_search(args, costs, params)


if __name__ == "__main__":
    main()
