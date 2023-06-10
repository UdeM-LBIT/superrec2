"""Compute reconciliations by iterating the complete search space."""
from itertools import product
from typing import Generator
from ete3 import TreeNode
from ..model.reconciliation import ReconciliationInput, ReconciliationOutput
from ..utils.dynamic_programming import (
    Candidate,
    Entry,
    MergePolicy,
    RetentionPolicy,
)


def generate_all(
    rec_input: ReconciliationInput, node: TreeNode = None
) -> Generator[ReconciliationOutput, None, None]:
    """
    Generate all possible outputs for a reconciliation input.

    :param rec_input: reconciliation input
    :returns: all possible reconciliations
    """
    if node is None:
        node = rec_input.object_tree

    if node.is_leaf():
        yield ReconciliationOutput(
            rec_input, {node: rec_input.leaf_object_species[node]}
        )
        return

    species_lca = rec_input.species_lca
    left_node, right_node = node.children

    for map_left, map_right in product(
        generate_all(rec_input, left_node),
        generate_all(rec_input, right_node),
    ):
        left_species = map_left.object_species[left_node]
        right_species = map_right.object_species[right_node]
        lca = species_lca(left_species, right_species)

        parent_species = lca
        while parent_species is not None:
            yield ReconciliationOutput(
                rec_input,
                {
                    node: parent_species,
                    **map_left.object_species,
                    **map_right.object_species,
                },
            )
            parent_species = parent_species.up

        for transfer_target, other_target in (
            (left_species, right_species),
            (right_species, left_species),
        ):
            if species_lca.is_ancestor_of(other_target, transfer_target):
                continue

            transfer_species = transfer_target

            while transfer_species != lca:
                yield ReconciliationOutput(
                    rec_input,
                    {
                        node: transfer_species,
                        **map_left.object_species,
                        **map_right.object_species,
                    },
                )
                transfer_species = transfer_species.up


def reconcile_exhaustive(
    rec_input: ReconciliationInput,
    policy: RetentionPolicy,
) -> Entry[int, ReconciliationOutput]:
    """
    Compute a minimum-cost reconciliation by enumerating all possible
    reconciliations.

    :param rec_input: reconciliation input
    :returns: any minimum-cost reconciliation, if there is one
    """
    results: Entry[int, ReconciliationOutput] = Entry(MergePolicy.MIN, policy)

    for output in generate_all(rec_input):
        results.update(Candidate(output.cost(), output))

    return results.infos()
