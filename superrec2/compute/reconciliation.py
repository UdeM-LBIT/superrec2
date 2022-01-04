"""Compute reconciliations with arbitrary event costs."""
from itertools import product
from typing import Generator, NamedTuple
from ete3 import TreeNode
from ..utils.trees import LowestCommonAncestor
from ..utils.dynamic_programming import (
    Candidate, DictDimension, MergePolicy, RetentionPolicy, Table
)
from ..model.reconciliation import (
    ReconciliationInput,
    ReconciliationOutput,
    NodeEvent,
    EdgeEvent,
    CostValues,
)


def reconcile_lca(rec_input: ReconciliationInput) -> ReconciliationOutput:
    """
    Compute the minimum-cost reconciliation between two trees,
    allowing for duplications and losses and where a duplication
    costs the same as a loss.

    :param rec_input: input for the reconciliation
    :returns: reconciliation result
    """
    rec = {}

    for node in rec_input.object_tree.traverse("postorder"):
        if node.is_leaf():
            species = rec_input.leaf_object_species[node]
            rec[node] = species
        else:
            left, right = node.children
            rec[node] = rec_input.species_lca(rec[left], rec[right])

    return ReconciliationOutput(rec_input, rec)


class MappingInfo(NamedTuple):
    """Information about an assignment of a node’s children to species."""
    # Species assigned to the left child
    left: TreeNode

    # Species assigned to the right child
    right: TreeNode


THLTable = Table[MappingInfo, int]


def _compute_thl_try_speciation(
    species_lca: LowestCommonAncestor,
    root_species: TreeNode,
    root_node: TreeNode,
    table: THLTable,
    costs: CostValues,
) -> None:
    loss_cost = costs[EdgeEvent.FULL_LOSS]

    left_species, right_species = root_species.children
    left_node, right_node = root_node.children

    # Optimal costs obtained by mapping the left or right node below
    # the left or right species
    min_ltl = table.entry()
    min_rtl = table.entry()
    min_ltr = table.entry()
    min_rtr = table.entry()

    for left_child in left_species.traverse():
        min_ltl.update(Candidate(
            table[left_node][left_child].value(),
            left_child,
        ))
        min_rtl.update(Candidate(
            table[right_node][left_child].value(),
            left_child,
        ))

    for right_child in right_species.traverse():
        min_ltr.update(Candidate(
            table[left_node][right_child].value(),
            right_child
        ))
        min_rtr.update(Candidate(
            table[right_node][right_child].value(),
            right_child
        ))

    spe_combinator = lambda left, right: Candidate(
        left.value + right.value
        + loss_cost * (
            species_lca.distance(root_species, left.info)
            + species_lca.distance(root_species, right.info)
            - 2
        ),
        MappingInfo(left.info, right.info),
    )

    table[root_node][root_species].update(
        *min_ltl.combine(min_rtr, spe_combinator),
        *min_ltr.combine(min_rtl, spe_combinator),
    )


def _compute_thl_try_duplication_transfer(
    species_lca: LowestCommonAncestor,
    root_species: TreeNode,
    root_node: TreeNode,
    table: THLTable,
    costs: CostValues,
) -> None:
    dup_cost = costs[NodeEvent.DUPLICATION]
    hgt_cost = costs[NodeEvent.HORIZONTAL_TRANSFER]
    loss_cost = costs[EdgeEvent.FULL_LOSS]

    left_node, right_node = root_node.children

    # Optimal costs obtained by mapping the left or right node inside
    # root_species’ subtree or outside of it
    min_ltc = table.entry()
    min_lts = table.entry()
    min_rtc = table.entry()
    min_rts = table.entry()

    for other_species in species_lca.tree.traverse():
        if species_lca.is_ancestor_of(root_species, other_species):
            min_ltc.update(Candidate(
                table[left_node][other_species].value(),
                other_species
            ))
            min_rtc.update(Candidate(
                table[right_node][other_species].value(),
                other_species
            ))
        elif not species_lca.is_ancestor_of(other_species, root_species):
            min_lts.update(Candidate(
                table[left_node][other_species].value(),
                other_species
            ))
            min_rts.update(Candidate(
                table[right_node][other_species].value(),
                other_species
            ))

    # Try mapping as a duplication
    dup_combin = lambda left, right: Candidate(
        dup_cost + left.value + right.value
        + loss_cost * (
            species_lca.distance(root_species, left.info)
            + species_lca.distance(root_species, right.info)
        ),
        MappingInfo(left.info, right.info),
    )


    # Try mapping as a horizontal transfer
    hgt_l_combin = lambda left, right: Candidate(
        hgt_cost + left.value + right.value
        + loss_cost * species_lca.distance(root_species, left.info),
        MappingInfo(left.info, right.info),
    )

    hgt_r_combin = lambda left, right: Candidate(
        hgt_cost + left.value + right.value
        + loss_cost * species_lca.distance(root_species, right.info),
        MappingInfo(left.info, right.info),
    )

    table[root_node][root_species].update(
        *min_ltc.combine(min_rtc, dup_combin),
        *min_lts.combine(min_rtc, hgt_r_combin),
        *min_ltc.combine(min_rts, hgt_l_combin),
    )


def _compute_thl_table(
    rec_input: ReconciliationInput,
    retention_policy: RetentionPolicy,
) -> THLTable:
    table = Table(
        (DictDimension(), DictDimension()),
        MergePolicy.MIN,
        retention_policy,
    )

    for root_node in rec_input.object_tree.traverse("postorder"):
        if root_node.is_leaf():
            root_species = rec_input.leaf_object_species[root_node]
            table[root_node][root_species] = Candidate(0)
        else:
            for root_species in rec_input.species_lca.tree.traverse("postorder"):
                if not root_species.is_leaf():
                    _compute_thl_try_speciation(
                        rec_input.species_lca,
                        root_species,
                        root_node,
                        table,
                        rec_input.costs
                    )

                _compute_thl_try_duplication_transfer(
                    rec_input.species_lca,
                    root_species,
                    root_node,
                    table,
                    rec_input.costs
                )

    return table


def _decode_thl_table(
    root_object: TreeNode,
    root_species: TreeNode,
    rec_input: ReconciliationInput,
    table: THLTable,
) -> Generator[ReconciliationOutput, None, None]:
    if not table[root_object][root_species].infos():
        yield ReconciliationOutput(
            rec_input, {root_object: root_species}
        )
        return

    for info in table[root_object][root_species].infos():
        left_object, right_object = root_object.children
        mappings = product(
            _decode_thl_table(
                left_object,
                info.left,
                rec_input,
                table,
            ),
            _decode_thl_table(
                right_object,
                info.right,
                rec_input,
                table,
            ),
        )

        for map_left, map_right in mappings:
            yield ReconciliationOutput(rec_input, {
                root_object: root_species,
                **map_left.object_species,
                **map_right.object_species,
            })


def reconcile_thl_any(rec_input: ReconciliationInput) -> ReconciliationOutput:
    """
    Find any minimum-cost reconciliation between two trees,
    allowing for duplications, horizontal transfers and losses.

    If you wish to disallow horizontal transfers, using `reconcile_lca` is
    faster than setting the cost of transfers to infinity.

    :param rec_input: reconciliation input
    :returns: any minimum-cost reconciliation
    """
    table = _compute_thl_table(rec_input, RetentionPolicy.ANY)

    root_object = rec_input.object_tree
    min_value, min_species = min(
        (table[root_object][root_species].value(), root_species)
        for root_species in rec_input.species_lca.tree.traverse()
    )

    for result in _decode_thl_table(root_object, min_species, rec_input, table):
        return result


def reconcile_thl_all(rec_input: ReconciliationInput) \
-> Generator[ReconciliationOutput, None, None]:
    """
    Find all minimum-cost reconciliations between two trees,
    allowing for duplications, horizontal transfers and losses.

    If you wish to disallow horizontal transfers, using `reconcile_lca` is
    faster than setting the cost of transfers to infinity.

    :param rec_input: reconciliation input
    :returns: generates all minimum-cost reconciliations
    """
    table = _compute_thl_table(rec_input, RetentionPolicy.ALL)

    root_object = rec_input.object_tree
    min_value = min(
        table[root_object][root_species].value()
        for root_species in rec_input.species_lca.tree.traverse()
    )

    for root_species in rec_input.species_lca.tree.traverse():
        if table[root_object][root_species].value() == min_value:
            yield from _decode_thl_table(
                root_object, root_species, rec_input, table
            )
