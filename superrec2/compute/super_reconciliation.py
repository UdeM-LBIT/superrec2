"""Compute and represent synteny-labeled reconciliations."""
from itertools import product
from typing import (
    Any,
    Dict,
    Generator,
    Optional,
    NamedTuple,
    Sequence,
    Set,
)
from ete3 import TreeNode
from ..utils.toposort import toposort_all
from ..utils.trees import LowestCommonAncestor
from ..utils.subsequences import (
    subseq_complete,
    mask_from_subseq,
    subseq_from_mask,
    subseq_segment_dist,
)
from ..model.synteny import Synteny, SyntenyMapping
from ..model.reconciliation import (
    ReconciliationOutput,
    SuperReconciliationInput,
    SuperReconciliationOutput,
    NodeEvent,
    EdgeEvent,
    CostValues,
)
from ..utils.dynamic_programming import (
    Candidate, DictDimension, Entry, MergePolicy, RetentionPolicy, Table
)


class LabelingInfo(NamedTuple):
    """Information about an assignment of syntenies."""
    # Synteny assigned to the left child (as a mask)
    left: int

    # Synteny assigned to the right child (as a mask)
    right: int


SPFSTable = Table[LabelingInfo, int]


def _compute_spfs_entry(
    root_node: TreeNode,
    root_synt_mask: int,
    species_lca: LowestCommonAncestor,
    rec_output: ReconciliationOutput,
    table: SPFSTable,
    costs: CostValues,
) -> None:
    """Find optimal syntenies for the node’s children given a parent synteny."""
    event = rec_output.node_event(root_node)

    min_full = (table.entry(), table.entry())
    min_seg = (table.entry(), table.entry())

    for child_index in (0, 1):
        child_node = root_node.children[child_index]

        for child_synt in table[child_node]:
            full = subseq_segment_dist(child_synt, root_synt_mask, edges=True)

            if full == -1:
                # Not a subsequence of the parent synteny
                continue

            segm = subseq_segment_dist(child_synt, root_synt_mask, edges=False)

            full *= costs[EdgeEvent.SEGMENTAL_LOSS]
            segm *= costs[EdgeEvent.SEGMENTAL_LOSS]
            cost = table[child_node][child_synt].value()

            min_full[child_index].update(Candidate(cost + full, child_synt))
            min_seg[child_index].update(Candidate(cost + segm, child_synt))

    sum_combinator = lambda left, right: Candidate(
        left.value + right.value,
        LabelingInfo(left.info, right.info),
    )

    if event == NodeEvent.SPECIATION:
        table[root_node][root_synt_mask].update(*min_full[0].combine(
            min_full[1],
            sum_combinator
        ))
    elif event == NodeEvent.DUPLICATION:
        table[root_node][root_synt_mask].update(
            *min_full[0].combine(min_seg[1], sum_combinator),
            *min_seg[0].combine(min_full[1], sum_combinator),
        )
    elif event == NodeEvent.HORIZONTAL_TRANSFER:
        keep_left = species_lca.is_comparable(
            rec_output.object_species[root_node],
            rec_output.object_species[root_node.children[0]]
        )

        if keep_left:
            table[root_node][root_synt_mask].update(*min_full[0].combine(
                min_seg[1],
                sum_combinator
            ))
        else:
            table[root_node][root_synt_mask].update(*min_seg[0].combine(
                min_full[1],
                sum_combinator
            ))
    else:
        raise ValueError("Invalid event")


def _compute_spfs_table(
    srec_input: SuperReconciliationInput,
    root_synt: Synteny,
    rec_output: ReconciliationOutput,
    retention_policy: RetentionPolicy,
) -> SPFSTable:
    subseq_count = 2 ** len(root_synt)
    table = Table(
        (DictDimension(), DictDimension()),
        MergePolicy.MIN,
        retention_policy,
    )

    for root_node in srec_input.object_tree.traverse("postorder"):
        if root_node.is_leaf():
            mask = mask_from_subseq(
                srec_input.leaf_syntenies[root_node],
                root_synt
            )
            table[root_node][mask] = Candidate(0)
        else:
            root_synt_masks = (
                (subseq_complete(root_synt),)
                if root_node == srec_input.object_tree
                else range(1, subseq_count)
            )

            # Test all possible subsequences for the current node’s synteny
            for root_synt_mask in root_synt_masks:
                _compute_spfs_entry(
                    root_node,
                    root_synt_mask,
                    srec_input.species_lca,
                    rec_output,
                    table,
                    srec_input.costs,
                )

    return table


def _decode_spfs_table(
    current_node: TreeNode,
    current_mask: int,
    root_synteny: Sequence[Any],
    srec_input: SuperReconciliationInput,
    rec_output: ReconciliationOutput,
    table: SPFSTable,
) -> Generator[SuperReconciliationOutput, None, None]:
    resolv_synteny = subseq_from_mask(current_mask, root_synteny)

    if not table[current_node][current_mask].infos():
        yield SuperReconciliationOutput(
            input=srec_input,
            object_species=rec_output.object_species,
            syntenies={current_node: resolv_synteny}
        )
        return

    for info in table[current_node][current_mask].infos():
        left_gene, right_gene = current_node.children
        mappings = product(
            _decode_spfs_table(
                left_gene,
                info.left,
                root_synteny,
                srec_input,
                rec_output,
                table,
            ),
            _decode_spfs_table(
                right_gene,
                info.right,
                root_synteny,
                srec_input,
                rec_output,
                table,
            ),
        )

        for map_left, map_right in mappings:
            yield SuperReconciliationOutput(
                input=srec_input,
                object_species=rec_output.object_species,
                syntenies={
                    current_node: resolv_synteny,
                    **map_left.syntenies,
                    **map_right.syntenies,
                }
            )


def _make_prec_graph(leaf_syntenies: SyntenyMapping):
    prec: Dict[Any, Set[Any]] = {}

    for leaf_synteny in leaf_syntenies.values():
        for gene_1, gene_2 in zip(leaf_synteny[0:-1], leaf_synteny[1:]):
            if gene_1 not in prec:
                prec[gene_1] = set()
            prec[gene_1].add(gene_2)

        if leaf_synteny[-1] not in prec:
            prec[leaf_synteny[-1]] = set()

    return prec


def _label_syntenies(
    srec_input: SuperReconciliationInput,
    rec_output: ReconciliationOutput,
    retention_policy: RetentionPolicy,
):
    assert rec_output.input.object_tree == srec_input.object_tree
    assert rec_output.input.species_lca == srec_input.species_lca
    assert rec_output.input.leaf_object_species == srec_input.leaf_object_species

    synteny_tree = srec_input.object_tree
    leaf_syntenies = srec_input.leaf_syntenies

    if synteny_tree not in leaf_syntenies:
        root_synts = toposort_all(_make_prec_graph(leaf_syntenies))
    else:
        root_synts = (leaf_syntenies[synteny_tree],)

    results = Entry(MergePolicy.MIN, retention_policy)

    for root_synt in root_synts:
        table = _compute_spfs_table(
            srec_input,
            root_synt,
            rec_output,
            retention_policy,
        )

        results.update(*map(
            lambda output: Candidate(output.cost(), output),
            _decode_spfs_table(
                synteny_tree,
                subseq_complete(root_synt),
                root_synt,
                srec_input,
                rec_output,
                table,
            )
        ))

    return results.infos()


def label_syntenies_any(
    srec_input: SuperReconciliationInput,
    rec_output: ReconciliationOutput,
) -> Optional[SuperReconciliationOutput]:
    """
    Find any minimum-cost super-reconciliation based on a reconciliation.

    The `leaf_syntenies` attribute of :param:`rec_input` must assign a synteny
    to each leaf of the gene tree. It can also assign a synteny to the root
    node. If no assignment of the root synteny is defined, this function will
    consider all possible orderings of the root synteny.

    The cost of a labeling is the total number of segments that are lost from
    each synteny to its children.

    :param srec_input: objects of the super-reconciliation
    :param rec_output: computed reconciliation for the input
    :returns: any minimum-cost super-reconciliation, if there is one
    """
    results = _label_syntenies(srec_input, rec_output, RetentionPolicy.ANY)

    if results:
        return min(results)
    else:
        return None


def label_syntenies_all(
    srec_input: SuperReconciliationInput,
    rec_output: ReconciliationOutput,
) -> Set[SuperReconciliationOutput]:
    """
    Find all minimum-cost super-reconciliation based on a reconciliation.

    The `leaf_syntenies` attribute of :param:`rec_input` must assign a synteny
    to each leaf of the gene tree. It can also assign a synteny to the root
    node. If no assignment of the root synteny is defined, this function will
    consider all possible orderings of the root synteny.

    The cost of a labeling is the total number of segments that are lost from
    each synteny to its children.

    :param srec_input: objects of the super-reconciliation
    :param rec_output: computed reconciliation for the input
    :returns: generates all minimum-cost super-reconciliations
    """
    return _label_syntenies(srec_input, rec_output, RetentionPolicy.ALL)
