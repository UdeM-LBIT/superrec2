"""Compute and represent synteny-labeled reconciliations."""
from collections import defaultdict
from typing import (
    Any,
    DefaultDict,
    Dict,
    NamedTuple,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)
from ete3 import PhyloTree, PhyloNode
from infinity import inf
from ..utils.toposort import toposort_all
from ..utils.min_sequence import ExtendedIntegral, MinSequence
from ..utils.trees import LowestCommonAncestor
from ..utils.subsequences import (
    subseq_complete,
    mask_from_subseq,
    subseq_from_mask,
    subseq_segment_dist,
)
from ..model.synteny import SyntenyMapping
from ..model.tree_mapping import TreeMapping
from .tools import Event, get_event


class LabelingInfo(NamedTuple):
    """Information about an assignment of a synteny to a node."""

    # Cost of this assignment
    cost: ExtendedIntegral = inf

    # Synteny assigned to the left child, if applicable
    left_mask: Optional[int] = None

    # Synteny assigned to the right child, if applicable
    right_mask: Optional[int] = None


SPFSTable = Dict[PhyloNode, DefaultDict[int, LabelingInfo]]
LRSubsynt = Tuple[MinSequence[int], MinSequence[int]]


def _compute_spfs_entry(  # pylint:disable=too-many-locals
    sub_gene: PhyloNode,
    species_lca: LowestCommonAncestor,
    rec: TreeMapping,
    sub_synt_mask: int,
    table: SPFSTable,
) -> None:
    """Find optimal syntenies for the node’s children given a parent synteny."""
    event = get_event(sub_gene, species_lca, rec)
    min_synt_full: LRSubsynt = (MinSequence(1), MinSequence(1))
    min_synt_segm: LRSubsynt = (MinSequence(1), MinSequence(1))

    for i in (0, 1):
        gene = sub_gene.children[i]

        for synt, info in table[gene].items():
            full = subseq_segment_dist(synt, sub_synt_mask, edges=True)

            if full == -1:
                # Not a subsequence of the parent synteny
                continue

            segm = subseq_segment_dist(synt, sub_synt_mask, edges=False)

            min_synt_full[i].update((info.cost + full, synt))
            min_synt_segm[i].update((info.cost + segm, synt))

    if any(min_synt.min == inf for min_synt in min_synt_full):
        return

    keep_left_cost = min_synt_full[0].min + min_synt_segm[1].min
    keep_right_cost = min_synt_segm[0].min + min_synt_full[1].min

    if event == Event.SPECIATION:
        table[sub_gene][sub_synt_mask] = LabelingInfo(
            cost=min_synt_full[0].min + min_synt_full[1].min,
            left_mask=min_synt_full[0][0],
            right_mask=min_synt_full[1][0],
        )
    elif event == Event.DUPLICATION:
        if keep_left_cost <= keep_right_cost:
            table[sub_gene][sub_synt_mask] = LabelingInfo(
                cost=keep_left_cost,
                left_mask=min_synt_full[0][0],
                right_mask=min_synt_segm[1][0],
            )
        else:
            table[sub_gene][sub_synt_mask] = LabelingInfo(
                cost=keep_right_cost,
                left_mask=min_synt_segm[0][0],
                right_mask=min_synt_full[1][0],
            )
    elif event == Event.HORIZONTAL_GENE_TRANSFER:
        keep_left = species_lca.is_comparable(
            rec[sub_gene], rec[sub_gene.children[0]]
        )
        if keep_left:
            table[sub_gene][sub_synt_mask] = LabelingInfo(
                cost=keep_left_cost,
                left_mask=min_synt_full[0][0],
                right_mask=min_synt_segm[1][0],
            )
        else:
            table[sub_gene][sub_synt_mask] = LabelingInfo(
                cost=keep_right_cost,
                left_mask=min_synt_segm[0][0],
                right_mask=min_synt_full[1][0],
            )
    else:
        raise ValueError("Invalid event")


def _compute_spfs_table(
    synteny_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: TreeMapping,
    known_syntenies: SyntenyMapping,
) -> SPFSTable:
    root_synt = known_syntenies[synteny_tree]
    subseq_count = 2 ** len(root_synt)
    table: SPFSTable = {}

    for sub_gene in synteny_tree.traverse("postorder"):
        table[sub_gene] = defaultdict(LabelingInfo)

        if sub_gene.is_leaf():
            sub_synt = known_syntenies[sub_gene]
            table[sub_gene][
                mask_from_subseq(sub_synt, root_synt)
            ] = LabelingInfo(cost=0)
        else:
            sub_synt_masks = (
                (subseq_complete(root_synt),)
                if sub_gene == synteny_tree
                else range(1, subseq_count)
            )

            # Test all possible subsequences for the current node’s synteny
            for sub_synt_mask in sub_synt_masks:
                _compute_spfs_entry(
                    sub_gene, species_lca, rec, sub_synt_mask, table
                )

    return table


def _decode_spfs_table(
    sub_gene: PhyloTree,
    root_synteny: Sequence[Any],
    sub_synteny: int,
    costs: SPFSTable,
) -> SyntenyMapping:
    resolv_synteny = subseq_from_mask(sub_synteny, root_synteny)

    if sub_gene.is_leaf():
        return {sub_gene: resolv_synteny}

    left_gene, right_gene = sub_gene.children
    info = costs[sub_gene][sub_synteny]
    assert info.left_mask is not None
    assert info.right_mask is not None

    return {
        sub_gene: resolv_synteny,
        **_decode_spfs_table(left_gene, root_synteny, info.left_mask, costs),
        **_decode_spfs_table(right_gene, root_synteny, info.right_mask, costs),
    }


def _label_with_root_order(
    synteny_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: TreeMapping,
    known_syntenies: SyntenyMapping,
) -> Tuple[int, SyntenyMapping]:
    costs = _compute_spfs_table(synteny_tree, species_lca, rec, known_syntenies)

    root_synteny = known_syntenies[synteny_tree]
    all_mask = subseq_complete(root_synteny)
    return (
        costs[synteny_tree][all_mask][0],
        _decode_spfs_table(synteny_tree, root_synteny, all_mask, costs),
    )


def label_ancestral_syntenies(
    synteny_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: TreeMapping,
    syntenies: Union[SyntenyMapping, Mapping[str, Sequence[Any]]],
) -> Tuple[ExtendedIntegral, List[SyntenyMapping]]:
    """
    Find a minimum-cost ancestral synteny labeling for a super-reconciliation.

    The :param:`leaf_labeling` must assign a synteny to each leaf of the
    gene tree. It can also assign a synteny to the root node. If no assignment
    of the root synteny is defined, this function will consider all possible
    orderings of the root synteny.

    The cost of a labeling is the total number of segments that are lost from
    each synteny to its children.

    :param synteny_tree: synteny tree to label
    :param species_lca: species ancestry information
    :param rec: mapping of the synteny tree onto the species tree
    :param syntenies: syntenies assigned to the leaves of the synteny tree
        (keys can either by synteny names or node instances)
    :returns: a tuple containing the minimum cost of a labeling and a
        list of all labelings with such a cost
    """
    results: MinSequence[Labeling] = MinSequence()
    leaf_labeling: Labeling = {
        (synteny_tree & key if isinstance(key, str) else key): value
        for key, value in syntenies.items()
    }

    if synteny_tree not in leaf_labeling:
        prec: Dict[Any, Set[Any]] = {}

        for leaf_synteny in leaf_labeling.values():
            for gene_1, gene_2 in zip(leaf_synteny[0:-1], leaf_synteny[1:]):
                if gene_1 not in prec:
                    prec[gene_1] = set()
                prec[gene_1].add(gene_2)

            if leaf_synteny[-1] not in prec:
                prec[leaf_synteny[-1]] = set()

        for order in toposort_all(prec):
            results.update(
                _label_with_root_order(
                    synteny_tree,
                    species_lca,
                    rec,
                    {**leaf_labeling, synteny_tree: order},
                )
            )
    else:
        results.update(
            _label_with_root_order(
                synteny_tree,
                species_lca,
                rec,
                leaf_labeling,
            )
        )

    return (results.min, list(results))
