from ..utils.toposort import toposort_all
from ..utils.min_sequence import ExtendedIntegral, MinSequence
from ..utils.lowest_common_ancestor import LowestCommonAncestor
from ..utils.subsequences import (
    subseq_complete,
    mask_from_subseq,
    subseq_from_mask,
    subseq_segment_dist,
)
from .tools import Labeling, Event, get_event, Reconciliation
from collections import defaultdict
from ete3 import PhyloTree, PhyloNode
from infinity import inf
from typing import (
    Any,
    DefaultDict,
    Dict,
    NamedTuple,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
)


class LabelingInfo(NamedTuple):
    cost: ExtendedIntegral = inf
    left_mask: Optional[int] = None
    right_mask: Optional[int] = None


def _compute_spfs_table(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
    known_syntenies: Labeling,
) -> Dict[PhyloNode, DefaultDict[int, LabelingInfo]]:
    root_synt = known_syntenies[gene_tree]
    subseq_count = 2 ** len(root_synt)
    costs: Dict[PhyloNode, DefaultDict[int, LabelingInfo]] = {}

    for sub_gene in gene_tree.traverse("postorder"):
        costs[sub_gene] = defaultdict(LabelingInfo)

        if sub_gene.is_leaf():
            sub_synt = known_syntenies[sub_gene]
            costs[sub_gene][
                mask_from_subseq(sub_synt, root_synt)
            ] = LabelingInfo(cost=0)
        else:
            left_gene, right_gene = sub_gene.children
            event = get_event(sub_gene, species_lca, rec)
            sub_synt_masks = (
                (subseq_complete(root_synt),)
                if sub_gene == gene_tree
                else range(1, subseq_count)
            )

            # Test all possible subsequences for the current node’s synteny
            for sub_synt_mask in sub_synt_masks:
                # Find optimal syntenies for the node’s two children
                LRSubsynt = Tuple[MinSequence[int], MinSequence[int]]
                min_synt_full: LRSubsynt = (MinSequence(1), MinSequence(1))
                min_synt_segm: LRSubsynt = (MinSequence(1), MinSequence(1))

                for i in (0, 1):
                    gene = sub_gene.children[i]

                    for synt, info in costs[gene].items():
                        full = subseq_segment_dist(
                            synt, sub_synt_mask, edges=True
                        )

                        if full == -1:
                            # Not a subsequence of the parent synteny
                            continue

                        segm = subseq_segment_dist(
                            synt, sub_synt_mask, edges=False
                        )

                        min_synt_full[i].update((info.cost + full, synt))
                        min_synt_segm[i].update((info.cost + segm, synt))

                if any(min_synt.min == inf for min_synt in min_synt_full):
                    continue

                keep_left_cost = min_synt_full[0].min + min_synt_segm[1].min
                keep_right_cost = min_synt_segm[0].min + min_synt_full[1].min

                if event == Event.Speciation:
                    costs[sub_gene][sub_synt_mask] = LabelingInfo(
                        cost=min_synt_full[0].min + min_synt_full[1].min,
                        left_mask=min_synt_full[0][0],
                        right_mask=min_synt_full[1][0],
                    )
                elif event == Event.Duplication:
                    if keep_left_cost <= keep_right_cost:
                        costs[sub_gene][sub_synt_mask] = LabelingInfo(
                            cost=keep_left_cost,
                            left_mask=min_synt_full[0][0],
                            right_mask=min_synt_segm[1][0],
                        )
                    else:
                        costs[sub_gene][sub_synt_mask] = LabelingInfo(
                            cost=keep_right_cost,
                            left_mask=min_synt_segm[0][0],
                            right_mask=min_synt_full[1][0],
                        )
                elif event == Event.HorizontalGeneTransfer:
                    keep_left = species_lca.is_comparable(
                        rec[sub_gene], rec[left_gene]
                    )
                    if keep_left:
                        costs[sub_gene][sub_synt_mask] = LabelingInfo(
                            cost=keep_left_cost,
                            left_mask=min_synt_full[0][0],
                            right_mask=min_synt_segm[1][0],
                        )
                    else:
                        costs[sub_gene][sub_synt_mask] = LabelingInfo(
                            cost=keep_right_cost,
                            left_mask=min_synt_segm[0][0],
                            right_mask=min_synt_full[1][0],
                        )
                else:
                    raise ValueError("Invalid event")

    return costs


def _decode_spfs_table(
    sub_gene: PhyloTree,
    root_synteny: Sequence[Any],
    sub_synteny: int,
    costs: Dict[PhyloNode, DefaultDict[int, LabelingInfo]],
) -> Labeling:
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
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
    known_syntenies: Labeling,
) -> Tuple[int, Labeling]:
    costs = _compute_spfs_table(gene_tree, species_lca, rec, known_syntenies)

    root_synteny = known_syntenies[gene_tree]
    all_mask = subseq_complete(root_synteny)
    return (
        costs[gene_tree][all_mask][0],
        _decode_spfs_table(gene_tree, root_synteny, all_mask, costs),
    )


def label_ancestral_syntenies(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
    leaf_labeling: Labeling,
) -> Tuple[ExtendedIntegral, List[Labeling]]:
    """
    Find a minimum-cost ancestral synteny labeling for a super-reconciliation.

    The :param:`leaf_labeling` must assign a synteny to each leaf of the
    gene tree. It can also assign a synteny to the root node. If no assignment
    of the root synteny is defined, this function will consider all possible
    orderings of the root synteny.

    The cost of a labeling is the total number of segments that are lost from
    each synteny to its children.

    :param gene_tree: gene tree to reconcile
    :param species_lca: species ancestry information
    :param rec: mapping of the gene tree onto the species tree
    :param leaf_labeling: syntenies assigned to the leaves of the gene tree
    :returns: a tuple containing the minimum cost of a labeling and a
        list of all labelings with such a cost
    """
    results: MinSequence[Labeling] = MinSequence()

    if gene_tree not in leaf_labeling:
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
                    gene_tree,
                    species_lca,
                    rec,
                    {**leaf_labeling, gene_tree: order},
                )
            )
    else:
        results.update(
            _label_with_root_order(
                gene_tree,
                species_lca,
                rec,
                leaf_labeling,
            )
        )

    return (results.min, list(results))
