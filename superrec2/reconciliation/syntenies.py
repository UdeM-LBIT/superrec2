from ..utils.toposort import toposort_all
from ..utils.min_sequence import ExtendedIntegral, MinSequence
from ..utils.lowest_common_ancestor import LowestCommonAncestor
from ..utils.subsequences import (
    subseq_complete,
    mask_from_subseq,
    subseq_from_mask,
    subseq_segment_dist,
)
from .tools import Event, get_event, Reconciliation
from collections import defaultdict
from ete3 import PhyloTree, PhyloNode
from infinity import inf
from typing import Any, DefaultDict, Dict, List, Mapping, Sequence, Set, Tuple


Syntenies = Mapping[PhyloNode, Sequence[Any]]


def _compute_spfs_table(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
    known_syntenies: Syntenies,
) -> Dict[PhyloNode, DefaultDict[int, ExtendedIntegral]]:
    root_synteny = known_syntenies[gene_tree]
    subseq_count = 2 ** len(root_synteny)
    costs: Dict[DefaultDict[int, ExtendedIntegral]] = {}

    for sub_gene in gene_tree.traverse("postorder"):
        costs[sub_gene] = defaultdict(lambda: (inf, None, None))

        if sub_gene.is_leaf():
            sub_synteny = known_syntenies[sub_gene]
            costs[sub_gene][mask_from_subseq(sub_synteny, root_synteny)] = (
                0, None, None
            )
        else:
            left_gene, right_gene = sub_gene.children
            # event = get_event(sub_gene, species_lca, rec)
            sub_syntenies = (
                (subseq_count - 1,) if sub_gene == gene_tree
                else range(1, subseq_count)
            )

            for sub_synteny in sub_syntenies:
                min_syntenies = [None, None]
                min_costs = [inf, inf]

                for i in (0, 1):
                    gene = sub_gene.children[i]

                    for synteny, (cost, _, _) in costs[gene].items():
                        dist = subseq_segment_dist(
                            synteny, sub_synteny, edges=True
                        )

                        if dist == -1:
                            continue

                        cost += dist

                        if cost < min_costs[i]:
                            min_syntenies[i] = synteny
                            min_costs[i] = cost

                if inf not in min_costs:
                    costs[sub_gene][sub_synteny] = (
                        sum(min_costs),
                        *min_syntenies,
                    )

    return costs


def _decode_spfs_table(
    sub_gene: PhyloTree,
    root_synteny: Sequence[Any],
    sub_synteny: int,
    costs: Dict[PhyloNode, DefaultDict[int, ExtendedIntegral]],
) -> Syntenies:
    resolv_synteny = subseq_from_mask(sub_synteny, root_synteny)

    if sub_gene.is_leaf():
        return {sub_gene: resolv_synteny}

    left_gene, right_gene = sub_gene.children
    _, left_synteny, right_synteny = costs[sub_gene][sub_synteny]

    return {
        sub_gene: resolv_synteny,
        **_decode_spfs_table(left_gene, root_synteny, left_synteny, costs),
        **_decode_spfs_table(right_gene, root_synteny, right_synteny, costs),
    }


def _label_with_root_order(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
    known_syntenies: Syntenies,
) -> Tuple[int, Syntenies]:
    costs = _compute_spfs_table(gene_tree, species_lca, rec, known_syntenies)

    root_synteny = known_syntenies[gene_tree]
    all_mask = subseq_complete(root_synteny)
    return (
        costs[gene_tree][all_mask][0],
        _decode_spfs_table(
            gene_tree, root_synteny, all_mask, costs
        )
    )


def get_labeling_cost(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
    labeling: Syntenies,
) -> int:
    total_cost = 0
    root_syn = labeling[gene_tree]
    masks = {
        gene_tree: subseq_complete(root_syn)
    }

    for sub_gene in gene_tree.traverse("preorder"):
        if not sub_gene.is_leaf():
            event = get_event(sub_gene, species_lca, rec)
            sub_mask = masks[sub_gene]
            left_gene, right_gene = sub_gene.children

            left_syn = labeling[left_gene]
            left_mask = masks[left_gene] = \
                mask_from_subseq(left_syn, root_syn)

            right_syn = labeling[right_gene]
            right_mask = masks[right_gene] = \
                mask_from_subseq(right_syn, root_syn)

            if event == Event.Speciation:
                total_cost += (
                    subseq_segment_dist(left_mask, sub_mask, True)
                    + subseq_segment_dist(right_mask, sub_mask, True)
                )
            elif event == Event.Duplication:
                total_cost += min(
                    (
                        subseq_segment_dist(left_mask, sub_mask, True)
                        + subseq_segment_dist(left_mask, sub_mask, False)
                    ),
                    (
                        subseq_segment_dist(left_mask, sub_mask, False)
                        + subseq_segment_dist(left_mask, sub_mask, True)
                    )
                )
            elif event == Event.HorizontalGeneTransfer:
                left_cons = species_lca.is_comparable(
                    rec[sub_gene], rec[left_gene]
                )
                total_cost += (
                    subseq_segment_dist(left_mask, sub_mask, left_cons)
                    + subseq_segment_dist(right_mask, sub_mask, not left_cons)
                )
            else:
                raise ValueError("Invalid event")

    return total_cost


def label_ancestral_syntenies(
    gene_tree: PhyloTree,
    species_lca: LowestCommonAncestor,
    rec: Reconciliation,
    input_labeling: Syntenies,
) -> Tuple[ExtendedIntegral, List[Syntenies]]:
    results = MinSequence()

    if gene_tree not in input_labeling:
        prec: Dict[Any, Set[Any]] = {}

        for leaf_synteny in input_labeling.values():
            for gene_1, gene_2 in zip(leaf_synteny[0:-1], leaf_synteny[1:]):
                if gene_1 not in prec: prec[gene_1] = set()
                if gene_2 not in prec: prec[gene_2] = set()
                prec[gene_1].add(gene_2)

        for order in toposort_all(prec):
            results.update(
                _label_with_root_order(
                    gene_tree,
                    species_lca,
                    rec,
                    {
                        **input_labeling,
                        gene_tree: order
                    },
                )
            )
    else:
        results.update(
            _label_with_root_order(
                gene_tree,
                species_lca,
                rec,
                input_labeling,
            )
        )

    return (results.min, list(results))