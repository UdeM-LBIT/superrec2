from collections import defaultdict, namedtuple
from itertools import product
from typing import List, Tuple
from ete3 import PhyloTree
from ..utils.lowest_common_ancestor import LowestCommonAncestor
from ..utils.min_sequence import MinSequence
from .tools import CostType, CostVector, ExtendedIntegral, Reconciliation


MappingInfo = namedtuple('MappingInfo', ['species', 'left', 'right'])


def reconcile_lca(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
) -> Reconciliation:
    """
    Compute the minimum-cost reconciliation between species tree and a gene
    tree, allowing for duplications and losses and where a duplication
    costs the same as loss.

    :param gene_tree: gene tree to reconcile
    :param species_tree: species tree to reconcile
    :returns: the reconciliation result
    """
    rec = {}
    species_lca = LowestCommonAncestor(species_tree)

    for gene in gene_tree.traverse("postorder"):
        if gene.is_leaf():
            species = species_tree.get_leaves_by_name(gene.species)[0]
            rec[gene] = species
        else:
            left, right = gene.children
            rec[gene] = species_lca(rec[left], rec[right])

    return rec


def _compute_thl_table(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
    costs: CostVector
) -> MinSequence[MappingInfo]:
    dup_cost = costs[CostType.Duplication]
    hgt_cost = costs[CostType.HorizontalGeneTransfer]
    loss_cost = costs[CostType.Loss]

    min_costs = defaultdict(MinSequence)
    species_lca = LowestCommonAncestor(species_tree)

    for v in gene_tree.traverse("postorder"):
        if v.is_leaf():
            s = species_tree.get_leaves_by_name(v.species)[0]
            min_costs[(v, s)].update((0, MappingInfo(
                species=s,
                left=None,
                right=None
            )))
        else:
            vl, vr = v.children

            for u in species_tree.traverse("postorder"):
                options = min_costs[(v, u)]

                # Try mapping as a speciation
                if not u.is_leaf():
                    ul, ur = u.children
                    min_vl_ul_ch = MinSequence()
                    min_vr_ul_ch = MinSequence()
                    min_vl_ur_ch = MinSequence()
                    min_vr_ur_ch = MinSequence()

                    for ul_ch in ul.traverse():
                        min_vl_ul_ch.update((min_costs[(vl, ul_ch)].min, ul_ch))
                        min_vr_ul_ch.update((min_costs[(vr, ul_ch)].min, ul_ch))

                    for ur_ch in ur.traverse():
                        min_vl_ur_ch.update((min_costs[(vl, ur_ch)].min, ur_ch))
                        min_vr_ur_ch.update((min_costs[(vr, ur_ch)].min, ur_ch))

                    # Map left gene to left species and right gene to right
                    for ul_ch, ur_ch in product(min_vl_ul_ch, min_vr_ur_ch):
                        options.update((
                            min_vl_ul_ch.min + min_vr_ur_ch.min
                            + loss_cost * (
                                species_lca.distance(u, ul_ch)
                                + species_lca.distance(u, ur_ch) - 2
                            ),
                            MappingInfo(
                                species=u,
                                left=ul_ch,
                                right=ur_ch,
                            )
                        ))

                    # Map right gene to left species and left gene to right
                    for ul_ch, ur_ch in product(min_vr_ul_ch, min_vl_ur_ch):
                        options.update((
                            min_vr_ul_ch.min + min_vl_ur_ch.min
                            + loss_cost * (
                                species_lca.distance(u, ul_ch)
                                + species_lca.distance(u, ur_ch) - 2
                            ),
                            MappingInfo(
                                species=u,
                                left=ur_ch,
                                right=ul_ch,
                            )
                        ))

                min_vl_u_ch = MinSequence()
                min_vl_u_inc = MinSequence()
                min_vr_u_ch = MinSequence()
                min_vr_u_inc = MinSequence()

                for w in species_tree.traverse():
                    if species_lca.is_ancestor_of(u, w):
                        min_vl_u_ch.update((min_costs[(vl, w)].min, w))
                        min_vr_u_ch.update((min_costs[(vr, w)].min, w))
                    elif not species_lca.is_ancestor_of(w, u):
                        min_vl_u_inc.update((min_costs[(vl, w)].min, w))
                        min_vr_u_inc.update((min_costs[(vr, w)].min, w))

                # Try mapping as a duplication
                for l_u_ch, r_u_ch in product(min_vl_u_ch, min_vr_u_ch):
                    options.update((
                        dup_cost + min_vl_u_ch.min + min_vr_u_ch.min
                        + loss_cost * (
                            species_lca.distance(u, l_u_ch)
                            + species_lca.distance(u, r_u_ch)
                        ),
                        MappingInfo(
                            species=u,
                            left=l_u_ch,
                            right=r_u_ch,
                        )
                    ))

                # Try mapping as a horizontal gene transfer
                # Map left gene as a transfer and right gene as conserved
                for l_u_inc, r_u_ch in product(min_vl_u_inc, min_vr_u_ch):
                    options.update((
                        hgt_cost + min_vl_u_inc.min + min_vr_u_ch.min
                        + loss_cost * species_lca.distance(u, r_u_ch),
                        MappingInfo(
                            species=u,
                            left=l_u_inc,
                            right=r_u_ch,
                        )
                    ))

                # Map left gene as conserved and right gene as a transfer
                for l_u_ch, r_u_inc in product(min_vl_u_ch, min_vr_u_inc):
                    options.update((
                        hgt_cost + min_vl_u_ch.min + min_vr_u_inc.min
                        + loss_cost * species_lca.distance(u, l_u_ch),
                        MappingInfo(
                            species=u,
                            left=l_u_ch,
                            right=r_u_inc,
                        )
                    ))

    return min_costs


def _decode_thl_table(root, solutions, min_costs):
    results = []

    for solution in solutions:
        if solution.left is None:
            assert solution.right is None
            results.append({root: solution.species})
        else:
            left, right = root.children
            mappings = product(
                _decode_thl_table(
                    left,
                    min_costs[(left, solution.left)],
                    min_costs,
                ),
                _decode_thl_table(
                    right,
                    min_costs[(right, solution.right)],
                    min_costs,
                ),
            )

            for map_left, map_right in mappings:
                results.append({
                    root: solution.species,
                    **map_left, **map_right,
                })

    return results


def reconcile_thl(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
    costs: CostVector
) -> Tuple[ExtendedIntegral, List[Reconciliation]]:
    """
    Find all minimum-cost reconciliations between a species tree and a gene
    tree, allowing for duplications, horizontal gene transfers and losses.

    If you wish to disallow horizontal gene transfers, using
    `reconcile_lca` is faster than setting `costs[HorizontalGeneTransfer]`
    to infinity.

    :param gene_tree: gene tree to reconcile
    :param species_tree: species tree to reconcile
    :param costs: cost of each evolutionary event
    :returns: minimum cost of each reconciliation and list of reconciliations
        with such costs
    """
    min_costs = _compute_thl_table(gene_tree, species_tree, costs or {})
    solutions = MinSequence()

    for u in species_tree.traverse("postorder"):
        for option in min_costs[(gene_tree, u)]:
            solutions.update((min_costs[(gene_tree, u)].min, option))

    return solutions.min, _decode_thl_table(gene_tree, solutions, min_costs)
