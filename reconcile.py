#!/usr/bin/env python3
"""Test reconciliation algorithms on exhaustive supertree generation."""
from ete3 import PhyloTree, Tree
from superrec2.reconciliation.tools import (
    SuperReconciliation,
    CostType,
    get_species_name,
    reconcile_all,
    get_reconciliation_cost
)
from superrec2.reconciliation.syntenies import label_ancestral_syntenies
from superrec2.utils.trees import LowestCommonAncestor, all_supertrees
from superrec2.utils.min_sequence import MinSequence


def main():  # pylint:disable=missing-function-docstring
    gene_trees = [
        Tree("((Y_1,Z_1),(X_2,X_3));"),
        Tree("((Z_1,X_1),X_2);"),
        Tree("(X_2,X_3);"),
    ]

    species_tree = Tree("(X,(Y,Z));")
    species_lca = LowestCommonAncestor(species_tree)

    leaf_labeling = {
        "X_1": list("b"),
        "X_2": list("abc"),
        "X_3": list("ac"),
        "Y_1": list("a"),
        "Z_1": list("ab"),
    }

    costs = {
        CostType.DUPLICATION: 1,
        CostType.HORIZONTAL_GENE_TRANSFER: 1,
        CostType.FULL_LOSS: 1,
        CostType.SEGMENTAL_LOSS: 1,
    }

    results: MinSequence[SuperReconciliation] = MinSequence()

    # Give names to internal species nodes
    for node in species_tree.traverse("postorder"):
        if not node.name:
            node.name = "".join(child.name for child in node.children)

    for synteny_tree_raw in all_supertrees(gene_trees):
        synteny_tree = PhyloTree(
            synteny_tree_raw.write(format=1),
            sp_naming_function=get_species_name,
        )

        # Give names to internal synteny nodes
        next_node = 0

        for node in synteny_tree.traverse():
            if not node.name:
                node.name = str(next_node)
                next_node += 1

        for rec in reconcile_all(synteny_tree, species_lca):
            rec_cost = get_reconciliation_cost(
                synteny_tree, species_lca, rec, costs
            )

            labeling_cost, labelings = label_ancestral_syntenies(
                synteny_tree, species_lca, rec, leaf_labeling
            )

            results.update(
                (
                    rec_cost + costs[CostType.SEGMENTAL_LOSS] * labeling_cost,
                    SuperReconciliation(
                        synteny_tree=synteny_tree,
                        reconciliation=rec,
                        labeling=labelings[0],
                    ),
                )
            )

    print("Optimal cost:", results.min)
    print("Solutions:")

    for solution in results:
        print(solution)


if __name__ == "__main__":
    main()
