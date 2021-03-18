#!/usr/bin/env python3
from itertools import product
from string import ascii_lowercase
from infinity import inf
from ete3 import PhyloNode
from superrec2 import reconciliation as rec
from superrec2.utils.lowest_common_ancestor import LowestCommonAncestor


def trees_all(labels):
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
            root.add_child(left)
            root.add_child(right)
            root.name = left.name + right.name
            yield root


def map_all(domain, values, next_value=0):
    if not domain:
        yield {}
        return

    for value in range(next_value):
        for rest in map_all(domain[1:], values, next_value):
            yield {domain[0]: value, **rest}

    if next_value < values:
        for rest in map_all(domain[1:], values, next_value + 1):
            yield {domain[0]: next_value, **rest}


def subseq_all(sequence):
    if len(sequence) == 0:
        yield []
        return

    for subsequence in subseq_all(sequence[1:]):
        yield subsequence
        yield [sequence[0]] + subsequence


def inputs_all(species_count, synteny_count, synteny_size):
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
                    synteny_tree.set_species_naming_function(
                        rec.tools.get_species_name
                    )
                    leaf_labeling = {
                        synteny_tree & synteny_names[i]: synteny_value
                        for i, synteny_value in enumerate(leaf_labeling_tuple)
                    }

                    yield (
                        species_tree,
                        species_lca,
                        synteny_tree,
                        leaf_labeling,
                    )


def main():  # pylint:disable=missing-function-docstring
    costs = {
        rec.tools.CostType.DUPLICATION: 1,
        rec.tools.CostType.HORIZONTAL_GENE_TRANSFER: 1,
        rec.tools.CostType.LOSS: 1,
    }

    for species_tree, species_lca, synteny_tree, leaf_labeling in inputs_all(
        species_count=2,
        synteny_count=4,
        synteny_size=2,
    ):
        print("Input:\n")
        print(species_tree.write(format=8, format_root_node=True))
        print(synteny_tree.write(format=8, format_root_node=True))

        for synteny_node, synteny_value in leaf_labeling.items():
            print(f"{synteny_node.name}:{''.join(synteny_value)}", end=",")

        print("\n\nOutput:\n")

        thl_cost, thl_results = rec.compute.reconcile_thl(
            synteny_tree, species_lca, costs
        )

        thl_label_cost = inf

        print("THL cost:", thl_cost)

        for thl_result in thl_results:
            result_cost = rec.tools.get_reconciliation_cost(
                synteny_tree, species_lca, thl_result, costs
            )
            assert result_cost == thl_cost

            label_cost, labelings = rec.syntenies.label_ancestral_syntenies(
                synteny_tree, species_lca, thl_result, leaf_labeling
            )

            for labeling in labelings:
                local_label_cost = rec.tools.get_labeling_cost(
                    synteny_tree, species_lca, thl_result, labeling
                )
                assert local_label_cost == label_cost

            thl_label_cost = min(thl_label_cost, result_cost + label_cost)

            for synteny_node, species_node in thl_result.items():
                print(f"{synteny_node.name}:{species_node.name}", end=",")

            print()

            for synteny_node, synteny_value in labelings[0].items():
                print(f"{synteny_node.name}:{''.join(synteny_value)}", end=",")

            print()

        print("\nLabeled THL cost:", thl_label_cost)
        all_label_cost = inf
        all_mini_rec = None
        all_mini_label = None

        for possible_rec in rec.tools.reconcile_all(
            synteny_tree, species_tree, species_lca
        ):
            rec_cost = rec.tools.get_reconciliation_cost(
                synteny_tree, species_lca, possible_rec, costs
            )
            assert rec_cost >= thl_cost

            label_cost, labelings = rec.syntenies.label_ancestral_syntenies(
                synteny_tree, species_lca, possible_rec, leaf_labeling
            )

            for labeling in labelings:
                local_label_cost = rec.tools.get_labeling_cost(
                    synteny_tree, species_lca, possible_rec, labeling
                )
                assert local_label_cost == label_cost

            if all_label_cost > rec_cost + label_cost:
                all_label_cost = rec_cost + label_cost
                all_mini_rec = possible_rec
                all_mini_label = labelings[0]

        if all_mini_rec:
            for synteny_node, species_node in all_mini_rec.items():
                print(f"{synteny_node.name}:{species_node.name}", end=",")
            print()

        if all_mini_label:
            for synteny_node, synteny_value in all_mini_label.items():
                print(f"{synteny_node.name}:{''.join(synteny_value)}", end=",")
            print()

        print("Optimal cost:", all_label_cost, "\n")
        assert thl_label_cost == all_label_cost

if __name__ == "__main__":
    main()
