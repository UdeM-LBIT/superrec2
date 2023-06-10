"""Compute unordered synteny-labeled reconciliations."""
from collections import defaultdict
from enum import Enum, auto
from itertools import product
from typing import (
    Callable,
    Dict,
    Generator,
    Iterable,
    NamedTuple,
    Optional,
    Set,
)
from ete3 import Tree, TreeNode
from infinity import inf
from tqdm import tqdm
from .reconciliation import reconcile_lca
from ..utils.trees import LowestCommonAncestor
from ..model.synteny import GeneFamily, UnorderedSynteny, sort_synteny
from ..model.reconciliation import (
    SuperReconciliationInput,
    SuperReconciliationOutput,
    NodeEvent,
    EdgeEvent,
    CostValues,
)
from ..utils.dynamic_programming import (
    Candidate,
    DictDimension,
    Entry,
    MergePolicy,
    RetentionPolicy,
    Table,
)


class SyntenyAssignment(Enum):
    """Kinds of assignments of a synteny set to an object."""

    # Assigned to the set of all families in the object’s subtree
    LCA = auto()

    # Assigned to the set of all families in the object’s subtree,
    # and inherits extra families from the parent objects
    INHERIT = auto()


class ObjectAssignment(NamedTuple):
    """Assignment of an object to a species and synteny set."""

    species: TreeNode
    synteny: SyntenyAssignment


class ChildrenAssignment(NamedTuple):
    """Assignment of an object’s children to species and syntenies."""

    left: ObjectAssignment
    right: ObjectAssignment


class MappingChoices(NamedTuple):
    """Available mapping choices for an object."""

    # Mapping to a species of the left species subtree
    left: Entry[ObjectAssignment, int]

    # Mapping to a species of the right species subtree
    right: Entry[ObjectAssignment, int]

    # Mapping to any species in the current subtree and a conserved synteny
    conserved: Entry[ObjectAssignment, int]

    # Mapping to any species in the current subtree and a segment synteny
    segment: Entry[ObjectAssignment, int]

    # Mapping to any species in a separate subtree (i.e. neither a descendant
    # nor an ancestor of the current species) with a segment synteny
    separate: Entry[ObjectAssignment, int]


USPFSTable = Table[ChildrenAssignment, int]


def _compute_gain_sets(
    srec_input: SuperReconciliationInput,
) -> Dict[TreeNode, UnorderedSynteny]:
    """
    Given a super-reconciliation input, compute the set of gene families
    gained at each internal node of the object tree.

    :param srec_input: input super-reconciliation
    :param gains: set of families gained at each node of the object tree
    """
    leaves_by_family: Dict[GeneFamily, Set[TreeNode]] = defaultdict(set)
    object_lca = LowestCommonAncestor(srec_input.object_tree)

    for leaf, synteny in srec_input.leaf_syntenies.items():
        for family in synteny:
            leaves_by_family[family].add(leaf)

    result: Dict[TreeNode, UnorderedSynteny] = {
        node: set() for node in srec_input.object_tree.traverse()
    }

    for family, leaves in leaves_by_family.items():
        result[object_lca(*leaves)].add(family)

    return result


def _compute_lca_sets(
    srec_input: SuperReconciliationInput,
    gain_sets: Dict[TreeNode, UnorderedSynteny],
) -> Dict[TreeNode, UnorderedSynteny]:
    """
    Given a super-reconciliation input, compute the minimal set of gene
    families that must be contained in the synteny of each internal node
    of the object tree.

    :param srec_input: input super-reconciliation
    :param gains: set of families gained at each node of the object tree
    """
    result: Dict[TreeNode, UnorderedSynteny] = {}

    for object_node in srec_input.object_tree.traverse("postorder"):
        if object_node.is_leaf():
            result[object_node] = set(srec_input.leaf_syntenies[object_node])
        else:
            result[object_node] = (
                set()
                .union(*(result[child] for child in object_node.children))
                .difference(*(gain_sets[child] for child in object_node.children))
            )

    return result


def _make_event_combinator(event_cost: int):
    """Combine the assignment of both children into a single assignment."""
    return lambda left, right: Candidate(
        event_cost + left.value + right.value,
        ChildrenAssignment(left.info, right.info),
    )


def _compute_uspfs_entry(
    species_lca: LowestCommonAncestor,
    root_species: TreeNode,
    root_object: TreeNode,
    lca_sets: Dict[TreeNode, UnorderedSynteny],
    table: USPFSTable,
    costs: CostValues,
) -> None:
    """
    Compute the assignments leading to minimum-cost unordered
    super-reconciliations of the object subtree at :param:`root_object` such
    that the root object is assigned to :param:`root_species`. Results are
    stored in :param:`table`.
    """
    sloss_cost = costs[EdgeEvent.SEGMENTAL_LOSS]
    floss_cost = costs[EdgeEvent.FULL_LOSS]
    lca = SyntenyAssignment.LCA
    inh = SyntenyAssignment.INHERIT

    subprobs = tuple(
        dict(
            (
                kind,
                MappingChoices._make(
                    table.entry() for _ in range(len(MappingChoices._fields))
                ),
            )
            for kind in SyntenyAssignment
        )
        for _ in range(2)
    )

    for child_index in range(2):
        child_object = root_object.children[child_index]

        if lca_sets[root_object] <= lca_sets[child_object]:
            lca_lca_dist = 0
            lca_inh_dist = inf
        else:
            lca_lca_dist = sloss_cost
            lca_inh_dist = 0

        for desc_species in species_lca.tree.traverse():
            child_entry = table[child_object][desc_species]

            lca_cost = child_entry[lca].value()
            lca_assign = ObjectAssignment(desc_species, lca)

            inh_cost = child_entry[inh].value()
            inh_assign = ObjectAssignment(desc_species, inh)

            if species_lca.is_ancestor_of(root_species, desc_species):
                above_species_dist = (
                    species_lca.distance(root_species, desc_species) * floss_cost
                )

                subprobs[child_index][inh].conserved.update(
                    Candidate(
                        value=above_species_dist + lca_cost + sloss_cost,
                        info=lca_assign,
                    ),
                    Candidate(
                        value=above_species_dist + inh_cost,
                        info=inh_assign,
                    ),
                )
                subprobs[child_index][lca].conserved.update(
                    Candidate(
                        value=above_species_dist + lca_cost + lca_lca_dist,
                        info=lca_assign,
                    ),
                    Candidate(
                        value=above_species_dist + inh_cost + lca_inh_dist,
                        info=inh_assign,
                    ),
                )

                subprobs[child_index][inh].segment.update(
                    Candidate(
                        value=above_species_dist + lca_cost,
                        info=lca_assign,
                    ),
                    Candidate(
                        value=above_species_dist + inh_cost,
                        info=inh_assign,
                    ),
                )
                subprobs[child_index][lca].segment.update(
                    Candidate(
                        value=above_species_dist + lca_cost,
                        info=lca_assign,
                    ),
                    Candidate(
                        value=above_species_dist + inh_cost + lca_inh_dist,
                        info=inh_assign,
                    ),
                )

                if not root_species.is_leaf():
                    left_species, right_species = root_species.children
                    species_dist = above_species_dist - floss_cost
                    inh_candidates = (
                        Candidate(
                            value=species_dist + lca_cost + sloss_cost,
                            info=lca_assign,
                        ),
                        Candidate(
                            value=species_dist + inh_cost,
                            info=inh_assign,
                        ),
                    )
                    lca_candidates = (
                        Candidate(
                            value=species_dist + lca_cost + lca_lca_dist,
                            info=lca_assign,
                        ),
                        Candidate(
                            value=species_dist + inh_cost + lca_inh_dist,
                            info=inh_assign,
                        ),
                    )

                    if species_lca.is_ancestor_of(left_species, desc_species):
                        subprobs[child_index][inh].left.update(*inh_candidates)
                        subprobs[child_index][lca].left.update(*lca_candidates)
                    elif species_lca.is_ancestor_of(right_species, desc_species):
                        subprobs[child_index][inh].right.update(*inh_candidates)
                        subprobs[child_index][lca].right.update(*lca_candidates)
            elif not species_lca.is_ancestor_of(desc_species, root_species):
                subprobs[child_index][inh].separate.update(
                    Candidate(value=lca_cost, info=lca_assign),
                    Candidate(value=inh_cost, info=inh_assign),
                )
                subprobs[child_index][lca].separate.update(
                    Candidate(value=lca_cost, info=lca_assign),
                    Candidate(value=inh_cost + lca_inh_dist, info=inh_assign),
                )

    spe_comb = _make_event_combinator(costs[NodeEvent.SPECIATION])
    dup_comb = _make_event_combinator(costs[NodeEvent.DUPLICATION])
    hgt_comb = _make_event_combinator(costs[NodeEvent.HORIZONTAL_TRANSFER])

    for kind in SyntenyAssignment:
        table[root_object][root_species][kind].update(
            *subprobs[0][kind].left.combine(subprobs[1][kind].right, spe_comb),
            *subprobs[0][kind].right.combine(subprobs[1][kind].left, spe_comb),
            *subprobs[0][kind].conserved.combine(subprobs[1][kind].segment, dup_comb),
            *subprobs[0][kind].segment.combine(subprobs[1][kind].conserved, dup_comb),
            *subprobs[0][kind].conserved.combine(subprobs[1][kind].separate, hgt_comb),
            *subprobs[0][kind].separate.combine(subprobs[1][kind].conserved, hgt_comb),
        )


def _compute_uspfs_table(
    srec_input: SuperReconciliationInput,
    lca_sets: Dict[TreeNode, UnorderedSynteny],
    allowed_species: Callable[[Tree, TreeNode], Iterable[TreeNode]],
    retention_policy: RetentionPolicy,
) -> USPFSTable:
    """
    Compute an assignment table that can be used to construct minimum-cost
    unordered super-reconciliations. Each entry `table[object][species][kind]`
    of the resulting table contains the possible assignments of `object`’s
    children that lead to minimum-cost unordered super-reconciliations of the
    subtrees rooted at `object` such that the root is mapped to `species` and
    to the LCA synteny or to a larger synteny depending on the value of `kind`.

    :param srec_input: objects of the unordered super-reconciliation
    :param lca_sets: set of genes in the LCA synteny of each object
    :param allowed_species: callable that gives the set of allowed species for
        each object
    :param retention_policy: whether to keep any minimal assignment in each
        entry or all minimal assignments
    :returns: computed assignment table
    """
    table = Table(
        (DictDimension(), DictDimension(), DictDimension()),
        MergePolicy.MIN,
        retention_policy,
    )
    lca = SyntenyAssignment.LCA

    for root_object in tqdm(
        srec_input.object_tree.traverse("postorder"),
        desc="Table entries",
        total=sum(1 for _ in srec_input.object_tree.traverse()),
        ascii=True,
        leave=False,
    ):
        if root_object.is_leaf():
            species = srec_input.leaf_object_species[root_object]
            table[root_object][species][lca] = Candidate(0)
        else:
            for root_species in tqdm(
                allowed_species(srec_input.species_lca.tree, root_object),
                desc="Object assignments",
                total=sum(
                    1 for _ in allowed_species(srec_input.species_lca.tree, root_object)
                ),
                ascii=True,
                leave=False,
            ):
                _compute_uspfs_entry(
                    srec_input.species_lca,
                    root_species,
                    root_object,
                    lca_sets,
                    table,
                    srec_input.costs,
                )

    return table


def _decode_uspfs_table(
    root_object: TreeNode,
    root_species: TreeNode,
    root_kind: SyntenyAssignment,
    ancestor_synteny: Optional[UnorderedSynteny],
    srec_input: SuperReconciliationInput,
    gain_sets: Dict[TreeNode, UnorderedSynteny],
    lca_sets: Dict[TreeNode, UnorderedSynteny],
    table: USPFSTable,
) -> Generator[SuperReconciliationOutput, None, None]:
    """
    Recursively reconstruct minimum-cost unordered super-reconciliations from a
    pre-computed assignment table.

    :param root_ordering: complete sequence of families of which other
        syntenies are subsequences
    :param root_object: root of the current object subtree
    :param root_species: species to which the current root is mapped
    :param root_kind: kind of synteny to which the current root is mapped
    :param ancestor_synteny: synteny to inherit from if the kind is INHERIT
    :param srec_input: objects of the unordered super-reconciliation
    :param table: table of optimal assignments as computed by
        :func:`_compute_uspfs_table`
    :returns: yields minimum-cost unordered super-reconciliations
    """
    if root_kind == SyntenyAssignment.LCA:
        ancestor_synteny = lca_sets[root_object]
        root_synteny = sort_synteny(lca_sets[root_object])
    else:
        ancestor_synteny |= gain_sets[root_object]
        root_synteny = sort_synteny(ancestor_synteny)

    if (
        root_object.is_leaf()
        and not table[root_object][root_species][root_kind].is_infinite()
    ):
        # pylint has trouble seeing the attributes from the descendant dataclass
        yield SuperReconciliationOutput(  # pylint: disable=unexpected-keyword-arg
            input=srec_input,
            object_species={root_object: root_species},
            syntenies={root_object: root_synteny},
            ordered=False,
        )
        return

    for info in table[root_object][root_species][root_kind].infos():
        left_object, right_object = root_object.children
        mappings = product(
            _decode_uspfs_table(
                left_object,
                info.left.species,
                info.left.synteny,
                ancestor_synteny,
                srec_input,
                gain_sets,
                lca_sets,
                table,
            ),
            _decode_uspfs_table(
                right_object,
                info.right.species,
                info.right.synteny,
                ancestor_synteny,
                srec_input,
                gain_sets,
                lca_sets,
                table,
            ),
        )

        for map_left, map_right in mappings:
            # pylint has trouble seeing the attributes from the descendant dataclass
            yield SuperReconciliationOutput(  # pylint: disable=unexpected-keyword-arg
                input=srec_input,
                object_species={
                    root_object: root_species,
                    **map_left.object_species,
                    **map_right.object_species,
                },
                syntenies={
                    root_object: root_synteny,
                    **map_left.syntenies,
                    **map_right.syntenies,
                },
                ordered=False,
            )


def _uspfs(
    srec_input: SuperReconciliationInput,
    policy: RetentionPolicy,
    allowed_species: Callable[[TreeNode], Iterable[TreeNode]],
) -> Set[SuperReconciliationOutput]:
    results: Entry[int, SuperReconciliationOutput] = Entry(MergePolicy.MIN, policy)

    for srec_input_bin in tqdm(
        srec_input.binarize(),
        desc="Multifurcation resolutions",
        total=sum(1 for _ in srec_input.binarize()),
        ascii=True,
    ):
        srec_input_bin.label_internal()
        gain_sets = _compute_gain_sets(srec_input_bin)
        lca_sets = _compute_lca_sets(srec_input_bin, gain_sets)
        synteny_tree = srec_input_bin.object_tree

        table = _compute_uspfs_table(
            srec_input_bin,
            lca_sets,
            allowed_species,
            policy,
        )

        for root_species in tqdm(
            srec_input_bin.species_lca.tree.traverse(),
            desc="Generate solutions",
            total=sum(1 for _ in srec_input.species_lca.tree.traverse()),
            ascii=True,
            leave=False,
        ):
            results.update(
                *map(
                    lambda output: Candidate(output.cost(), output),
                    _decode_uspfs_table(
                        synteny_tree,
                        root_species,
                        SyntenyAssignment.LCA,
                        lca_sets[synteny_tree],
                        srec_input_bin,
                        gain_sets,
                        lca_sets,
                        table,
                    ),
                )
            )

    return results.infos()


def usreconcile_base_uspfs(
    srec_input: SuperReconciliationInput,
    policy: RetentionPolicy,
) -> Set[SuperReconciliationOutput]:
    """
    Compute a minimum-cost unordered super-reconciliation using the original
    Unordered Small-Phylogeny-for-Syntenies (USPFS) algorithm. This cannot
    generate horizontal gene transfer events and assumes that the cost of
    duplications is equal to the cost of losses.

    :param srec_input: objects of the super-reconciliation
    :param policy: whether to generate any minimal solution or all possible
        minimal solutions
    :returns: any minimum-cost super-reconciliation, if there is one
    """
    rec_output = reconcile_lca(srec_input)
    return _uspfs(
        srec_input,
        policy,
        allowed_species=lambda _, obj: [rec_output.object_species[obj]],
    )


def usreconcile_extended_uspfs(
    srec_input: SuperReconciliationInput,
    policy: RetentionPolicy,
) -> Set[SuperReconciliationOutput]:
    """
    Compute a minimum-cost unordered super-reconciliation using the
    SuperDTL algorithm. This is an extension of the Unordered
    Small-Phylogeny-for-Syntenies (EUSPFS) algorithm that can handle
    all kinds of events (including transfers) and cost values.

    :param srec_input: objects of the super-reconciliation
    :param policy: whether to generate any minimal solution or all possible
        minimal solutions
    :returns: any minimum-cost super-reconciliation, if there is one
    """
    return _uspfs(
        srec_input,
        policy,
        allowed_species=lambda species, _: species.traverse("postorder"),
    )
