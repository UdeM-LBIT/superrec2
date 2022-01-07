"""Compute and represent synteny-labeled reconciliations."""
from itertools import product
from typing import (
    Any,
    Callable,
    Dict,
    Generator,
    Iterable,
    Optional,
    NamedTuple,
    Sequence,
    Set,
)
from ete3 import Tree, TreeNode
from .reconciliation import reconcile_lca
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
    Candidate,
    DictDimension,
    Entry,
    MergePolicy,
    RetentionPolicy,
    Table,
)


class ObjectAssignment(NamedTuple):
    """Assignment of an object to a species and synteny."""

    species: TreeNode
    synteny: int


class ChildrenAssignment(NamedTuple):
    """Assignment of an object’s children to species and syntenies."""

    left: ObjectAssignment
    right: ObjectAssignment


SPFSTable = Table[ChildrenAssignment, int]


class SpeciesMappingChoices(NamedTuple):
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


def _make_event_combinator(event_cost: int):
    """Combine the assignment of both children into a single assignment."""
    return lambda left, right: Candidate(
        event_cost + left.value + right.value,
        ChildrenAssignment(left.info, right.info),
    )


def _compute_spfs_entry(
    species_lca: LowestCommonAncestor,
    root_species: TreeNode,
    root_synteny: int,
    root_object: TreeNode,
    table: SPFSTable,
    costs: CostValues,
) -> None:
    """
    Compute the assignments leading to minimum-cost super-reconciliations of
    the object subtree at :param:`root_object` such that the root object is
    assigned to :param:`root_species` and :param:`root_synteny`. Results are
    stored in :param:`table`.
    """
    sloss_cost = costs[EdgeEvent.SEGMENTAL_LOSS]
    floss_cost = costs[EdgeEvent.FULL_LOSS]

    subprobs = tuple(
        SpeciesMappingChoices._make(
            table.entry() for _ in range(len(SpeciesMappingChoices._fields))
        )
        for _ in range(2)
    )

    for child_index in range(2):
        child_object = root_object.children[child_index]

        for desc_species in species_lca.tree.traverse():
            for child_synteny in table[child_object][desc_species]:
                conserv_dist = (
                    subseq_segment_dist(
                        child_synteny,
                        root_synteny,
                        edges=True,
                    )
                    * sloss_cost
                )

                if conserv_dist < 0:
                    # Not a subsequence of the parent synteny
                    continue

                segment_dist = (
                    subseq_segment_dist(
                        child_synteny, root_synteny, edges=False
                    )
                    * sloss_cost
                )

                sub_cost = table[child_object][desc_species][
                    child_synteny
                ].value()
                assignment = ObjectAssignment(desc_species, child_synteny)

                if species_lca.is_ancestor_of(root_species, desc_species):
                    above_species_dist = (
                        species_lca.distance(root_species, desc_species)
                        * floss_cost
                    )

                    subprobs[child_index].conserved.update(
                        Candidate(
                            value=above_species_dist + sub_cost + conserv_dist,
                            info=assignment,
                        )
                    )
                    subprobs[child_index].segment.update(
                        Candidate(
                            value=above_species_dist + sub_cost + segment_dist,
                            info=assignment,
                        )
                    )

                    if not root_species.is_leaf():
                        left_species, right_species = root_species.children
                        is_left_desc = species_lca.is_ancestor_of(
                            left_species, desc_species
                        )
                        is_right_desc = species_lca.is_ancestor_of(
                            right_species, desc_species
                        )
                        species_dist = above_species_dist - floss_cost

                        if is_left_desc:
                            subprobs[child_index].left.update(
                                Candidate(
                                    value=species_dist
                                    + sub_cost
                                    + conserv_dist,
                                    info=assignment,
                                )
                            )
                        elif is_right_desc:
                            subprobs[child_index].right.update(
                                Candidate(
                                    value=species_dist
                                    + sub_cost
                                    + conserv_dist,
                                    info=assignment,
                                )
                            )
                elif not species_lca.is_ancestor_of(desc_species, root_species):
                    subprobs[child_index].separate.update(
                        Candidate(
                            value=sub_cost + segment_dist,
                            info=assignment,
                        )
                    )

    spe_comb = _make_event_combinator(costs[NodeEvent.SPECIATION])
    dup_comb = _make_event_combinator(costs[NodeEvent.DUPLICATION])
    hgt_comb = _make_event_combinator(costs[NodeEvent.HORIZONTAL_TRANSFER])

    table[root_object][root_species][root_synteny].update(
        *subprobs[0].left.combine(subprobs[1].right, spe_comb),
        *subprobs[0].right.combine(subprobs[1].left, spe_comb),
        *subprobs[0].conserved.combine(subprobs[1].segment, dup_comb),
        *subprobs[0].segment.combine(subprobs[1].conserved, dup_comb),
        *subprobs[0].conserved.combine(subprobs[1].separate, hgt_comb),
        *subprobs[0].separate.combine(subprobs[1].conserved, hgt_comb),
    )


def _compute_spfs_table(
    srec_input: SuperReconciliationInput,
    root_ordering: Synteny,
    allowed_species: Callable[[Tree, TreeNode], Iterable[TreeNode]],
    allowed_syntenies: Callable[[Synteny, TreeNode], Iterable[int]],
    retention_policy: RetentionPolicy,
) -> SPFSTable:
    """
    Compute an assignment table that can be used to construct minimum-cost
    super-reconciliations. Each entry `table[object][species][synteny]` of the
    resulting table contains the possible assignments of `object`’s children
    that lead to minimum-cost super-reconciliations of the subtrees rooted at
    `object` such that the root is mapped to `species` and `synteny`.

    :param srec_input: objects of the super-reconciliation
    :param root_ordering: complete sequence of families of which other
        syntenies are subsequences
    :param allowed_species: callable that gives the set of allowed species for
        each object
    :param allowed_syntenies: callable that gives the set of allowed syntenies
        for each object (specified as a binary mask of the complete sequence)
    :param retention_policy: whether to keep any minimal assignment in each
        entry or all minimal assignments
    :returns: computed assignment table
    """
    table = Table(
        (DictDimension(), DictDimension(), DictDimension()),
        MergePolicy.MIN,
        retention_policy,
    )

    for root_object in srec_input.object_tree.traverse("postorder"):
        if root_object.is_leaf():
            synteny = mask_from_subseq(
                srec_input.leaf_syntenies[root_object], root_ordering
            )
            species = srec_input.leaf_object_species[root_object]
            table[root_object][species][synteny] = Candidate(0)
        else:
            # Test all possible species mapping and synteny subsequences for
            # the current node’s synteny
            for root_species, root_synteny in product(
                allowed_species(srec_input.species_lca.tree, root_object),
                allowed_syntenies(root_ordering, root_object),
            ):
                _compute_spfs_entry(
                    srec_input.species_lca,
                    root_species,
                    root_synteny,
                    root_object,
                    table,
                    srec_input.costs,
                )

    return table


def _decode_spfs_table(
    root_ordering: Sequence[Any],
    root_object: TreeNode,
    root_species: TreeNode,
    root_synteny: int,
    srec_input: SuperReconciliationInput,
    table: SPFSTable,
) -> Generator[SuperReconciliationOutput, None, None]:
    """
    Recursively reconstruct minimum-cost super-reconciliations from a
    pre-computed assignment table.

    :param root_ordering: complete sequence of families of which other
        syntenies are subsequences
    :param root_object: root of the current object subtree
    :param root_species: species to which the current root is mapped
    :param root_synteny: synteny to which the current root is mapped
    :param srec_input: objects of the super-reconciliation
    :param table: table of optimal assignments as computed by
        :func:`_compute_spfs_table`
    :returns: yields minimum-cost super-reconciliations
    """
    resolv_synteny = subseq_from_mask(root_synteny, root_ordering)

    if (
        root_object.is_leaf()
        and not table[root_object][root_species][root_synteny].is_infinite()
    ):
        yield SuperReconciliationOutput(
            input=srec_input,
            object_species={root_object: root_species},
            syntenies={root_object: resolv_synteny},
        )
        return

    for info in table[root_object][root_species][root_synteny].infos():
        left_object, right_object = root_object.children
        mappings = product(
            _decode_spfs_table(
                root_ordering,
                left_object,
                info.left.species,
                info.left.synteny,
                srec_input,
                table,
            ),
            _decode_spfs_table(
                root_ordering,
                right_object,
                info.right.species,
                info.right.synteny,
                srec_input,
                table,
            ),
        )

        for map_left, map_right in mappings:
            yield SuperReconciliationOutput(
                input=srec_input,
                object_species={
                    root_object: root_species,
                    **map_left.object_species,
                    **map_right.object_species,
                },
                syntenies={
                    root_object: resolv_synteny,
                    **map_left.syntenies,
                    **map_right.syntenies,
                },
            )


def _make_prec_graph(leaf_syntenies: SyntenyMapping):
    """
    Construct a precedence graph of gene families encoding the relative
    ordering of families in the given mapping.
    """
    prec: Dict[Any, Set[Any]] = {}

    for leaf_synteny in leaf_syntenies.values():
        for gene_1, gene_2 in zip(leaf_synteny[0:-1], leaf_synteny[1:]):
            if gene_1 not in prec:
                prec[gene_1] = set()
            prec[gene_1].add(gene_2)

        if leaf_synteny[-1] not in prec:
            prec[leaf_synteny[-1]] = set()

    return prec


def _spfs(
    srec_input: SuperReconciliationInput,
    policy: RetentionPolicy,
    allowed_species: Callable[[TreeNode], Iterable[TreeNode]],
    allowed_syntenies: Callable[[TreeNode], Iterable[int]],
) -> Set[SuperReconciliationOutput]:
    rec_output = reconcile_lca(srec_input)
    synteny_tree = srec_input.object_tree
    leaf_syntenies = srec_input.leaf_syntenies

    if synteny_tree not in leaf_syntenies:
        root_orderings = toposort_all(_make_prec_graph(leaf_syntenies))
    else:
        root_orderings = (leaf_syntenies[synteny_tree],)

    results: Entry[int, SuperReconciliationOutput] = Entry(
        MergePolicy.MIN, policy
    )

    for root_species, root_ordering in product(
        srec_input.species_lca.tree.traverse(), root_orderings
    ):
        table = _compute_spfs_table(
            srec_input,
            root_ordering,
            allowed_species,
            allowed_syntenies,
            policy,
        )

        results.update(
            *map(
                lambda output: Candidate(output.cost(), output),
                _decode_spfs_table(
                    root_ordering,
                    synteny_tree,
                    root_species,
                    subseq_complete(root_ordering),
                    srec_input,
                    table,
                ),
            )
        )

    return results.infos()


def sreconcile_base_spfs(
    srec_input: SuperReconciliationInput,
    policy: RetentionPolicy,
) -> Set[SuperReconciliationOutput]:
    """
    Compute a minimum-cost super-reconciliation using the original
    Small-Phylogeny-for-Syntenies (SPFS) algorithm. This cannot generate
    horizontal gene transfer events and assumes that the cost of duplications
    is equal to the cost of losses. If the input does not include a root
    synteny, this algorithm will consider all possible orderings for the root
    (this can take a lot of time!).

    :param srec_input: objects of the super-reconciliation
    :param policy: whether to generate any minimal solution or all possible
        minimal solutions
    :returns: any minimum-cost super-reconciliation, if there is one
    """
    rec_output = reconcile_lca(srec_input)
    return _spfs(
        srec_input,
        policy,
        allowed_species=lambda _, obj: [rec_output.object_species[obj]],
        allowed_syntenies=lambda ordering, obj: (
            (subseq_complete(ordering),)
            if obj == srec_input.object_tree
            else range(2 ** len(ordering))
        ),
    )


def sreconcile_extended_spfs(
    srec_input: SuperReconciliationInput,
    policy: RetentionPolicy,
) -> Set[SuperReconciliationOutput]:
    """
    Compute a minimum-cost super-reconciliation using the Extended
    Small-Phylogeny-for-Syntenies (ESPFS) algorithm. This can handle all kinds
    of events and costs. If the input does not include a root synteny, this
    algorithm will consider all possible orderings for the root (this can take
    a lot of time!).

    :param srec_input: objects of the super-reconciliation
    :param policy: whether to generate any minimal solution or all possible
        minimal solutions
    :returns: any minimum-cost super-reconciliation, if there is one
    """
    return _spfs(
        srec_input,
        policy,
        allowed_species=lambda species, _: species.traverse("postorder"),
        allowed_syntenies=lambda ordering, obj: (
            (subseq_complete(ordering),)
            if obj == srec_input.object_tree
            else range(2 ** len(ordering))
        ),
    )