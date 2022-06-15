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
    """Synteny sets that can be assigned to an object."""

    # Minimal set of families that must be present in the object’s synteny
    MINIMAL = auto()

    # Any set of families strictly larger set than the minimal set
    INHERIT = auto()


class ObjectAssignment(NamedTuple):
    """Assignment of an object to a species and synteny set."""

    species: TreeNode
    synteny: SyntenyAssignment


class ChildrenAssignment(NamedTuple):
    """Assignment of an object’s children to species and syntenies."""

    left: ObjectAssignment
    right: ObjectAssignment


AssignmentTable = Table[ChildrenAssignment, int]
SpeciesTable = Table[TreeNode, int]


class SuperReconciliationTables(NamedTuple):
    """Dynamic programming tables used in the super-reconciliation algorithm."""

    # `all[object][species][kind]` contains the minimum-cost assignment of
    # `object`’s children such that `object` is assigned to `species` and to
    # a synteny according to `kind`
    all: AssignmentTable

    # `inside[object][species][kind]` contains the minimum-cost assignment
    # of `object` to species `s` descending from `species`, including the cost
    # for going from `species` to `s`
    inside: SpeciesTable

    # `insideSkip[object][species][kind]` is similar to `inside` but does not
    # include the cost for going from `species` to `s`
    insideSkip: SpeciesTable

    # `outside[object][species][kind]` contains the minimum-cost assignment
    # of `object` to species `s` separate from `species`
    outside: SpeciesTable

    @classmethod
    def make_empty(cls, retention_policy: RetentionPolicy):
        tables = {}

        for table in ("all", "inside", "insideSkip", "outside"):
            tables[table] = Table(
                (DictDimension(), DictDimension(), DictDimension()),
                MergePolicy.MIN,
                retention_policy,
            )

        return cls(**tables)


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


def _compute_minimal_sets(
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
                .difference(
                    *(gain_sets[child] for child in object_node.children)
                )
            )

    return result


def _initialize_tables(
    srec_input: SuperReconciliationInput,
    retention_policy: RetentionPolicy,
) -> SuperReconciliationTables:
    """Initialize dynamic programming tables for object leaves."""
    tables = SuperReconciliationTables.make_empty(retention_policy)

    minimal = SyntenyAssignment.MINIMAL
    floss_cost = srec_input.costs[EdgeEvent.FULL_LOSS]
    species_lca = srec_input.species_lca

    for obj in srec_input.object_tree:
        leaf = srec_input.leaf_object_species[obj]
        tables.all[obj][leaf][minimal] = Candidate(0)

        for target in species_lca.tree.traverse("postorder"):
            if species_lca.is_ancestor_of(target, leaf):
                tables.inside[obj][target][minimal] = Candidate(
                    (
                        species_lca.distance(target, leaf)
                        * floss_cost
                    ),
                    leaf
                )
                tables.insideSkip[obj][target][minimal] = Candidate(0, leaf)
            else:
                tables.outside[obj][target][minimal] = Candidate(0, leaf)

    return tables


def _compute_segmental_cost(
    minimal_sets,
    orig_object,
    child_object,
    kind,
    kind_child,
    partial
) -> Optional[int]:
    if kind == SyntenyAssignment.MINIMAL and \
            minimal_sets[orig_object] <= minimal_sets[child_object]:
        if kind_child == SyntenyAssignment.MINIMAL:
            # Both parent and child have the minimal synteny, but all
            # minimal families of the parent are present in the child
            return 0

        # Parent has the minimal synteny while child has a larger synteny,
        # but all minimal families of the parent are present in the child,
        # so the child cannot actually have a larger synteny!
        return None

    if kind_child == SyntenyAssignment.MINIMAL:
        # Child has the minimal synteny while either:
        # (1) Parent has a larger synteny
        # (2) Parent has the minimal synteny and not all minimal families of
        #     the parent are present in the child
        return 1 if not partial else 0

    # Child has a larger synteny while either:
    # (1) Parent also has a larger synteny
    # (2) Parent has the minimal synteny and not all minimal families of
    #     the parent are present in the child
    return 0


def _event_cost(base_cost, count):
    if count == 0:
        return 0

    return base_cost * count


def _compute_tables_entry(
    species_lca: LowestCommonAncestor,
    orig_object: TreeNode,
    minimal_sets: Dict[TreeNode, UnorderedSynteny],
    tables: SuperReconciliationTables,
    costs: CostValues,
) -> None:
    """
    Compute the assignments leading to minimum-cost unordered
    super-reconciliations of the object subtree at :param:`orig_object`.
    """
    floss_cost = costs[EdgeEvent.FULL_LOSS]
    sloss_cost = costs[EdgeEvent.SEGMENTAL_LOSS]
    tloss_cost = costs[EdgeEvent.TRANSFER_LOSS]
    xgt_cost = costs[EdgeEvent.UNSAMPLED_TRANSFER]

    spe_cost = costs[NodeEvent.SPECIATION]
    dup_cost = costs[NodeEvent.DUPLICATION]
    hgt_cost = costs[NodeEvent.HORIZONTAL_TRANSFER]

    minimal = SyntenyAssignment.MINIMAL
    inherit = SyntenyAssignment.INHERIT

    left, right = orig_object.children

    # Compute the optimal cost of assigning `orig_object` to any species
    for (
        orig_species,
        (kind, kind_left, kind_right),
        (table_name_left, table_name_right)
    ) in product(
        species_lca.tree.traverse("postorder"),
        product(SyntenyAssignment, repeat=3),
        product(("inside", "insideSkip", "outside"), repeat=2),
    ):
        table_left = getattr(tables, table_name_left)
        table_right = getattr(tables, table_name_right)

        xgt_left = 1 if table_name_left == "insideSkip" else 0
        xgt_right = 1 if table_name_right == "insideSkip" else 0

        tloss_left = 1 if table_name_left == "outside" else 0
        tloss_right = 1 if table_name_right == "outside" else 0

        partial_left = 1 if table_name_left != "inside" else 0
        partial_right = 1 if table_name_right != "inside" else 0

        segmental_left = _compute_segmental_cost(
            minimal_sets,
            orig_object,
            left,
            kind,
            kind_left,
            partial_left,
        )

        segmental_right = _compute_segmental_cost(
            minimal_sets,
            orig_object,
            right,
            kind,
            kind_right,
            partial_right,
        )

        if segmental_left is None or segmental_right is None:
            continue

        # Try to map as a speciation event
        if not orig_species.is_leaf():
            spe_left, spe_right = orig_species.children
            bounds = ((spe_left, spe_right), (spe_right, spe_left))

            for (bound_left, bound_right) in bounds:
                tables.all[orig_object][orig_species][kind].update(
                    *table_left[left][bound_left][kind_left].combine(
                        table_right[right][bound_right][kind_right],
                        lambda target_left, target_right: Candidate(
                            value=(
                                spe_cost
                                + target_left.value + target_right.value
                                + _event_cost(sloss_cost, segmental_left + segmental_right)
                                + _event_cost(xgt_cost, xgt_left + xgt_right)
                                + _event_cost(tloss_cost, tloss_left + tloss_right)
                            ),
                            info=ChildrenAssignment(
                                left=ObjectAssignment(
                                    species=target_left.info,
                                    synteny=kind_left,
                                ),
                                right=ObjectAssignment(
                                    species=target_right.info,
                                    synteny=kind_right
                                ),
                            )
                        )
                    )
                )

        # Try to map as a duplication event
        tables.all[orig_object][orig_species][kind].update(
            *table_left[left][orig_species][kind_left].combine(
                table_right[right][orig_species][kind_right],
                lambda target_left, target_right: Candidate(
                    value=(
                        dup_cost
                        + target_left.value + target_right.value
                        + _event_cost(sloss_cost, min(1, segmental_left + segmental_right))
                        + _event_cost(xgt_cost, xgt_left + xgt_right)
                        + _event_cost(tloss_cost, tloss_left + tloss_right)
                    ),
                    info=ChildrenAssignment(
                        left=ObjectAssignment(
                            species=target_left.info,
                            synteny=kind_left,
                        ),
                        right=ObjectAssignment(
                            species=target_right.info,
                            synteny=kind_right
                        ),
                    )
                )
            )
        )

        # Try to map as a transfer event (transfer from the right)
        tables.all[orig_object][orig_species][kind].update(
            *table_left[left][orig_species][kind_left].combine(
                tables.outside[right][orig_species][kind_right],
                lambda target_left, target_right: Candidate(
                    value=(
                        hgt_cost
                        + target_left.value + target_right.value
                        + _event_cost(sloss_cost, segmental_left)
                        + _event_cost(xgt_cost, xgt_left)
                        + _event_cost(tloss_cost, tloss_left)
                    ),
                    info=ChildrenAssignment(
                        left=ObjectAssignment(
                            species=target_left.info,
                            synteny=kind_left,
                        ),
                        right=ObjectAssignment(
                            species=target_right.info,
                            synteny=kind_right
                        ),
                    )
                )
            )
        )

        # (transfer from the left)
        tables.all[orig_object][orig_species][kind].update(
            *tables.outside[left][orig_species][kind_left].combine(
                table_right[right][orig_species][kind_right],
                lambda target_left, target_right: Candidate(
                    value=(
                        hgt_cost
                        + target_left.value + target_right.value
                        + _event_cost(sloss_cost, segmental_right)
                        + _event_cost(xgt_cost, xgt_right)
                        + _event_cost(tloss_cost, tloss_right)
                    ),
                    info=ChildrenAssignment(
                        left=ObjectAssignment(
                            species=target_left.info,
                            synteny=kind_left,
                        ),
                        right=ObjectAssignment(
                            species=target_right.info,
                            synteny=kind_right
                        ),
                    )
                )
            )
        )

    # Compute auxiliary tables `inside` and `insideSkip` based on the
    # newly-computed optimal assignment of `orig_object` in `all`
    for target in species_lca.tree.traverse("postorder"):
        for kind in SyntenyAssignment:
            tables.inside[orig_object][target][kind] = Candidate(
                tables.all[orig_object][target][kind].value(),
                target,
            )

            tables.insideSkip[orig_object][target][kind] = Candidate(
                tables.all[orig_object][target][kind].value(),
                target,
            )

            for child in target.children:
                in_val = tables.inside[orig_object][child][kind].value()
                in_infos = tables.inside[orig_object][child][kind].infos()

                for target_child in in_infos:
                    tables.inside[orig_object][target][kind] = Candidate(
                        in_val + floss_cost,
                        target_child,
                    )

                tables.insideSkip[orig_object][target][kind].update(
                    *tables.insideSkip[orig_object][child][kind]
                )

    # Compute auxiliary tables `outside` based on the newly-computed
    # values in `insideSkip`
    for orig_species in species_lca.tree.traverse("preorder"):
        for kind in SyntenyAssignment:
            if orig_species.is_leaf():
                continue

            spe_left, spe_right = orig_species.children

            tables.outside[orig_object][spe_left][kind].update(
                *tables.outside[orig_object][orig_species][kind],
                *tables.insideSkip[orig_object][spe_right][kind],
            )

            tables.outside[orig_object][spe_right][kind].update(
                *tables.outside[orig_object][orig_species][kind],
                *tables.insideSkip[orig_object][spe_left][kind],
            )


def _compute_table(
    srec_input: SuperReconciliationInput,
    minimal_sets: Dict[TreeNode, UnorderedSynteny],
    allowed_species: Callable[[Tree, TreeNode], Iterable[TreeNode]],
    retention_policy: RetentionPolicy,
) -> AssignmentTable:
    """
    Compute an assignment table that can be used to construct minimum-cost
    unordered super-reconciliations. Each entry `table[object][species][kind]`
    of the resulting table contains the possible assignments of `object`’s
    children that lead to minimum-cost unordered super-reconciliations of the
    subtrees rooted at `object` such that the root is mapped to `species` and
    to the minimal synteny or to a larger synteny depending on the value of
    `kind`.

    :param srec_input: objects of the unordered super-reconciliation
    :param minimal_sets: set of genes in the minimal synteny of each object
    :param allowed_species: callable that gives the set of allowed species for
        each object
    :param retention_policy: whether to keep any minimal assignment in each
        entry or all minimal assignments
    :returns: computed assignment table
    """
    tables = _initialize_tables(srec_input, retention_policy)

    for orig_object in tqdm(
        srec_input.object_tree.traverse("postorder"),
        desc="Table entries",
        total=sum(1 for _ in srec_input.object_tree.traverse()),
        ascii=True,
        leave=False,
    ):
        if orig_object.is_leaf():
            continue

        _compute_tables_entry(
            srec_input.species_lca,
            orig_object,
            minimal_sets,
            tables,
            srec_input.costs,
        )

    return tables.all


def _decode_table(
    root_object: TreeNode,
    root_species: TreeNode,
    root_kind: SyntenyAssignment,
    ancestor_synteny: Optional[UnorderedSynteny],
    srec_input: SuperReconciliationInput,
    gain_sets: Dict[TreeNode, UnorderedSynteny],
    minimal_sets: Dict[TreeNode, UnorderedSynteny],
    table: AssignmentTable,
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
        :func:`_compute_table`
    :returns: yields minimum-cost unordered super-reconciliations
    """
    if root_kind == SyntenyAssignment.MINIMAL:
        ancestor_synteny = minimal_sets[root_object]
        root_synteny = sort_synteny(minimal_sets[root_object])
    else:
        ancestor_synteny |= gain_sets[root_object]
        root_synteny = sort_synteny(ancestor_synteny)

    if (
        root_object.is_leaf()
        and not table[root_object][root_species][root_kind].is_infinite()
    ):
        # pylint has trouble seeing the attributes from the descendant dataclass
        yield SuperReconciliationOutput( # pylint: disable=unexpected-keyword-arg
            input=srec_input,
            object_species={root_object: root_species},
            syntenies={root_object: root_synteny},
            ordered=False,
        )
        return

    for info in table[root_object][root_species][root_kind].infos():
        left_object, right_object = root_object.children
        mappings = product(
            _decode_table(
                left_object,
                info.left.species,
                info.left.synteny,
                ancestor_synteny,
                srec_input,
                gain_sets,
                minimal_sets,
                table,
            ),
            _decode_table(
                right_object,
                info.right.species,
                info.right.synteny,
                ancestor_synteny,
                srec_input,
                gain_sets,
                minimal_sets,
                table,
            ),
        )

        for map_left, map_right in mappings:
            # pylint has trouble seeing the attributes from the descendant dataclass
            yield SuperReconciliationOutput( # pylint: disable=unexpected-keyword-arg
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


def _superrec(
    srec_input: SuperReconciliationInput,
    policy: RetentionPolicy,
    allowed_species: Callable[[TreeNode], Iterable[TreeNode]],
) -> Set[SuperReconciliationOutput]:
    results: Entry[int, SuperReconciliationOutput] = Entry(
        MergePolicy.MIN, policy
    )

    for srec_input_bin in tqdm(
        srec_input.binarize(),
        desc="Multifurcation resolutions",
        total=sum(1 for _ in srec_input.binarize()),
        ascii=True,
    ):
        srec_input_bin.label_internal()
        gain_sets = _compute_gain_sets(srec_input_bin)
        minimal_sets = _compute_minimal_sets(srec_input_bin, gain_sets)
        synteny_tree = srec_input_bin.object_tree

        table = _compute_table(
            srec_input_bin,
            minimal_sets,
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
                    _decode_table(
                        synteny_tree,
                        root_species,
                        SyntenyAssignment.MINIMAL,
                        minimal_sets[synteny_tree],
                        srec_input_bin,
                        gain_sets,
                        minimal_sets,
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
    return _superrec(
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
    return _superrec(
        srec_input,
        policy,
        allowed_species=lambda species, _: species.traverse("postorder"),
    )
