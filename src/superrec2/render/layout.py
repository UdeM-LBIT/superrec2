"""Compute layouts for reconciliations."""
from typing import Any, Dict
from ete3 import Tree, TreeNode
from .tikz import measure_nodes
from .model import (
    Branch,
    GeneAnchor,
    DrawParams,
    Orientation,
    Layout,
    SubtreeLayout,
    PseudoGene,
)
from ..model.reconciliation import (
    NodeEvent,
    EdgeEvent,
    ReconciliationOutput,
    SuperReconciliationOutput,
)
from ..model.synteny import format_synteny
from ..utils import tex
from ..utils.geometry import Position, Rect, Size


def _add_losses(
    layout_state: Dict[TreeNode, Dict],
    gene: TreeNode,
    start_species: TreeNode,
    end_species: TreeNode,
) -> GeneAnchor:
    """
    Insert virtual gene loss nodes between
    a parent species and a child species.

    :param gene: parent gene that is lost
    :param start_species: lower species in which the gene is conserved
    :param end_species: parent of the species from which the
        gene originated
    :returns: first virtual child created in the process
    """
    prev_gene = gene
    color = getattr(gene, "color", None)
    prev_species = start_species
    start_species = start_species.up

    while start_species != end_species:
        is_left = prev_species == start_species.children[0]
        is_right = prev_species == start_species.children[1]
        state = layout_state[start_species]
        cur_gene = PseudoGene()

        state["anchor_nodes"].add(cur_gene)
        state["branches"][cur_gene] = {
            "kind": EdgeEvent.FULL_LOSS,
            "name": "",
            "left": prev_gene if is_left else None,
            "right": prev_gene if is_right else None,
        }

        if color is not None:
            state["branches"][cur_gene]["color"] = color

        prev_gene = cur_gene
        prev_species = start_species
        start_species = start_species.up

    return prev_gene


def _compute_branches(  # pylint:disable=too-many-locals
    layout_state: Dict[TreeNode, Dict],
    rec: ReconciliationOutput,
    params: DrawParams,
) -> None:
    """Create the branching nodes for each species."""
    gene_tree = rec.input.object_tree
    species_lca = rec.input.species_lca
    species_tree = species_lca.tree
    mapping = rec.object_species
    syntenies = rec.syntenies if isinstance(rec, SuperReconciliationOutput) else {}

    # Propagate color feature downwards in the tree
    last_color = None
    last_color_node = None

    for root_gene in gene_tree.traverse("preorder"):
        if hasattr(root_gene, "color"):
            last_color = root_gene.color
            last_color_node = root_gene
        elif last_color_node is not None:
            if last_color_node in root_gene.iter_ancestors():
                root_gene.add_feature("color", last_color)
            else:
                last_color = None
                last_color_node = None

    # Find gene tree nodes associated to each species and create branches
    for root_species in species_tree.traverse("postorder"):
        state: Dict[str, Any] = {
            "branches": {},
            "anchor_nodes": set(),
        }
        layout_state[root_species] = state

        for root_gene in gene_tree.traverse("postorder"):
            if mapping[root_gene] != root_species:
                continue

            synteny = (
                format_synteny(
                    map(tex.escape, syntenies[root_gene]),
                    params.event_label_width,
                ).replace("\n", "\\\\")
                if root_gene in syntenies
                else ""
            )
            equal_to_parent = syntenies.get(root_gene) == syntenies.get(root_gene.up)

            if root_gene.is_leaf():
                # Create branches even for leaf genes
                if synteny:
                    name = synteny
                elif root_gene.name:
                    species_name, gene_name = root_gene.name.rsplit("_", 1)
                    name = (
                        rf"{tex.escape(species_name)}"
                        rf"\textsubscript{{{tex.escape(gene_name)}}}"
                    )
                else:
                    name = ""

                state["anchor_nodes"].add(root_gene)
                state["branches"][root_gene] = {
                    "kind": NodeEvent.LEAF,
                    "name": name,
                }
            else:
                # Create branches for actual internal nodes
                left_gene, right_gene = root_gene.children
                event = rec.node_event(root_gene)
                name = synteny if not equal_to_parent else ""

                if event == NodeEvent.SPECIATION:
                    # Speciation nodes are located below the trunk
                    # and linked to child species’s gene anchors
                    left_species = root_species.children[0]

                    if species_lca.is_ancestor_of(left_species, mapping[right_gene]):
                        # Left gene and right gene are swapped relative
                        # to the left and right species
                        left_gene, right_gene = right_gene, left_gene

                    left_gene = _add_losses(
                        layout_state,
                        left_gene,
                        mapping[left_gene],
                        root_species,
                    )
                    right_gene = _add_losses(
                        layout_state,
                        right_gene,
                        mapping[right_gene],
                        root_species,
                    )

                    state["anchor_nodes"].add(root_gene)
                    state["branches"][root_gene] = {
                        "kind": NodeEvent.SPECIATION,
                        "name": name,
                        "left": left_gene,
                        "right": right_gene,
                    }
                elif event == NodeEvent.DUPLICATION:
                    # Duplications are located in the trunk and linked
                    # to other nodes in the same species
                    left_gene = _add_losses(
                        layout_state,
                        left_gene,
                        mapping[left_gene],
                        root_species.up,
                    )
                    right_gene = _add_losses(
                        layout_state,
                        right_gene,
                        mapping[right_gene],
                        root_species.up,
                    )

                    state["anchor_nodes"].add(root_gene)
                    state["anchor_nodes"].remove(left_gene)
                    state["anchor_nodes"].remove(right_gene)
                    state["branches"][root_gene] = {
                        "kind": NodeEvent.DUPLICATION,
                        "name": name,
                        "left": left_gene,
                        "right": right_gene,
                    }
                elif event == NodeEvent.HORIZONTAL_TRANSFER:
                    # Transfers are located in the trunk, like duplications,
                    # but are linked to a node outside the current subtree
                    conserv_gene, foreign_gene = (
                        (left_gene, right_gene)
                        if species_lca.is_ancestor_of(root_species, mapping[left_gene])
                        else (right_gene, left_gene)
                    )
                    conserv_gene = _add_losses(
                        layout_state,
                        conserv_gene,
                        mapping[conserv_gene],
                        root_species.up,
                    )

                    state["anchor_nodes"].add(root_gene)
                    state["anchor_nodes"].remove(conserv_gene)
                    state["branches"][root_gene] = {
                        "kind": NodeEvent.HORIZONTAL_TRANSFER,
                        "name": name,
                        "left": conserv_gene,
                        "right": foreign_gene,
                    }
                else:
                    raise ValueError("Invalid event")

            if hasattr(root_gene, "color"):
                state["branches"][root_gene]["color"] = root_gene.color


def _layout_branches(  # pylint:disable=too-many-locals
    layout_state: Dict[TreeNode, Dict],
    species_tree: Tree,
    params: DrawParams,
):
    """Compute the size and relative position of each branch."""
    branch_nodes = []
    branches = []

    for layout in layout_state.values():
        for node, branch in layout["branches"].items():
            branch_nodes.append(node)
            branches.append((branch["kind"], branch["name"]))

    measures = dict(zip(branch_nodes, measure_nodes(branches, params)))

    for root_species in species_tree.traverse():
        next_pos_across: float = 0
        next_pos_sequence = params.species_branch_padding
        layout = layout_state[root_species]
        layout["anchors"] = {}

        for root_gene, branch in layout["branches"].items():
            size = (
                measures[root_gene].overall_size()
                if root_gene in measures
                else Size(0, 0)
            )

            if branch["kind"] == NodeEvent.LEAF:
                if params.orientation == Orientation.VERTICAL:
                    next_pos_across -= size.w
                    pos = Position(next_pos_across, -size.h)
                else:
                    next_pos_across -= size.h
                    pos = Position(-size.w, next_pos_across)

                next_pos_across -= params.gene_branch_spacing
            elif branch["kind"] in (
                NodeEvent.SPECIATION,
                EdgeEvent.FULL_LOSS,
            ):
                if params.orientation == Orientation.VERTICAL:
                    next_pos_across -= size.w
                    pos = Position(next_pos_across, next_pos_sequence)
                    next_pos_sequence += size.h
                else:
                    next_pos_across -= size.h
                    pos = Position(next_pos_sequence, next_pos_across)
                    next_pos_sequence += size.w

                next_pos_across -= params.gene_branch_spacing
                next_pos_sequence += params.gene_branch_spacing
            elif branch["kind"] == NodeEvent.DUPLICATION:
                left_rect = layout["branches"][branch["left"]]["rect"]
                right_rect = layout["branches"][branch["right"]]["rect"]

                if params.orientation == Orientation.VERTICAL:
                    across = ((left_rect.center() + right_rect.center()).x - size.w) / 2
                    sequence = (
                        min(
                            params.species_branch_padding,
                            left_rect.y,
                            right_rect.y,
                        )
                        - params.species_branch_padding
                        - size.h
                    )
                    pos = Position(across, sequence)
                else:
                    across = ((left_rect.center() + right_rect.center()).y - size.h) / 2
                    sequence = (
                        min(
                            params.species_branch_padding,
                            left_rect.x,
                            right_rect.x,
                        )
                        - params.species_branch_padding
                        - size.w
                    )
                    pos = Position(sequence, across)
            elif branch["kind"] == NodeEvent.HORIZONTAL_TRANSFER:
                cons_rect = layout["branches"][branch["left"]]["rect"]

                if params.orientation == Orientation.VERTICAL:
                    across = cons_rect.center().x - size.w / 2
                    sequence = (
                        min(params.species_branch_padding, cons_rect.y)
                        - params.species_branch_padding
                        - size.h
                    )
                    pos = Position(across, sequence)
                else:
                    across = cons_rect.center().y - size.h / 2
                    sequence = (
                        min(params.species_branch_padding, cons_rect.x)
                        - params.species_branch_padding
                        - size.w
                    )
                    pos = Position(sequence, across)
            else:
                raise ValueError("Invalid node type")

            rect = Rect.make_from(pos, size)
            branch["rect"] = rect

            if root_gene in layout["anchor_nodes"]:
                if params.orientation == Orientation.VERTICAL:
                    layout["anchors"][root_gene] = Position(rect.center().x, 0)
                else:
                    layout["anchors"][root_gene] = Position(0, rect.center().y)

        del layout["anchor_nodes"]

        # Shift all nodes to the left or up to make room
        # for the initial padding
        if layout["branches"]:
            if params.orientation == Orientation.VERTICAL:
                padding_shift = Position(
                    x=(
                        min(
                            -branch["rect"].right().x
                            for branch in layout["branches"].values()
                        )
                        - params.species_branch_padding
                    ),
                    y=0,
                )
            else:
                padding_shift = Position(
                    x=0,
                    y=(
                        min(
                            -branch["rect"].bottom().y
                            for branch in layout["branches"].values()
                        )
                        - params.species_branch_padding
                    ),
                )

            for root_gene in layout["branches"]:
                branch = layout["branches"][root_gene]
                branch["rect"] += padding_shift

            for root_gene in layout["anchors"]:
                layout["anchors"][root_gene] += padding_shift


def _layout_subtrees(
    layout_state: Dict[TreeNode, Dict],
    species_tree: Tree,
    params: DrawParams,
):
    """Compute the size and absolute position of each subtree."""
    # Compute the size of each subtree
    for root_species in species_tree.traverse("postorder"):
        state = layout_state[root_species]

        if state["branches"]:
            if params.orientation == Orientation.VERTICAL:
                trunk_width = (
                    max(
                        -branch["rect"].top_left().x
                        for branch in state["branches"].values()
                    )
                    + params.species_branch_padding
                )
                trunk_height = (
                    max(
                        0,
                        max(
                            -branch["rect"].top_left().y
                            for branch in state["branches"].values()
                        ),
                    )
                    + params.trunk_overhead
                )
                fork_thickness = (
                    max(
                        0,
                        max(
                            branch["rect"].bottom_right().y
                            for branch in state["branches"].values()
                        ),
                    )
                    + params.species_branch_padding
                )
            else:
                trunk_width = (
                    max(
                        0,
                        max(
                            -branch["rect"].top_left().x
                            for branch in state["branches"].values()
                        ),
                    )
                    + params.trunk_overhead
                )
                trunk_height = (
                    max(
                        -branch["rect"].top_left().y
                        for branch in state["branches"].values()
                    )
                    + params.species_branch_padding
                )
                fork_thickness = (
                    max(
                        0,
                        max(
                            branch["rect"].bottom_right().x
                            for branch in state["branches"].values()
                        ),
                    )
                    + params.species_branch_padding
                )
        else:
            # Empty subtree
            fork_thickness = 0

            if params.orientation == Orientation.VERTICAL:
                trunk_width = 0
                trunk_height = params.trunk_overhead
            else:
                trunk_width = params.trunk_overhead
                trunk_height = 0

        trunk_size = Size(trunk_width, trunk_height)

        if root_species.is_leaf():
            # Extant species
            state["size"] = trunk_size
            state["trunk"] = Rect.make_from(Position(0, 0), trunk_size)
            state["fork_thickness"] = 0
        else:
            # Ancestral species
            left_species, right_species = root_species.children
            left_info = layout_state[left_species]
            right_info = layout_state[right_species]

            if params.orientation == Orientation.VERTICAL:
                subtree_span = (
                    max(left_info["size"].h, right_info["size"].h) + trunk_height
                )
            else:
                subtree_span = (
                    max(left_info["size"].w, right_info["size"].w) + trunk_width
                )

            subtree_span += params.level_spacing + fork_thickness

            if params.orientation == Orientation.VERTICAL:
                left_trunk_dist = left_info["size"].w - left_info["trunk"].right().x
                right_trunk_dist = right_info["trunk"].left().x
                subtree_spacing = max(
                    trunk_width - (left_trunk_dist + right_trunk_dist),
                    params.min_subtree_spacing,
                )

                state["size"] = Size(
                    left_info["size"].w + subtree_spacing + right_info["size"].w,
                    subtree_span,
                )
                state["left_pos"] = Position(
                    0,
                    subtree_span - left_info["size"].h,
                )
                state["right_pos"] = Position(
                    left_info["size"].w + subtree_spacing,
                    subtree_span - right_info["size"].h,
                )
                trunk_pos = Position(
                    left_info["size"].w + (subtree_spacing - trunk_width) / 2,
                    0,
                )
            else:
                left_trunk_dist = left_info["size"].h - left_info["trunk"].bottom().y
                right_trunk_dist = right_info["trunk"].top().y
                subtree_spacing = max(
                    trunk_height - (left_trunk_dist + right_trunk_dist),
                    params.min_subtree_spacing,
                )

                state["size"] = Size(
                    subtree_span,
                    left_info["size"].h + subtree_spacing + right_info["size"].h,
                )
                state["left_pos"] = Position(
                    subtree_span - left_info["size"].w,
                    0,
                )
                state["right_pos"] = Position(
                    subtree_span - right_info["size"].w,
                    left_info["size"].h + subtree_spacing,
                )
                trunk_pos = Position(
                    0,
                    left_info["size"].h + (subtree_spacing - trunk_height) / 2,
                )

            state["trunk"] = Rect.make_from(trunk_pos, trunk_size)
            state["fork_thickness"] = fork_thickness

    # Compute the absolute position of each subtree
    layout_state[species_tree]["rect"] = Rect.make_from(
        position=Position(0, 0),
        size=layout_state[species_tree]["size"],
    )
    del layout_state[species_tree]["size"]

    for root_species in species_tree.traverse("preorder"):
        this_layout = layout_state[root_species]
        this_rect = this_layout["rect"]

        # Position child subtrees
        if not root_species.is_leaf():
            left_species, right_species = root_species.children

            layout_state[left_species]["rect"] = Rect.make_from(
                position=this_rect.top_left() + this_layout["left_pos"],
                size=layout_state[left_species]["size"],
            )
            del this_layout["left_pos"]
            del layout_state[left_species]["size"]

            layout_state[right_species]["rect"] = Rect.make_from(
                position=this_rect.top_left() + this_layout["right_pos"],
                size=layout_state[right_species]["size"],
            )
            del this_layout["right_pos"]
            del layout_state[right_species]["size"]

        # Make trunk, anchor, and branch nodes positions absolute
        # and compute branch anchors
        this_layout["trunk"] += this_rect.top_left()
        trunk_rect = this_layout["trunk"]

        for anchor in this_layout["anchors"]:
            if params.orientation == Orientation.VERTICAL:
                this_layout["anchors"][anchor] += trunk_rect.top_right()
            else:
                this_layout["anchors"][anchor] += trunk_rect.bottom_left()

        for branch in this_layout["branches"].values():
            branch["rect"] += trunk_rect.bottom_right()
            branch_rect = branch["rect"]

            if branch["kind"] == EdgeEvent.FULL_LOSS:
                branch["anchor_parent"] = branch_rect.center()
                branch["anchor_left"] = branch_rect.center()
                branch["anchor_right"] = branch_rect.center()
                branch["anchor_child"] = branch_rect.center()
            else:
                if params.orientation == Orientation.VERTICAL:
                    branch["anchor_parent"] = branch_rect.top()
                    branch["anchor_left"] = branch_rect.left()
                    branch["anchor_right"] = branch_rect.right()
                    branch["anchor_child"] = branch_rect.bottom()
                else:
                    branch["anchor_parent"] = branch_rect.left()
                    branch["anchor_left"] = branch_rect.top()
                    branch["anchor_right"] = branch_rect.bottom()
                    branch["anchor_child"] = branch_rect.right()


def _finalize_layout(
    layout_state: Dict[TreeNode, Dict],
    species_tree: Tree,
) -> Layout:
    """Turn a computed layout into final immutable structures."""
    result: Dict[TreeNode, SubtreeLayout] = {}

    for root_species in species_tree.traverse("preorder"):
        this_layout = layout_state[root_species]

        for anchor, branch in this_layout["branches"].items():
            this_layout["branches"][anchor] = Branch(**branch)

        result[root_species] = SubtreeLayout(**layout_state[root_species])

    return result


def compute(
    rec: ReconciliationOutput,
    params: DrawParams = DrawParams(),
) -> Layout:
    """
    Compute the layout of a gene tree embedded in a species tree.

    :param rec: reconciliation object defining the gene and species trees
        their embedding, and an optional synteny labelling
    :param params: layout parameters
    :returns: layout information for each species node
    """
    layout_state: Dict[TreeNode, Dict] = {}
    species_tree = rec.input.species_lca.tree
    _compute_branches(layout_state, rec, params)
    _layout_branches(layout_state, species_tree, params)
    _layout_subtrees(layout_state, species_tree, params)
    return _finalize_layout(layout_state, species_tree)
