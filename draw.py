#!/usr/bin/env python3
"""
Draw a representation of the (super-)reconciliation of two trees.

Input arguments OBJECT_TREE, SPECIES_TREE, RECONCILIATION and SYNTENIES can
be passed either as plain strings or as a path to a file.
"""
import argparse
import textwrap
import sys
from ete3 import Tree
from superrec2.draw import compute_layout, render_to_tikz
from superrec2.model.synteny import parse_synteny_mapping
from superrec2.model.tree_mapping import (
    get_species_mapping, parse_tree_mapping
)
from superrec2.model.reconciliation import (
    SuperReconciliationInput, SuperReconciliationOutput
)
from superrec2.utils.trees import LowestCommonAncestor
from superrec2.utils.tex import xelatex_compile, TeXError


def parse_arguments():
    """Retrieve arguments or show help."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "object_tree",
        metavar="OBJECT_TREE",
        help="embedded object tree (Newick format)",
    )
    parser.add_argument(
        "species_tree",
        metavar="SPECIES_TREE",
        help="host species tree (Newick format)",
    )
    parser.add_argument(
        "reconciliation",
        metavar="RECONCILIATION",
        help="mapping of internal genes onto species",
    )
    parser.add_argument(
        "syntenies",
        metavar="SYNTENIES",
        nargs="?",
        help="synteny labeling of gene nodes",
    )
    parser.add_argument(
        "--output",
        metavar="PATH",
        default="-",
        help="where to output the result (default: output to stdout)",
    )
    return parser.parse_args()


def generate_tikz(args):
    """Generate TikZ code corresponding to the given reconciliation."""
    object_tree = Tree(args.object_tree, format=1)
    species_tree = Tree(args.species_tree, format=1)

    mapping = {
        **get_species_mapping(object_tree, species_tree),
        **parse_tree_mapping(object_tree, species_tree, args.reconciliation),
    }
    syntenies = (
        parse_synteny_mapping(object_tree, args.syntenies)
        if args.syntenies is not None
        else {}
    )
    srec_input = SuperReconciliationInput(
        object_tree=object_tree,
        species_lca=LowestCommonAncestor(species_tree),
        leaf_object_species={},
        costs={},
        leaf_syntenies={},
    )
    srec_output = SuperReconciliationOutput(
        input=srec_input,
        object_species=mapping,
        syntenies=syntenies,
    )
    layout_info = compute_layout(srec_output)
    return render_to_tikz(srec_output, layout_info)


def output(args, tikz):
    """Generate output."""
    if args.output == "-":
        print(tikz, end="")
    elif args.output.endswith(".tex"):
        with open(args.output, "w") as tex_source:
            tex_source.write(tikz)
    elif args.output.endswith(".pdf"):
        try:
            xelatex_compile(
                source=textwrap.dedent(
                    r"""
                    \documentclass[crop, tikz, border=20pt]{standalone}
                    \usepackage{tikz}
                    \usetikzlibrary{arrows.meta}
                    \usetikzlibrary{shapes}
                    \begin{document}
                    \scrollmode
                    """
                ).lstrip()
                + tikz
                + textwrap.dedent(
                    r"""
                    \batchmode
                    \end{document}
                    """
                ).lstrip(),
                dest=args.output,
            )
        except TeXError as err:
            print(f"XeLaTeX returned an error (code: {err.code})")
            print("Output from the compiler:")

            for line in err.message.splitlines():
                print(f"> {line}")

            return 1
    else:
        raise RuntimeError(f"Unrecognized file extension: {args.output}")

    return 0


def main():  # pylint:disable=missing-function-docstring
    args = parse_arguments()
    tikz = generate_tikz(args)
    return output(args, tikz)


if __name__ == "__main__":
    sys.exit(main())
