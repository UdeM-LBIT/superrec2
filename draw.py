#!/usr/bin/env python3
"""Draw a representation of a reconciliation."""
import argparse
import textwrap
import sys
from ete3 import PhyloTree
from superrec2.reconciliation.draw import layout, render_to_tikz
from superrec2.reconciliation.tools import (
    get_species_name,
    reconcile_leaves,
    parse_reconciliation,
    parse_labeling,
)
from superrec2.utils.tex import xelatex_compile, TeXError

# Retrieve arguments or show help
parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument(
    "gene_tree",
    metavar="GENE_TREE",
    help="embedded gene tree in Newick format",
)

parser.add_argument(
    "species_tree",
    metavar="SPECIES_TREE",
    help="host species tree in Newick format",
)

parser.add_argument(
    "reconciliation",
    metavar="RECONCILIATION",
    help="mapping of internal genes onto species",
)

parser.add_argument(
    "labeling",
    metavar="LABELING",
    nargs="?",
    help="synteny labeling of gene nodes",
)

parser.add_argument(
    "--output",
    metavar="PATH",
    default="-",
    help="where to output the result (default: output to stdout)",
)

args = parser.parse_args()

# Parse trees, reconciliation, and labeling
gene_tree = PhyloTree(
    args.gene_tree,
    sp_naming_function=get_species_name, format=1
)
species_tree = PhyloTree(args.species_tree, format=1)

rec = {
    **reconcile_leaves(gene_tree, species_tree),
    **parse_reconciliation(gene_tree, species_tree, args.reconciliation),
}

labeling = (
    parse_labeling(gene_tree, args.labeling)
    if args.labeling is not None
    else {}
)

# Generate TikZ code
layout_info = layout(gene_tree, species_tree, rec)
tikz = render_to_tikz(species_tree, rec, layout_info, labeling)

# Generate output
if args.output == "-":
    print(tikz, end="")
elif args.output.endswith(".tex"):
    with open(args.output, "w") as tex_source:
        tex_source.write(tikz)
elif args.output.endswith(".pdf"):
    try:
        xelatex_compile(
            source=textwrap.dedent(r"""
                \documentclass[crop, tikz, border=20pt]{standalone}
                \usepackage{tikz}
                \usetikzlibrary{arrows.meta}
                \usetikzlibrary{shapes}
                \begin{document}
                \scrollmode
            """).lstrip()
            + tikz
            + textwrap.dedent(r"""
                \batchmode
                \end{document}
            """).lstrip(),
            dest=args.output,
        )
    except TeXError as err:
        print(f"XeLaTeX returned an error (code: {err.code})")
        print("Output from the compiler:")

        for line in err.message.splitlines():
            print(f"> {line}")

        sys.exit(1)
else:
    raise RuntimeError(f"Unrecognized file extension: {args.output}")
