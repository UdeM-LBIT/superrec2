#!/usr/bin/env python3
"""Draw a representation of a reconciliation."""
import argparse
import os
import subprocess
import tempfile
import textwrap
import shutil
from ete3 import PhyloTree
from superrec2.reconciliation_draw import layout, render_to_tikz
from superrec2.reconciliation import (
    get_species_name,
    reconcile_leaves,
    parse_reconciliation,
)

# Retrieve arguments or show help
parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument(
    "gene_tree",
    metavar="GENE_TREE",
    help="embedded gene tree in Newick format"
)

parser.add_argument(
    "species_tree",
    metavar="SPECIES_TREE",
    help="host species tree in Newick format"
)

parser.add_argument(
    "reconciliation",
    metavar="REC",
    help="mapping of internal genes onto species"
)

parser.add_argument(
    "--output",
    metavar="PATH",
    default="-",
    help="where to output the result (default: output to stdout)",
)

args = parser.parse_args()

# Parse trees and reconciliation
gene_tree = PhyloTree(
    args.gene_tree,
    sp_naming_function=get_species_name, format=1
)
species_tree = PhyloTree(args.species_tree, format=1)
rec = {
    **reconcile_leaves(gene_tree, species_tree),
    **parse_reconciliation(gene_tree, species_tree, args.reconciliation),
}

# Generate TikZ code
layout_info = layout(gene_tree, species_tree, rec)
tikz = render_to_tikz(species_tree, rec, layout_info)

# Generate output
if args.output == "-":
    print(tikz, end="")
elif args.output.endswith(".tex"):
    with open(args.output, "w") as tex_source:
        tex_source.write(tikz)
elif args.output.endswith(".pdf"):
    with tempfile.TemporaryDirectory() as tmpdir:
        tex_source_path = os.path.join(tmpdir, "rec.tex")
        pdf_gen_path = os.path.join(tmpdir, "rec.pdf")

        with open(tex_source_path, "w") as tex_source:
            tex_source.write(textwrap.dedent(r"""
                \documentclass[crop, tikz, border=20pt]{standalone}

                \usepackage{tikz}
                \usetikzlibrary{arrows.meta}
                \usetikzlibrary{shapes}
                \usetikzlibrary{decorations.pathreplacing}

                \begin{document}
            """).lstrip())
            tex_source.write(tikz)
            tex_source.write("\\end{document}\n")

        result = subprocess.run(
            ["xelatex", tex_source_path],
            cwd=tmpdir,
        )

        if result.returncode == 0:
            shutil.move(pdf_gen_path, args.output)
else:
    raise RuntimeError(f"Unrecognized file extension: {args.output}")
