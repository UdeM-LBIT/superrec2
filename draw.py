#!/usr/bin/env python3
"""Draw a representation of the (super-)reconciliation of two trees."""
import argparse
import json
import textwrap
import sys
from ete3 import Tree
from superrec2.draw import compute_layout, render_to_tikz, DrawParams, Orientation
from superrec2.model.synteny import parse_synteny_mapping
from superrec2.model.tree_mapping import get_species_mapping, parse_tree_mapping
from superrec2.model.reconciliation import (
    SuperReconciliationInput,
    SuperReconciliationOutput,
)
from superrec2.utils.trees import LowestCommonAncestor
from superrec2.utils.tex import xelatex_compile, TeXError


def parse_arguments():
    """Retrieve arguments or show help."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output_type",
        metavar="TYPE",
        nargs="?",
        choices=("tikz", "pdf"),
        help="kind of output to generate (default: guess based on output file \
extension, or 'tikz' for output to stdout)",
    )
    parser.add_argument(
        "--input",
        metavar="PATH",
        default="-",
        help="path to a file defining the reconciliation input \
(default: read from stdin)",
    )
    parser.add_argument(
        "--output",
        metavar="PATH",
        default="-",
        help="path where the result will be stored (default: output to stdout)",
    )
    parser.add_argument(
        "--orientation",
        choices=("vertical", "horizontal"),
        default="vertical",
        help="growing direction of the tree (default: %(default)s)",
    )
    return parser.parse_args()


def generate_tikz(args):
    """Generate TikZ code corresponding to the given reconciliation."""
    infile = open(args.input, "r") if args.input != "-" else sys.stdin
    data = json.load(infile)
    params = DrawParams(
        orientation=Orientation[args.orientation.upper()],
    )
    rec_output = SuperReconciliationOutput.from_dict(data)
    layout_info = compute_layout(rec_output, params)
    return render_to_tikz(rec_output, layout_info, params)


def output(args, tikz):
    """Generate output."""
    outfile = (
        open(args.output, "wb") if args.output != "-" else sys.stdout.buffer
    )
    output_type = args.output_type

    if output_type is None:
        if args.output == "-" or args.output.endswith(".tex"):
            output_type = "tikz"
        elif args.output.endswith(".pdf"):
            output_type = "pdf"
        else:
            print(
                "Error: Unknown file extension, please specify output \
type explicitly",
                file=sys.stderr,
            )
            return 1

    if output_type == "tikz":
        outfile.write(tikz.encode())
    elif output_type == "pdf":
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
                dest=outfile,
            )
        except TeXError as err:
            print(f"XeLaTeX returned an error (code: {err.code})")
            print("Output from the compiler:")

            for line in err.message.splitlines():
                print(f"> {line}")

            return 1

    return 0


def main():  # pylint:disable=missing-function-docstring
    args = parse_arguments()
    tikz = generate_tikz(args)
    return output(args, tikz)


if __name__ == "__main__":
    sys.exit(main())
