"""Draw a representation of the (super-)reconciliation of two trees."""
import json
import textwrap
import sys
from .util import add_arg_input, add_arg_output, open_std
from superrec2.render import layout, tikz
from superrec2.render.layout import (
    DrawParams,
    Orientation,
)
from superrec2.model.reconciliation import (
    ReconciliationOutput,
    SuperReconciliationOutput,
)
from superrec2.utils.tex import xelatex_compile, TeXError


def generate_tikz(args):
    """Generate TikZ code corresponding to the given reconciliation."""
    with open_std(args.input, "r", encoding="utf8") as infile_handle:
        data = json.load(infile_handle)

    params = DrawParams(
        orientation=Orientation[args.orientation.upper()],
    )

    if "syntenies" in data:
        rec_output = SuperReconciliationOutput.from_dict(data)
    else:
        rec_output = ReconciliationOutput.from_dict(data)

    layout_info = layout.compute(rec_output, params)
    return tikz.render(rec_output, layout_info, params)


def output(args, tikz_code):
    """Generate output."""
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
        with open_std(args.output, "w", encoding="utf8") as outfile:
            outfile.write(tikz_code)
    elif output_type == "pdf":
        try:
            with open_std(args.output, "wb") as outfile:
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
                    + tikz_code
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


def draw(args):
    """Run the drawing subcommand with the given arguments."""
    tikz_code = generate_tikz(args)
    return output(args, tikz_code)


def add_args(parser):
    """Add the drawing subcommand to a command-line argument parser."""
    subparser = parser.add_parser("draw", description=__doc__)
    add_arg_input(subparser, "a file defining the reconciliation result to draw")
    add_arg_output(subparser, "where the resulting drawing will be stored")
    subparser.add_argument(
        "output_type",
        metavar="TYPE",
        nargs="?",
        choices=("tikz", "pdf"),
        help="kind of output to generate (default: guess based on output file \
extension, or 'tikz' for output to stdout)",
    )
    subparser.add_argument(
        "--orientation",
        choices=("vertical", "horizontal"),
        default="vertical",
        help="growing direction of the tree (default: %(default)s)",
    )
    subparser.set_defaults(func=draw)
