"""Draw a representation of the (super-)reconciliation of two trees."""
import json
import textwrap
import sys
from .util import add_arg_input, add_arg_output
from superrec2.render import layout, tikz
from superrec2.render.layout import (
    DrawParams,
    Orientation,
)
from superrec2.model.reconciliation import (
    ReconciliationOutput,
    SuperReconciliationOutput,
)
from superrec2.utils.tex import tex_compile, TeXError


def generate_tikz(args):
    """Generate TikZ code corresponding to the given reconciliation."""
    data = json.load(args.input)
    params = DrawParams(orientation=Orientation[args.orientation.upper()])

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
        if args.output.name == "-" or args.output.name.endswith(".tex"):
            output_type = "tikz"
        elif args.output.name.endswith(".pdf"):
            output_type = "pdf"
        else:
            print(
                "Error: Unknown file extension, please specify output \
type explicitly",
                file=sys.stderr,
            )
            return 1

    if output_type == "tikz":
        args.output.write(tikz_code.encode())
    elif output_type == "pdf":
        try:
            tex_compile(
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
                dest=args.output,
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
    add_arg_input(subparser, "file defining the history to draw")
    add_arg_output(subparser, "file where the resulting drawing will be stored", "wb")
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
        default="horizontal",
        help="growing direction of the tree (default: %(default)s)",
    )
    subparser.set_defaults(func=draw)
