"""Draw a representation of the (super-)reconciliation of two trees."""

import json
import textwrap
import sys
from sowing import traversal
from .util import add_arg_input, add_arg_output
from ..model.history import History
from ..render import layout, tikz
from ..render.model import DrawParams, Orientation
from ..utils.tex import tex_compile, TeXError


def generate_tex(args) -> str:
    """Generate TeX code corresponding to the given reconciliation."""
    for line in args.input:
        try:
            history = History.from_mapping(json.loads(line))
            break
        except json.JSONDecodeError:
            pass
    else:
        raise RuntimeError("No valid history found in input")

    history.validate()

    params = DrawParams(
        debug=args.debug,
        orientation=Orientation[args.orientation.title()],
    )

    rects = tikz.measure_events(
        events=(
            (cursor.node.data, history.host_index[cursor.node.data.host].node.data)
            for cursor in traversal.depth(history.event_tree)
        ),
        params=params,
    )

    result = layout.compute(history, rects, params)
    return tikz.render(result, params)


def output(args, tex_code) -> int:
    """Generate output."""
    output_type = args.output_type

    if output_type is None:
        if args.output.name == "-" or args.output.name.endswith(".tex"):
            output_type = "tex"
        elif args.output.name.endswith(".pdf"):
            output_type = "pdf"
        else:
            print(
                "Error: Unknown file extension, please specify output "
                "type explicitly",
                file=sys.stderr,
            )
            return 1

    if output_type == "tex":
        args.output.write(tex_code.encode())
    elif output_type == "pdf":
        try:
            tex_compile(source=tex_code, dest=args.output)
        except TeXError as err:
            print(f"TeX compiler returned an error (code: {err.code})")
            print("Output from the compiler:")

            for line in err.message.splitlines():
                print(f"> {line}")

            return 1

    return 0


def draw(args):
    """Run the drawing subcommand with the given arguments."""
    tex_code = generate_tex(args)
    return output(args, tex_code)


def add_args(parser):
    """Add the drawing subcommand to a command-line argument parser."""
    subparser = parser.add_parser("draw", description=__doc__)

    add_arg_input(subparser, "file defining the history to draw")
    add_arg_output(subparser, "file where the resulting drawing will be stored", "wb")

    subparser.add_argument(
        "output_type",
        metavar="TYPE",
        nargs="?",
        choices=("tex", "pdf"),
        help="kind of output to generate (default: guess based on output file \
extension, or 'tex' for output to stdout)",
    )

    subparser.add_argument(
        "--orientation",
        choices=("vertical", "horizontal"),
        default="horizontal",
        help="growing direction of the tree (default: %(default)s)",
    )

    subparser.add_argument(
        "--debug",
        action="store_true",
        help="draw rectangles around computed areas to aid in debugging",
    )

    subparser.set_defaults(func=draw)
