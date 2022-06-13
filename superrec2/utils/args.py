"""Utilities for command-line arguments."""

def add_arg_input(parser, message):
    """Add argparse argument for specifying an input file."""
    parser.add_argument(
        "--input",
        metavar="PATH",
        default="-",
        help=f"path to {message} (default: read from stdin)",
    )


def add_arg_output(parser, message):
    """Add argparse argument for specifying an output file."""
    parser.add_argument(
        "--output",
        metavar="PATH",
        default="-",
        help=f"path to {message} (default: write to stdout)",
    )
