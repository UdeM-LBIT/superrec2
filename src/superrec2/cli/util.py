"""Utilities for the command-line interface."""
import argparse


def add_arg_input(parser, message, mode="r"):
    """Add argparse argument for specifying an input file."""
    parser.add_argument(
        "--input",
        metavar="PATH",
        type=argparse.FileType(mode),
        default="-",
        help=f"path to {message} (default: read from stdin)",
    )


def add_arg_output(parser, message, mode="w"):
    """Add argparse argument for specifying an output file."""
    parser.add_argument(
        "--output",
        metavar="PATH",
        type=argparse.FileType(mode),
        default="-",
        help=f"path to {message} (default: write to stdout)",
    )
