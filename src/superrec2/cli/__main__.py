import sys
import argparse
from . import reconcile, draw


def run():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(required=True)

    reconcile.add_args(subparsers)
    draw.add_args(subparsers)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(run())
