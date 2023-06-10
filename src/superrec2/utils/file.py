"""File management utilities."""
import sys


class _OpenStd:
    def __init__(self, path, mode, args):
        self.path = path
        self.mode = mode
        self.args = args
        self.file = None

    def __enter__(self):
        if self.path == "-":
            if self.mode == "r":
                return sys.stdin
            if self.mode == "rb":
                return sys.stdin.buffer
            if self.mode == "w":
                return sys.stdout
            if self.mode == "wb":
                return sys.stdout.buffer

            raise RuntimeError(f"Invalid mode '{self.mode}' for standard stream")

        self.file = open(
            self.path, self.mode, **self.args
        )  # pylint: disable=unspecified-encoding
        return self.file

    def __exit__(self, kind, value, traceback):
        if self.file is not None:
            self.file.close()
            self.file = None


def open_std(path, mode, **args):
    """
    Open a file in read or write mode with support for special file '-'.

    Opening '-' in read mode will return stdin, or stdout in write mode.
    All other paths behave normally.

    :param path: path of the file to open
    :param mode: mode to use
    :returns: context manager
    """
    return _OpenStd(path, mode, args)
