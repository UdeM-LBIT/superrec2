import unittest
import os
import sys
import tempfile
from superrec2.cli.util import open_std


class TestCliUtil(unittest.TestCase):
    def test_open_std(self):
        with open_std("-", "r") as handle:
            self.assertEqual(handle, sys.stdin)

        self.assertFalse(sys.stdin.closed)

        with open_std("-", "rb") as handle:
            self.assertEqual(handle, sys.stdin.buffer)

        self.assertFalse(sys.stdin.closed)

        with open_std("-", "w") as handle:
            self.assertEqual(handle, sys.stdout)

        self.assertFalse(sys.stdout.closed)

        # The following test does not work because unittest runner replaces sys.stdout
        # with a StringIO instance that does not support accessing a raw buffer
        #
        # with open_std("-", "wb") as handle:
        #    self.assertEqual(handle, sys.stdout.buffer)
        #
        # self.assertFalse(sys.stdout.closed)

        with self.assertRaises(RuntimeError, msg="Invalid mode 'r+' for standard stream"):
            with open_std("-", "r+") as handle:
                pass

        path = tempfile.mkstemp()[1]

        with open_std(path, "w", encoding="utf8") as handle:
            self.assertEqual(handle.name, path)
            self.assertFalse(handle.closed)

        self.assertTrue(handle.closed)
        os.remove(path)
