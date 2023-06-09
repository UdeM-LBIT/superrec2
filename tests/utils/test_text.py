import unittest
from superrec2.utils.text import balanced_wrap


class TestText(unittest.TestCase):
    def test_balanced_wrap(self):
        self.assertEqual(
            balanced_wrap("this is a short test with short words", 32),
            "this is a short test\nwith short words",
        )
        self.assertEqual(
            balanced_wrap("elongated beautiful unbreakable wording", 8),
            "elongated\nbeautiful\nunbreakable\nwording",
        )
        self.assertEqual(
            balanced_wrap("", 8),
            "",
        )
