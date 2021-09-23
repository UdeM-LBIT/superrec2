import unittest
from .escape import escape, unescape


class TestUtilsEscape(unittest.TestCase):
    def test_escape(self):
        escape_chars = r"\,()"
        self.assertEqual(escape("cas9", escape_chars), "cas9")
        self.assertEqual(escape("cas3'", escape_chars), "cas3'")
        self.assertEqual(escape('cas3"', escape_chars), 'cas3"')
        self.assertEqual(escape("gene(parts)", escape_chars), r"gene\(parts\)")
        self.assertEqual(escape(r"\\ba\ck\\", escape_chars), r"\\\\ba\\ck\\\\")

    def test_unescape(self):
        self.assertEqual(unescape("cas9"), "cas9")
        self.assertEqual(unescape("cas3'"), "cas3'")
        self.assertEqual(unescape('cas3"'), 'cas3"')
        self.assertEqual(unescape(r"gene\(parts\)"), "gene(parts)")
        self.assertEqual(unescape(r"\a\b\c\d\)"), "abcd)")
