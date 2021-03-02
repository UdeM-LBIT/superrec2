import unittest
from .tex import measure, MeasureBox


class TestTex(unittest.TestCase):
    def test_measure(self):
        self.assertEqual(
            measure([
                "abcdef",
                r"\(abcdef\)",
                r"\(a_1b_2c_3d_4e_5f_5\)",
            ]),
            [
                MeasureBox("28.34pt", "7.05pt", "0.10999pt"),
                MeasureBox("29.7385pt", "6.94444pt", "1.94444pt"),
                MeasureBox("55.57887pt", "6.94444pt", "1.94444pt"),
            ]
        )
