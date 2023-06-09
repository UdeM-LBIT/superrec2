import unittest
from superrec2.utils.tex import measure, MeasureBox


class TestTex(unittest.TestCase):
    def test_measure(self):
        self.assertEqual(
            measure(
                [
                    "abcdef",
                    r"\(abcdef\)",
                    r"\(a_1b_2c_3d_4e_5f_5\)",
                ]
            ),
            [
                MeasureBox(28.34, 7.05, 0.10999),
                MeasureBox(29.7385, 6.94444, 1.94444),
                MeasureBox(55.57887, 6.94444, 1.94444),
            ],
        )
