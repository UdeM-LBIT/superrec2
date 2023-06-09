import unittest
from superrec2.utils.subsequences import (
    subseq_complete,
    mask_from_subseq,
    subseq_from_mask,
    subseq_segment_dist as dist,
)


class TestUtilsSubsequences(unittest.TestCase):
    def test_subseq_complete(self):
        self.assertEqual(subseq_complete(""), 0)
        self.assertEqual(subseq_complete("a"), 1)
        self.assertEqual(subseq_complete("1234"), 0b1111)

    def test_mask_from_subseq(self):
        self.assertEqual(
            mask_from_subseq(
                "bcfijl",
                "abcdefghijkl",
            ),
            0b1011_0010_0110,
        )

        self.assertEqual(
            mask_from_subseq(
                "",
                "abcdefghijkl",
            ),
            0,
        )

        self.assertEqual(
            mask_from_subseq(
                "abcdefghijkl",
                "abcdefghijkl",
            ),
            0b1111_1111_1111,
        )

    def test_subseq_from_mask(self):
        self.assertEqual(
            subseq_from_mask(
                0b1011_0010_0110,
                "abcdefghijkl",
            ),
            list("bcfijl"),
        )

        self.assertEqual(
            subseq_from_mask(
                0,
                "abcdefghijkl",
            ),
            [],
        )

        self.assertEqual(
            subseq_from_mask(
                0b1111_1111_1111,
                "abcdefghijkl",
            ),
            list("abcdefghijkl"),
        )

    def test_subseq_segment_dist(self):
        self.assertEqual(dist(0b1111_0111, 0b1111_1111, True), 1)
        self.assertEqual(dist(0b1111_0111, 0b1111_1111, False), 1)

        self.assertEqual(dist(0b1110_0011, 0b1111_0111, True), 1)
        self.assertEqual(dist(0b1110_0011, 0b1111_0111, False), 1)

        self.assertEqual(dist(0b1100_0010, 0b1110_0011, True), 2)
        self.assertEqual(dist(0b1100_0010, 0b1110_0011, False), 1)

        self.assertEqual(dist(0b0100_0010, 0b1100_0010, True), 1)
        self.assertEqual(dist(0b0100_0010, 0b1100_0010, False), 0)

        self.assertEqual(dist(0b1010_1010, 0b0101_0101, True), -1)
        self.assertEqual(dist(0b1010_1010, 0b0101_0101, False), -1)

        self.assertEqual(dist(0b111, 0b110, True), -1)
        self.assertEqual(dist(0b111, 0b110, False), -1)

        mask_all = (1 << 64) - 1
        mask_s1 = mask_all & ~(1 << 52) & ~(1 << 53)

        self.assertEqual(dist(mask_s1, mask_all, True), 1)
        self.assertEqual(dist(mask_s1, mask_all, False), 1)

        mask_s2 = mask_s1 & ~(1 << 63) & ~1

        self.assertEqual(dist(mask_s2, mask_s1, True), 2)
        self.assertEqual(dist(mask_s2, mask_s1, False), 0)
