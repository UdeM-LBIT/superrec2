import unittest
from .subsequences import (
    mask_from_subsequence,
    subsequence_from_mask,
    subsequence_segment_dist as dist,
)


class TestUtilsSubsequences(unittest.TestCase):
    def test_mask_from_subsequence(self):
        self.assertEqual(
            mask_from_subsequence(
                "bcfijl",
                "abcdefghijkl",
            ),
            0b101100100110,
        )

        self.assertEqual(
            mask_from_subsequence(
                "",
                "abcdefghijkl",
            ),
            0,
        )

        self.assertEqual(
            mask_from_subsequence(
                "abcdefghijkl",
                "abcdefghijkl",
            ),
            0b111111111111,
        )

    def test_subsequence_from_mask(self):
        self.assertEqual(
            subsequence_from_mask(
                0b101100100110,
                "abcdefghijkl",
            ),
            list("bcfijl"),
        )

        self.assertEqual(
            subsequence_from_mask(
                0,
                "abcdefghijkl",
            ),
            [],
        )

        self.assertEqual(
            subsequence_from_mask(
                0b111111111111,
                "abcdefghijkl",
            ),
            list("abcdefghijkl"),
        )

    def test_dist(self):
        self.assertEqual(dist(0b11110111, 0b11111111, True), 1)
        self.assertEqual(dist(0b11110111, 0b11111111, False), 1)

        self.assertEqual(dist(0b11100011, 0b11110111, True), 1)
        self.assertEqual(dist(0b11100011, 0b11110111, False), 1)

        self.assertEqual(dist(0b11000010, 0b11100011, True), 2)
        self.assertEqual(dist(0b11000010, 0b11100011, False), 1)

        self.assertEqual(dist(0b01000010, 0b11000010, True), 1)
        self.assertEqual(dist(0b01000010, 0b11000010, False), 0)

        self.assertEqual(dist(0b10101010, 0b01010101, True), -1)
        self.assertEqual(dist(0b10101010, 0b01010101, False), -1)

        self.assertEqual(dist(0b111, 0b110, True), -1)
        self.assertEqual(dist(0b111, 0b110, False), -1)

        mask_all = (1 << 64) - 1
        mask_s1 = mask_all & ~(1 << 52) & ~(1 << 53)

        self.assertEqual(dist(mask_s1, mask_all, True), 1)
        self.assertEqual(dist(mask_s1, mask_all, False), 1)

        mask_s2 = mask_s1 & ~(1 << 63) & ~1

        self.assertEqual(dist(mask_s2, mask_s1, True), 2)
        self.assertEqual(dist(mask_s2, mask_s1, False), 0)
