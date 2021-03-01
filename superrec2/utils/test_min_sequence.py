import unittest
from infinity import inf
from .min_sequence import MinSequence


class TestUtilsMinSequence(unittest.TestCase):
    def test_empty(self):
        seq = MinSequence()
        self.assertEqual(seq.min, inf)
        self.assertEqual(list(seq), [])

    def test_replace(self):
        seq = MinSequence()

        seq.update((100, 'hundred'))
        self.assertEqual(seq.min, 100)
        self.assertEqual(list(seq), ['hundred'])

        seq.update((50, 'fifty'))
        self.assertEqual(seq.min, 50)
        self.assertEqual(list(seq), ['fifty'])

    def test_accumulate(self):
        seq = MinSequence()

        seq.update((100, 'hundred'))
        seq.update((50, 'fifty'))
        seq.update((50, 'second'))
        seq.update((50, 'third'))
        self.assertEqual(seq.min, 50)
        self.assertEqual(list(seq), ['fifty', 'second', 'third'])

        seq.update((25, 'twenty-five'))
        self.assertEqual(seq.min, 25)
        self.assertEqual(list(seq), ['twenty-five'])

    def test_accumulate_max(self):
        seq = MinSequence(max_keep=1)

        seq.update((75, 'twenty-five'))
        self.assertEqual(seq.min, 75)
        self.assertEqual(list(seq), ['twenty-five'])

        seq.update((100, 'hundred'))
        seq.update((50, 'fifty'))
        seq.update((50, 'second'))
        seq.update((50, 'third'))
        self.assertEqual(seq.min, 50)
        self.assertEqual(list(seq), ['fifty'])
