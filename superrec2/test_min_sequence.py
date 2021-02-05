import unittest
from .min_sequence import MinSequence

class TestMinSequence(unittest.TestCase):
    def test_empty(self):
        seq = MinSequence()
        self.assertEqual(seq.min, float('inf'))
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

        print(seq._items)
        seq.update((100, 'hundred'))
        print(seq._items)
        seq.update((50, 'fifty'))
        print(seq._items)
        seq.update((50, 'second'))
        print(seq._items)
        seq.update((50, 'third'))
        print(seq._items)
        self.assertEqual(seq.min, 50)
        self.assertEqual(list(seq), ['fifty', 'second', 'third'])

        seq.update((25, 'twenty-five'))
        self.assertEqual(seq.min, 25)
        self.assertEqual(list(seq), ['twenty-five'])
