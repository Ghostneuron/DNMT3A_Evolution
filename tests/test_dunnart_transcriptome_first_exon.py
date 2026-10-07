import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from dunnart_transcriptome_first_exon import (
    SEED_LENGTH,
    junction_seeds,
    sequence_entropy,
    spaced_seeds,
)


class DunnartTranscriptomeFirstExonTests(unittest.TestCase):
    def test_entropy(self):
        self.assertEqual(sequence_entropy("A" * 35), 0.0)
        self.assertGreater(sequence_entropy("ACGT" * 9), 1.9)

    def test_spaced_seed_length(self):
        sequence = "ACGT" * 60
        seeds = spaced_seeds(sequence)
        self.assertTrue(seeds)
        self.assertTrue(all(len(seed) == SEED_LENGTH for _, seed in seeds))

    def test_junction_seeds_cross_boundary(self):
        left = "ACGT" * 30
        right = "TGCA" * 30
        for start, seed in junction_seeds(left, right):
            self.assertLess(start, len(left))
            self.assertGreater(start + len(seed), len(left))


if __name__ == "__main__":
    unittest.main()
