import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from promoter_sequence_similarity import kmers


class PromoterSequenceSimilarityTests(unittest.TestCase):
    def test_kmers_are_unique_and_exclude_ambiguous(self):
        self.assertEqual(kmers("AAAA", 2), {"AA"})
        self.assertEqual(kmers("AANAA", 2), {"AA"})


if __name__ == "__main__":
    unittest.main()
