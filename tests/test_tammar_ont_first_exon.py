import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from tammar_ont_first_exon import classify_read


class TammarOntFirstExonTests(unittest.TestCase):
    def test_classify_read_counts_components_and_orientation(self):
        seeds = [
            {"component": "first_exon", "orientation": "forward", "seed": "AACCGG"},
            {"component": "first_junction", "orientation": "forward", "seed": "CCGGTT"},
            {"component": "first_exon", "orientation": "reverse_complement", "seed": "TTGGCC"},
        ]
        result = classify_read("AAAACCGGAAAACCGGTT", seeds)
        self.assertEqual(result[("first_exon", "forward")], 1)
        self.assertEqual(result[("first_junction", "forward")], 1)
        self.assertNotIn(("first_exon", "reverse_complement"), result)


if __name__ == "__main__":
    unittest.main()
