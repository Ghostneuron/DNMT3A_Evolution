import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from direct_tss_evidence import evidence_tier, point_interval_distance


class DirectTssEvidenceTests(unittest.TestCase):
    def test_half_open_interval_distance(self):
        self.assertEqual(point_interval_distance(10, 10, 20), 0)
        self.assertEqual(point_interval_distance(19, 10, 20), 0)
        self.assertEqual(point_interval_distance(20, 10, 20), 1)
        self.assertEqual(point_interval_distance(8, 10, 20), 2)

    def test_evidence_tiers(self):
        self.assertEqual(
            evidence_tier(0, 0), "direct_CAGE_peak_exact_or_near_exact"
        )
        self.assertEqual(
            evidence_tier(2, 7), "direct_CAGE_peak_near_annotated_TSS"
        )
        self.assertEqual(
            evidence_tier(100, 100), "direct_CAGE_peak_within_500bp"
        )


if __name__ == "__main__":
    unittest.main()
