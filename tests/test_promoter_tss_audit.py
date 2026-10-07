import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from promoter_tss_audit import downstream_distance, genomic_position


class PromoterTssAuditTests(unittest.TestCase):
    def test_genomic_coordinate_conversion(self):
        self.assertEqual(genomic_position(1001, 0), 1001)
        self.assertEqual(genomic_position(1001, 99), 1100)

    def test_downstream_distance_plus_strand(self):
        self.assertEqual(downstream_distance(1200, 1000, 1), 200)

    def test_downstream_distance_minus_strand(self):
        self.assertEqual(downstream_distance(800, 1000, -1), 200)


if __name__ == "__main__":
    unittest.main()
