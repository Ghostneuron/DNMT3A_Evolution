import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from dunnart_cerebellum_isoform_classification import classify_isoform


class DunnartCerebellumIsoformClassificationTests(unittest.TestCase):
    def test_isoform_evidence_precedence(self):
        self.assertEqual(
            classify_isoform(2, 0, 0), "full_length_XM_074300059"
        )
        self.assertEqual(
            classify_isoform(2, 1, 0), "internal_XM_074300067"
        )
        self.assertEqual(
            classify_isoform(0, 0, 0), "downstream_DNMT3A_only"
        )


if __name__ == "__main__":
    unittest.main()
