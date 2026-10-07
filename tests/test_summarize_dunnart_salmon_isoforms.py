import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from summarize_dunnart_salmon_isoforms import aggregate_families, family_for_model


class SummarizeDunnartSalmonIsoformsTests(unittest.TestCase):
    def test_family_assignment(self):
        self.assertEqual(
            family_for_model("XM_074300067.1"),
            "internal_066_067_family",
        )

    def test_family_aggregation(self):
        rows = [
            {"Name": "XM_074300066.1", "NumReads": "2.5", "TPM": "10"},
            {"Name": "XM_074300067.1", "NumReads": "3.5", "TPM": "20"},
        ]
        families = aggregate_families(rows)
        self.assertEqual(
            families["internal_066_067_family"]["NumReads"],
            6.0,
        )
        self.assertEqual(
            families["internal_066_067_family"]["TPM"],
            30.0,
        )


if __name__ == "__main__":
    unittest.main()
