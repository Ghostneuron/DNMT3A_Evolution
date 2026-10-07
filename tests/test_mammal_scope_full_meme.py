import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = (
    ROOT
    / "results/04_selection/mammals/mammal_scope/full_meme/summary.json"
)


class MammalScopeFullMemeTests(unittest.TestCase):
    def setUp(self):
        if not SUMMARY.exists():
            self.skipTest("Corrected mammal-scope full MEME scan is incomplete")

    def test_complete_correction_family(self):
        summary = json.loads(SUMMARY.read_text())
        self.assertEqual(summary["sequences"], 250)
        self.assertEqual(summary["codon_sites_tested"], 909)
        self.assertEqual(
            set(summary["retained_candidate_results"]), {"G34", "T12", "S97"}
        )


if __name__ == "__main__":
    unittest.main()
