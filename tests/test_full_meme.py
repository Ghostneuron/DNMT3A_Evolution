import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/04_selection/mammals/full_meme"


def load_module():
    path = ROOT / "scripts/summarize_full_meme.py"
    spec = importlib.util.spec_from_file_location("summarize_full_meme", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class FullMemeTests(unittest.TestCase):
    def test_bh_adjustment(self):
        module = load_module()
        observed = module.bh_adjust([0.01, 0.04, 0.03, 0.002])
        expected = [0.02, 0.04, 0.04, 0.008]
        for left, right in zip(observed, expected):
            self.assertAlmostEqual(left, right)

    def test_completed_scan_summary(self):
        path = OUTDIR / "summary.json"
        if not path.exists():
            self.skipTest("Alignment-wide MEME scan has not completed")
        data = json.loads(path.read_text())
        self.assertEqual(data["codon_sites_tested"], 901)
        self.assertEqual(data["sequences"], 250)
        self.assertIn("G34", data["targeted_candidate_full_scan_results"])


if __name__ == "__main__":
    unittest.main()
