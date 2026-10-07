import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/04_selection/mammals/full_meme"


def load_module():
    path = ROOT / "scripts/full_meme_sensitivity.py"
    spec = importlib.util.spec_from_file_location("full_meme_sensitivity", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class FullMemeSensitivityTests(unittest.TestCase):
    def test_holm_adjustment(self):
        module = load_module()
        observed = module.holm_adjust([0.01, 0.04, 0.03])
        expected = [0.03, 0.06, 0.06]
        for left, right in zip(observed, expected):
            self.assertAlmostEqual(left, right)

    def test_classification(self):
        module = load_module()
        self.assertIn("all_three", module.classify(0.01, 0.02))
        self.assertIn("prank_sensitive", module.classify(0.01, 0.2))
        self.assertIn("fails_mafft_and_prank", module.classify(0.2, 0.2))

    def test_completed_summary(self):
        path = OUTDIR / "discovery_alignment_sensitivity_summary.json"
        if not path.exists():
            self.skipTest("Discovery sensitivity analysis has not completed")
        data = json.loads(path.read_text())
        self.assertEqual(data["discovery_sites"], 3)
        self.assertEqual(data["meme_supported_all_three_alignments"], 0)


if __name__ == "__main__":
    unittest.main()
