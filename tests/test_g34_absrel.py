import csv
import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/04_selection/mammals/g34_absrel"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class G34AbsrelTests(unittest.TestCase):
    def test_holm_adjustment(self):
        module = load_module(
            "summarize_g34_absrel",
            ROOT / "scripts/summarize_g34_absrel.py",
        )
        self.assertEqual(module.holm_adjust([0.01, 0.04]), [0.02, 0.04])
        self.assertEqual(module.holm_adjust([0.8, 0.1]), [0.8, 0.2])

    def test_foreground_manifest_and_tree(self):
        manifest = OUTDIR / "foreground_manifest.tsv"
        tree = OUTDIR / "DNMT3A.G34_foreground.nwk"
        if not manifest.exists() or not tree.exists():
            self.skipTest("G34 aBSREL inputs have not been prepared")
        with manifest.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["branch"] for row in rows}, {"Node145", "Node288"})
        self.assertTrue(all(float(row["branch_length"]) > 1e-6 for row in rows))
        self.assertTrue(all(float(row["sh_alrt"]) == 100 for row in rows))
        self.assertTrue(all(float(row["ufboot"]) == 100 for row in rows))
        self.assertTrue(all(row["prank_state_concordant"] == "true" for row in rows))
        self.assertEqual(tree.read_text().count("{Foreground}"), 2)

    def test_result_summary_scope(self):
        summary = OUTDIR / "summary.json"
        if not summary.exists():
            self.skipTest("Targeted G34 aBSREL has not completed")
        data = json.loads(summary.read_text())
        self.assertEqual(data["foreground_branches_tested"], 2)
        self.assertIn("not an independent confirmation", data["interpretation"])


if __name__ == "__main__":
    unittest.main()
