import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = (
    ROOT
    / "results/04_selection/mammals/full_meme"
    / "consensus_masked"
)


class FullMemeConsensusMaskTests(unittest.TestCase):
    def test_mask_manifest(self):
        path = OUTDIR / "mask_manifest.tsv"
        if not path.exists():
            self.skipTest("Consensus masks have not been prepared")
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 6)
        self.assertEqual(
            {row["safe_id"] for row in rows},
            {"Puma_concolor", "Vombatus_ursinus"},
        )
        self.assertTrue(all(row["mask"] == "NNN" for row in rows))

    def test_mask_summary(self):
        path = OUTDIR / "mask_summary.json"
        if not path.exists():
            self.skipTest("Consensus masks have not been prepared")
        rows = json.loads(path.read_text())
        self.assertEqual({row["human_site"] for row in rows}, {"P5", "S6", "E18"})
        self.assertTrue(all(row["masked_taxa_count"] == 2 for row in rows))
        self.assertTrue(all(row["unmasked_taxa_count"] == 248 for row in rows))


if __name__ == "__main__":
    unittest.main()
