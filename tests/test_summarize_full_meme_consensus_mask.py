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


class SummarizeFullMemeConsensusMaskTests(unittest.TestCase):
    def test_summary(self):
        path = OUTDIR / "consensus_masked_summary.json"
        if not path.exists():
            self.skipTest("Consensus-masked MEME summary has not been built")
        data = json.loads(path.read_text())
        self.assertEqual(data["sites_tested"], 3)
        self.assertEqual(data["sites_supported_after_consensus_masking_p05"], 0)
        self.assertEqual(data["masked_taxa_per_site"], 2)

    def test_all_rows_are_nonsignificant(self):
        path = OUTDIR / "consensus_masked_summary.tsv"
        if not path.exists():
            self.skipTest("Consensus-masked MEME summary has not been built")
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["human_site"] for row in rows}, {"P5", "S6", "E18"})
        self.assertTrue(
            all(row["supported_after_consensus_masking_p05"] == "false" for row in rows)
        )


if __name__ == "__main__":
    unittest.main()
