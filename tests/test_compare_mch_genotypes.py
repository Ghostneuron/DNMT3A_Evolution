import csv
import tempfile
import unittest
from pathlib import Path

from scripts.compare_mch_genotypes import compare


def write_summary(path: Path, ca_fraction: float) -> None:
    methylated = round(ca_fraction * 1000)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["scope", "context", "methylated", "unmethylated", "sites", "fraction"]
        )
        writer.writerow(
            ["autosome", "CA", methylated, 1000 - methylated, 500, ca_fraction]
        )
        writer.writerow(["phage_lambda", "CA", 2, 998, 400, 0.002])


class CompareMchGenotypesTests(unittest.TestCase):
    def test_reports_genotype_effect_relative_to_wt(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            summaries = {}
            values = {
                "SRR16892995": 0.020,
                "SRR16892996": 0.022,
                "SRR16892997": 0.010,
                "SRR16892998": 0.012,
            }
            for run, value in values.items():
                path = root / f"{run}.tsv"
                write_summary(path, value)
                summaries[run] = path
            output = root / "out"
            report = compare(root, summaries, output)
            self.assertFalse(report["complete"])
            contrasts = (output / "genotype_vs_wt_contrasts.tsv").read_text()
            self.assertIn("Dnmt3a1_KO-WT", contrasts)
            self.assertIn("-0.01002004008", contrasts)


if __name__ == "__main__":
    unittest.main()
