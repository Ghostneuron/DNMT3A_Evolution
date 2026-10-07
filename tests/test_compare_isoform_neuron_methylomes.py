import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.compare_isoform_neuron_methylomes import (
    analyze,
    chromosome_key,
    classify_site,
)


class IsoformNeuronMethylomeTests(unittest.TestCase):
    def test_deposited_lexicographic_chromosome_order(self) -> None:
        self.assertLess(chromosome_key("chr10"), chromosome_key("chr2"))
        self.assertLess(chromosome_key("chr9"), chromosome_key("chrM"))

    def test_site_classification(self) -> None:
        label, a1_ordered, a2_ordered = classify_site(
            (0.8, 0.9), (0.4, 0.5), (0.82, 0.83)
        )
        self.assertEqual(label, "Dnmt3a1_KO_only_loss")
        self.assertTrue(a1_ordered)
        self.assertFalse(a2_ordered)

    def test_six_way_intersection_and_effects(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            raw = Path(temporary)
            groups = ["WT", "WT", "Dnmt3a1_KO", "Dnmt3a1_KO", "Dnmt3a2_KO", "Dnmt3a2_KO"]
            values = [
                [0.8, 0.6, 0.9],
                [0.8, 0.6, 0.9],
                [0.4, 0.2, 0.5],
                [0.4, 0.2, 0.5],
                [0.8, 0.3, 0.9],
                [0.8, 0.3, 0.9],
            ]
            rows = []
            for index, (group, sample_values) in enumerate(zip(groups, values), 1):
                filename = f"sample{index}.bedgraph.gz"
                with gzip.open(raw / filename, "wt") as handle:
                    if index == 6:
                        handle.write("chr1\t5\t6\t0.2\n")
                    for position, value in zip((10, 20, 30), sample_values):
                        handle.write(f"chr1\t{position}\t{position + 1}\t{value}\n")
                rows.append({
                    "sample_accession": f"S{index}",
                    "group": group,
                    "replicate": str((index - 1) % 2 + 1),
                    "local_file": filename,
                })

            result = analyze(rows, raw_dir=raw, window_size=100)
            self.assertEqual(result["common_sites"], 3)
            self.assertAlmostEqual(result["group_means"]["WT"], 2.3 / 3)
            self.assertAlmostEqual(result["group_means"]["Dnmt3a1_KO"], 1.1 / 3)
            self.assertAlmostEqual(result["group_means"]["Dnmt3a2_KO"], 2.0 / 3)
            self.assertEqual(
                result["mean_effect_class_counts"]["Dnmt3a1_KO_only_loss"], 2
            )
            self.assertEqual(result["mean_effect_class_counts"]["shared_loss"], 1)


if __name__ == "__main__":
    unittest.main()
