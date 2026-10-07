import csv
import gzip
import hashlib
import tempfile
import unittest
from pathlib import Path

import pyBigWig

from scripts.extract_p21_polycomb_signals import build


class P21PolycombSignalTests(unittest.TestCase):
    def test_build_summarizes_six_validated_tracks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genes = root / "genes.tsv.gz"
            with gzip.open(genes, "wt", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["chrom", "start", "end", "gene", "strand"], delimiter="\t")
                writer.writeheader(); writer.writerow({"chrom": "chr1", "start": 100, "end": 200, "gene": "GeneA", "strand": "+"})
            labels = ["wt_dnmt3a", "delta_n_dnmt3a", "input", "h3k4me3", "h3k27me3", "h2ak119ub"]
            rows = []
            for index, label in enumerate(labels):
                path = root / f"{label}.bw"
                bw = pyBigWig.open(str(path), "w"); bw.addHeader([("chr1", 1000)])
                bw.addEntries(["chr1"], [0], ends=[1000], values=[float(index + 1)]); bw.close()
                rows.append({"accession": f"GSM{index}", "label": label, "sample_scope": "test", "genotype": "test", "filename": path.name, "url": "test", "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
            manifest = root / "manifest.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t"); writer.writeheader(); writer.writerows(rows)
            report = build(genes, manifest, root, root / "out", 20)
            self.assertTrue(report["complete"]); self.assertEqual(report["genes"], 1)
            with gzip.open(root / "out/p21_polycomb_gene_signals.tsv.gz", "rt") as handle:
                result = next(csv.DictReader(handle, delimiter="\t"))
            self.assertAlmostEqual(float(result["h2ak119ub_gene_body_mean"]), 6.0)


if __name__ == "__main__":
    unittest.main()
