import csv
import gzip
import hashlib
import tempfile
import unittest
from pathlib import Path

from scripts.prepare_mch_cortex_expression import build


class PrepareMchCortexExpressionTests(unittest.TestCase):
    def test_builds_experiment_specific_wt_covariates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data = root / "data"
            data.mkdir()
            gene_info = (
                "#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\n"
                "10090\t1\tGeneA\t-\t-\tEnsembl:ENSMUSG00000000001\n"
            )
            with gzip.open(data / "Mus_musculus.gene_info.gz", "wt") as handle:
                handle.write(gene_info)
            manifest_rows = []
            experiments = [("Dnmt3a1_KO", 3), ("Dnmt3a2_KO", 2), ("Dnmt3a1_delta_N", 3)]
            accession = 0
            for experiment, replicates in experiments:
                for replicate in range(1, replicates + 1):
                    accession += 1
                    filename = f"s{accession}.gz"
                    with gzip.open(data / filename, "wt") as handle:
                        handle.write("ENSMUSG00000000001.4\t10\n")
                    digest = hashlib.sha256((data / filename).read_bytes()).hexdigest()
                    manifest_rows.append([f"s{accession}", "WT", experiment, replicate, filename, digest])
            manifest = root / "manifest.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["accession", "genotype", "experiment", "replicate", "filename", "sha256"])
                writer.writerows(manifest_rows)
            output = root / "out"
            report = build(data, manifest, output)
            self.assertTrue(report["complete"])
            self.assertEqual(report["mapped_unique_gene_symbols"], 1)
            self.assertIn("GeneA", (output / "cortex_WT_expression_covariates.tsv").read_text())


if __name__ == "__main__":
    unittest.main()
