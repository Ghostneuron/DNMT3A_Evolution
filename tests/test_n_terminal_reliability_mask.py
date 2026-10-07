import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/02_alignment/reliability_mask"


class NTerminalReliabilityMaskTests(unittest.TestCase):
    def setUp(self):
        if not (OUT / "summary.json").exists():
            self.skipTest("N-terminal reliability mask has not been built")

    def test_summary_and_mask_counts(self):
        summary = json.loads((OUT / "summary.json").read_text())
        with (OUT / "n_terminal_mask_manifest.tsv").open() as handle:
            manifest = list(csv.DictReader(handle, delimiter="\t"))
        changed = [row for row in manifest if row["action"] == "mask_primary_to_NNN"]
        self.assertEqual(len(manifest), summary["discordant_taxon_site_cells"])
        self.assertEqual(len(changed), summary["newly_masked_informative_codons"])
        self.assertTrue(all(row["mask"] == "NNN" for row in changed))

    def test_sites_are_n_terminal_and_unique(self):
        with (OUT / "site_reliability.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        positions = [int(row["human_DNMT3A1_aa"]) for row in rows]
        self.assertEqual(positions, sorted(set(positions)))
        self.assertGreaterEqual(min(positions), 1)
        self.assertLessEqual(max(positions), 277)

    def test_disputed_discovery_cells_are_captured(self):
        with (OUT / "n_terminal_mask_manifest.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        observed = {
            (row["safe_id"], row["human_site"])
            for row in rows
            if row["human_site"] in {"P5", "S6", "E18"}
        }
        expected = {
            (species, site)
            for species in {"Puma_concolor", "Vombatus_ursinus"}
            for site in {"P5", "S6", "E18"}
        }
        self.assertTrue(expected <= observed)

    def test_masked_alignment_shape(self):
        def fasta(path):
            records = {}
            current = None
            for line in path.read_text().splitlines():
                if line.startswith(">"):
                    current = line[1:].split()[0]
                    records[current] = ""
                else:
                    records[current] += line.strip()
            return records

        original = fasta(
            ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
        )
        masked = fasta(
            OUT / "DNMT3A_mammals_HyPhy_unique_N_terminal_consensus_masked.fasta"
        )
        self.assertEqual(set(original), set(masked))
        self.assertEqual(
            {identifier: len(sequence) for identifier, sequence in original.items()},
            {identifier: len(sequence) for identifier, sequence in masked.items()},
        )


if __name__ == "__main__":
    unittest.main()
