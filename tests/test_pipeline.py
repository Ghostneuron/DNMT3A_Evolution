#!/usr/bin/env python3
"""Offline smoke tests for curation and MACSE-cleaning logic."""

from __future__ import annotations

import subprocess
import shutil
import os
import tempfile
import textwrap
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/dnmt3a_pipeline.py"
DOMAINS = Path(__file__).resolve().parents[1] / "config/human_domains.tsv"


class PipelineSmokeTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "config").mkdir()
        (self.root / "data/raw/ncbi_dataset/data").mkdir(parents=True)
        (self.root / "config/human_domains.tsv").write_text(DOMAINS.read_text())
        self.config = self.root / "config/project.toml"
        self.config.write_text(textwrap.dedent("""
            [project]
            gene = "DNMT3A"
            taxon_scope = "test"
            human_species = "Homo sapiens"
            human_protein_accession = "NP_TEST.1"
            human_length_aa = 12
            minimum_length_fraction = 0.90
            maximum_length_fraction = 1.10
            minimum_clean_codon_occupancy = 0.66
            [paths]
            protein_fasta = "data/raw/ncbi_dataset/data/protein.faa"
            cds_fasta = "data/raw/ncbi_dataset/data/cds.fna"
            curated_dir = "results/01_curated"
            alignment_dir = "results/02_alignment"
            tree_dir = "results/03_tree"
            selection_dir = "results/04_selection"
            [tools]
            macse = "macse"
            mafft = "mafft"
            iqtree = "iqtree2"
            hyphy = "hyphy"
            [alignment]
            macse_refine_iterations = 2
            macse_optimization = 1
            macse_local_realign_initial_fraction = 0.10
            [tree]
            sequence_type = "DNA"
            model = "GTR+F+G4"
            bootstrap_replicates = 1000
            alrt_replicates = 1000
            threads = "AUTO"
            mammal_constraint = "results/03_tree/constraints/NCBI_mammal_orders_constraint.nwk"
            [selection]
            run_fubar = true
            run_meme = true
            run_absrel = true
            meme_candidate_sites = [9, 30, 107]
        """))
        seqs = {
            "Homo sapiens": ("NP_TEST.1", "ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"),
            "Mus musculus": ("NP_MOUSE.1", "ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCC"),
            "Canis lupus": ("XP_DOG.1", "ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCA"),
        }
        proteins, cds = [], []
        for species, (accession, nt) in seqs.items():
            aa = "M" + "A" * 11
            proteins.append(f">{accession} [gene=DNMT3A] [organism={species}] [isoform=1]\n{aa}\n")
            cds.append(f">CDS_{accession} [gene=DNMT3A] [organism={species}]\n{nt}\n")
        raw = self.root / "data/raw/ncbi_dataset/data"
        (raw / "protein.faa").write_text("".join(proteins))
        (raw / "cds.fna").write_text("".join(cds))

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def run_stage(self, stage: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["python3", str(SCRIPT), "--config", str(self.config), stage],
            text=True, capture_output=True,
        )

    def test_curate_and_clean(self) -> None:
        result = self.run_stage("curate")
        self.assertEqual(result.returncode, 0, result.stderr)
        curated = self.root / "results/01_curated/DNMT3A_representative_cds.fna"
        self.assertTrue(curated.exists())
        text = curated.read_text()
        self.assertEqual(text.count(">"), 3)

        alignment_dir = self.root / "results/02_alignment"
        alignment_dir.mkdir(parents=True)
        (alignment_dir / "DNMT3A_MACSE_NT.fasta").write_text(text)
        result = self.run_stage("clean")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((alignment_dir / "DNMT3A_MACSE_clean_HyPhy.fasta").exists())
        sequence_qc = (alignment_dir / "sequence_codon_qc.tsv").read_text()
        self.assertIn("retained_fraction_after_column_filter", sequence_qc)
        crosswalk = (alignment_dir / "human_coordinate_crosswalk.tsv").read_text()
        self.assertIn("human_DNMT3A1_aa", crosswalk)
        self.assertIn("\t12\t", crosswalk)

        result = subprocess.run(
            [
                "python3", str(SCRIPT), "--config", str(self.config),
                "--dataset", "mammals", "--dry-run", "tree",
            ],
            text=True, capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("DNMT3A_mammals.fasta", result.stderr)

    @unittest.skipUnless(
        os.environ.get("RUN_EXTERNAL_ALIGNMENT_TEST") == "1" and shutil.which("macse") and shutil.which("mafft"),
        "set RUN_EXTERNAL_ALIGNMENT_TEST=1 with MACSE and MAFFT installed",
    )
    def test_external_alignment_command(self) -> None:
        result = self.run_stage("curate")
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.run_stage("align")
        log = self.root / "results/02_alignment/alignment_commands.log"
        details = log.read_text() if log.exists() else "no alignment log"
        self.assertEqual(result.returncode, 0, result.stderr + "\n" + details)
        alignment = self.root / "results/02_alignment/DNMT3A_MACSE_NT.fasta"
        self.assertTrue(alignment.exists())
        self.assertGreater(alignment.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
