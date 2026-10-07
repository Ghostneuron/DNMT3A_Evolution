import unittest
from pathlib import Path

from scripts.run_mch_bismark import (
    align_command,
    extract_command,
    prepare_command,
    read_samples,
    reference_prepared,
    trim_command,
)


class MchBismarkPipelineTests(unittest.TestCase):
    def setUp(self):
        self.external = Path("/external/mch")
        self.samples = read_samples()

    def test_manifest_has_eight_paired_samples(self):
        self.assertEqual(len(self.samples), 8)
        for sample in self.samples.values():
            self.assertIn("read1", sample)
            self.assertIn("read2", sample)
        self.assertEqual(
            self.samples["SRR16893001"]["genotype"], "Dnmt3a1_delta_N"
        )

    def test_prepare_uses_composite_mm9_and_bowtie2(self):
        command = prepare_command(self.external, 4)
        self.assertIn("--bowtie2", command)
        self.assertEqual(command[-1], "/external/mch/reference/mm9_controls")

    def test_incomplete_reference_is_not_prepared(self):
        self.assertFalse(reference_prepared(self.external))

    def test_pilot_alignment_is_paired_and_limited(self):
        sample = self.samples["SRR16892995"]
        command = align_command(self.external, sample, 4, 100_000)
        self.assertIn("--upto", command)
        self.assertIn("100000", command)
        self.assertIn("-1", command)
        self.assertIn("-2", command)
        self.assertNotIn("--non_directional", command)
        # Bismark Rust 3.1.0 parses read-group flags but rejects them in v1.
        self.assertNotIn("--rg_tag", command)
        self.assertEqual(
            command[command.index("--temp_dir") + 1],
            "/external/mch/results/bismark/SRR16892995/alignment/tmp",
        )
        self.assertIn("trimmed.fastq.gz", command[-1])

    def test_fastp_only_trims_adapters_for_pilot(self):
        sample = self.samples["SRR16892995"]
        command = trim_command(self.external, sample, 4, 100_000)
        self.assertIn("--detect_adapter_for_pe", command)
        self.assertIn("--disable_quality_filtering", command)
        self.assertIn("--disable_trim_poly_g", command)
        self.assertEqual(command[command.index("--reads_to_process") + 1], "100000")

    def test_extraction_retains_contexts_and_overlap_safety(self):
        sample = self.samples["SRR16892995"]
        command = extract_command(self.external, sample, 4, (1, 2, 3, 4))
        self.assertIn("--CX", command)
        self.assertIn("--comprehensive", command)
        self.assertIn("--no_overlap", command)
        self.assertNotIn("--merge_non_CpG", command)
        self.assertNotIn("filter", command)
        self.assertEqual(command[command.index("--ignore_r2") + 1], "2")


if __name__ == "__main__":
    unittest.main()
