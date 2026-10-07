import tempfile
import unittest
from pathlib import Path

from scripts.prepare_mch_reference import build_composite, fasta_headers


class PrepareMchReferenceTests(unittest.TestCase):
    def test_builds_mouse_control_composite(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mouse = root / "mouse.fa"
            controls = root / "controls.fa"
            output = root / "composite.fa"
            mouse.write_text(">chr1\nACGT\n")
            controls.write_text(">phage_lambda\nCCAA\n>plasmid_puc19c\nCGCG\n")
            build_composite(mouse, controls, output)
            self.assertEqual(
                fasta_headers(output), ["chr1", "phage_lambda", "plasmid_puc19c"]
            )

    def test_rejects_missing_controls(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mouse = root / "mouse.fa"
            controls = root / "controls.fa"
            mouse.write_text(">chr1\nACGT\n")
            controls.write_text(">lambda\nCCAA\n")
            with self.assertRaises(ValueError):
                build_composite(mouse, controls, root / "out.fa")


if __name__ == "__main__":
    unittest.main()
