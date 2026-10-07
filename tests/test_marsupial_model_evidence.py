import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from marsupial_model_evidence import LONG_READS, RNA_COVERAGE, match_integer


class MarsupialModelEvidenceTests(unittest.TestCase):
    def test_note_count_parsing(self):
        note = "similarity to: 43 long SRA reads, with 98% coverage by RNAseq alignments"
        self.assertEqual(match_integer(LONG_READS, note), 43)
        self.assertEqual(match_integer(RNA_COVERAGE, note), 98)


if __name__ == "__main__":
    unittest.main()
