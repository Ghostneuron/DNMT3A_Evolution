import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from prepare_dunnart_salmon_reference import (
    automaton_shared_kmer_count,
    kmer_automaton,
    kmers,
    shared_kmer_count,
)


class PrepareDunnartSalmonReferenceTests(unittest.TestCase):
    def test_kmers_exclude_ambiguous_windows(self):
        self.assertEqual(kmers("AACNGT", length=3), {"AAC"})

    def test_shared_kmer_count_is_distinct(self):
        target = kmers("AACCGG", length=3)
        self.assertEqual(shared_kmer_count("AACCAACC", target, length=3), 2)
        self.assertEqual(
            automaton_shared_kmer_count(
                "AACCAACC",
                kmer_automaton(target),
            ),
            2,
        )


if __name__ == "__main__":
    unittest.main()
