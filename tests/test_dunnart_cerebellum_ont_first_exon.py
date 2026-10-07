import sys
import unittest
from pathlib import Path

import ahocorasick


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from dunnart_cerebellum_ont_first_exon import exact_matches


class DunnartCerebellumOntFirstExonTests(unittest.TestCase):
    def test_exact_matches_preserves_target_and_component(self):
        rows = {
            "AACCGG": [{
                "target_transcript": "XM_TEST",
                "component": "first_junction",
                "orientation": "forward",
                "seed": "AACCGG",
            }]
        }
        automaton = ahocorasick.Automaton()
        automaton.add_word("AACCGG", "AACCGG")
        automaton.make_automaton()
        result = exact_matches("TTTAACCGGTTT", automaton, rows)
        self.assertEqual(
            result[("XM_TEST", "first_junction", "forward")], 1
        )


if __name__ == "__main__":
    unittest.main()
