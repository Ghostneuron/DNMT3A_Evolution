import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import ahocorasick


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from dunnart_neocortex_junction_scan import (
    exact_matches,
    normalized_pair_id,
    transcript_seed_specificity,
)


class DunnartNeocortexJunctionScanTests(unittest.TestCase):
    def test_pair_id_normalization(self):
        self.assertEqual(normalized_pair_id("read42/1 extra"), "read42")
        self.assertEqual(normalized_pair_id("read42 1:N:0"), "read42")

    def test_exact_component_match(self):
        lookup = {
            "AACCGG": [{
                "component": "internal_067_first_junction",
                "model_transcript": "XM_TEST",
                "orientation": "forward",
                "seed": "AACCGG",
            }]
        }
        automaton = ahocorasick.Automaton()
        automaton.add_word("AACCGG", "AACCGG")
        automaton.make_automaton()
        result = exact_matches("TTTAACCGGTT", automaton, lookup)
        self.assertEqual(
            result[(
                "internal_067_first_junction", "XM_TEST", "forward"
            )],
            1,
        )

    def test_seed_specificity(self):
        class Part:
            def __init__(self, sequence):
                self.sequence = sequence

            def extract(self, _):
                return self.sequence

        feature = SimpleNamespace(
            type="mRNA",
            qualifiers={"transcript_id": ["XM_TEST"]},
            location=SimpleNamespace(parts=[Part("TTAACCGGAA")]),
        )
        record = SimpleNamespace(features=[feature], seq="")
        rows = transcript_seed_specificity(record, [{
            "component": "test_junction",
            "model_transcript": "XM_TEST",
            "orientation": "forward",
            "seed": "AACCGG",
        }])
        self.assertEqual(rows[0]["matching_annotated_transcripts"], "XM_TEST")
        self.assertEqual(rows[0]["matching_transcript_count"], 1)


if __name__ == "__main__":
    unittest.main()
