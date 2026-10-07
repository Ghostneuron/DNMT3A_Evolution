import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from internal_start_architecture import overlaps, transcriptional_start


class InternalStartArchitectureTests(unittest.TestCase):
    def test_half_open_overlap(self):
        self.assertTrue(overlaps((10, 20), (15, 25)))
        self.assertFalse(overlaps((10, 20), (20, 25)))

    def test_transcriptional_interval_start(self):
        self.assertEqual(transcriptional_start((10, 20), 1), 10)
        self.assertEqual(transcriptional_start((10, 20), -1), 19)


if __name__ == "__main__":
    unittest.main()
