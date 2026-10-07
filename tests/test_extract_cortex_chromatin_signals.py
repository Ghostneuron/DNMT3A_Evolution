import tempfile
import unittest
from pathlib import Path

import pyBigWig

from scripts.extract_cortex_chromatin_signals import mean_signal


class CortexChromatinSignalTests(unittest.TestCase):
    def test_mean_signal_in_synthetic_bigwig(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test.bw"
            writer = pyBigWig.open(str(path), "w")
            writer.addHeader([("chr1", 1000)])
            writer.addEntries(["chr1"], [100], ends=[200], values=[2.0])
            writer.close()
            reader = pyBigWig.open(str(path))
            try:
                self.assertAlmostEqual(mean_signal(reader, "chr1", 100, 200), 2.0)
                self.assertEqual(mean_signal(reader, "chr1", 300, 400), 0.0)
            finally:
                reader.close()


if __name__ == "__main__":
    unittest.main()
