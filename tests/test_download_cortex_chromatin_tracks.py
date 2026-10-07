import hashlib
import tempfile
import unittest
from pathlib import Path

from scripts.download_cortex_chromatin_tracks import validate


class CortexChromatinDownloadTests(unittest.TestCase):
    def test_validate_checks_size_and_sha256(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            content = b"bigwig-placeholder"
            (root / "track.bw").write_bytes(content)
            rows = [{
                "accession": "GSM1", "filename": "track.bw",
                "bytes": str(len(content)), "sha256": hashlib.sha256(content).hexdigest(),
            }]
            result = validate(rows, root)[0]
            self.assertTrue(result["size_matches"])
            self.assertTrue(result["sha256_matches"])


if __name__ == "__main__":
    unittest.main()
