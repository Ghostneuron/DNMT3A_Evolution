import hashlib
import tempfile
import unittest
from pathlib import Path

from scripts.download_mch_fastq import aria2_input_lines, md5, validate_completed


class DownloadMchFastqTests(unittest.TestCase):
    def test_aria2_input_includes_expected_checksum(self) -> None:
        rows = [{
            "fastq_url": "https://example.org/read.fastq.gz",
            "md5": "0123456789abcdef0123456789abcdef",
        }]
        self.assertEqual(
            aria2_input_lines(rows),
            [
                "https://example.org/read.fastq.gz",
                "  out=read.fastq.gz",
                "  checksum=md5=0123456789abcdef0123456789abcdef",
            ],
        )

    def test_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            content = b"example"
            (output / "read.fastq.gz").write_bytes(content)
            rows = [{
                "fastq_url": "https://example.org/read.fastq.gz",
                "bytes": str(len(content)),
                "md5": hashlib.md5(content).hexdigest(),
                "run_accession": "SRR1",
                "read_pair": "1",
            }]
            result = validate_completed(rows, output)
        self.assertTrue(result[0]["size_matches"])
        self.assertTrue(result[0]["md5_matches"])


if __name__ == "__main__":
    unittest.main()
