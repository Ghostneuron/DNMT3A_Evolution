import sys
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from inventory_marsupial_sra import classify_evidence, parse_package


class MarsupialSraInventoryTests(unittest.TestCase):
    def test_long_read_brain_is_first_exon_capable_not_exact_tss(self):
        evidence_class, rank, limit = classify_evidence(
            "RNA-Seq", "cDNA", "OXFORD_NANOPORE", "adult cerebellum"
        )
        self.assertEqual(evidence_class, "full_length_transcript_capable")
        self.assertEqual(rank, 2)
        self.assertIn("not exact capped TSS", limit)

    def test_cage_is_exact_5prime_capable(self):
        evidence_class, rank, _ = classify_evidence(
            "OTHER", "CAGE", "ILLUMINA", "brain CAGE"
        )
        self.assertEqual(evidence_class, "direct_5prime_capable")
        self.assertEqual(rank, 1)

    def test_minimal_package_parsing(self):
        package = ET.fromstring(
            """
            <EXPERIMENT_PACKAGE>
              <EXPERIMENT accession="ERX1">
                <STUDY_REF accession="ERP1"><IDENTIFIERS>
                  <EXTERNAL_ID namespace="BioProject">PRJ1</EXTERNAL_ID>
                </IDENTIFIERS></STUDY_REF>
                <DESIGN><LIBRARY_DESCRIPTOR>
                  <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
                  <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
                  <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
                  <LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT>
                </LIBRARY_DESCRIPTOR></DESIGN>
                <PLATFORM><OXFORD_NANOPORE>
                  <INSTRUMENT_MODEL>PromethION</INSTRUMENT_MODEL>
                </OXFORD_NANOPORE></PLATFORM>
              </EXPERIMENT>
              <STUDY><DESCRIPTOR><STUDY_TITLE>brain atlas</STUDY_TITLE>
              </DESCRIPTOR></STUDY>
              <SAMPLE><TITLE>adult cerebellum</TITLE>
                <SAMPLE_NAME><SCIENTIFIC_NAME>Test species</SCIENTIFIC_NAME></SAMPLE_NAME>
                <SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE>
                  <TAG>tissue_type</TAG><VALUE>Cerebellum</VALUE>
                </SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES>
              </SAMPLE>
              <RUN_SET bases="1234" spots="12"><RUN accession="ERR1"/></RUN_SET>
            </EXPERIMENT_PACKAGE>
            """
        )
        row = parse_package("test", package)[0]
        self.assertEqual(row["run_accession"], "ERR1")
        self.assertEqual(row["bases"], 1234)
        self.assertEqual(row["tissue"], "Cerebellum")
        self.assertEqual(row["priority_rank"], 2)


if __name__ == "__main__":
    unittest.main()
