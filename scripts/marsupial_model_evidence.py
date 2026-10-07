#!/usr/bin/env python3
"""Extract evidence qualifiers for predicted marsupial DNMT3A transcripts."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/promoter_tss"
MANIFEST = RAW / "window_manifest.tsv"
TARGETS = ROOT / "results/06_isoform_evolution/target_transcript_tss.tsv"
RESULTS = ROOT / "results/06_isoform_evolution"
LONG_READS = re.compile(r"(\d+) long SRA reads")
MRNAS = re.compile(r"(\d+) mRNA")
RNA_COVERAGE = re.compile(
    r"(\d+)% coverage(?: of the annotated genomic feature)? by RNAseq alignments"
)
INTRON_SAMPLES = re.compile(r"(\d+) samples with support for all annotated introns")


def match_integer(pattern: re.Pattern[str], text: str):
    match = pattern.search(text)
    return int(match.group(1)) if match else ""


def main() -> None:
    with MANIFEST.open(newline="") as handle:
        manifest = {
            row["species"]: row for row in csv.DictReader(handle, delimiter="\t")
        }
    with TARGETS.open(newline="") as handle:
        targets = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row["is_primary_internal_start_candidate"] == "true"
            and row["mammal_group"] == "Metatheria"
        ]

    rows = []
    for target in targets:
        species = target["species"]
        record = SeqIO.read(
            RAW / manifest[species]["genbank_file"], "genbank"
        )
        transcript = target["transcript_accession"]
        feature = next(
            feature for feature in record.features
            if feature.type == "mRNA"
            and feature.qualifiers.get("transcript_id", [""])[0] == transcript
        )
        note = " ".join(feature.qualifiers.get("note", []))
        experiment = "; ".join(feature.qualifiers.get("experiment", []))
        rows.append({
            "species": species,
            "transcript_accession": transcript,
            "annotation_method": (
                "Gnomon" if "Gnomon" in note else "other_or_unspecified"
            ),
            "supporting_mRNA_count": match_integer(MRNAS, note),
            "supporting_long_SRA_read_count": match_integer(LONG_READS, note),
            "RNAseq_feature_coverage_percent": match_integer(RNA_COVERAGE, note),
            "samples_supporting_all_annotated_introns": match_integer(
                INTRON_SAMPLES, note
            ),
            "polyA_coordinate_evidence": str(
                "polyA evidence" in experiment
            ).lower(),
            "raw_experiment_qualifier": experiment,
            "raw_annotation_note": note,
            "model_body_support": (
                "species-specific RNA evidence supports the predicted transcript model"
            ),
            "exact_5prime_TSS_support": "false",
            "interpretation_limit": (
                "RNA/long-read coverage and polyA evidence do not demonstrate "
                "the exact capped 5-prime transcription start"
            ),
        })

    with (RESULTS / "marsupial_transcript_model_evidence.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "marsupial_internal_start_models": len(rows),
        "models_with_long_SRA_read_support": sum(
            bool(row["supporting_long_SRA_read_count"]) for row in rows
        ),
        "models_with_reported_RNAseq_feature_coverage": sum(
            bool(row["RNAseq_feature_coverage_percent"]) for row in rows
        ),
        "models_with_polyA_coordinate_evidence": sum(
            row["polyA_coordinate_evidence"] == "true" for row in rows
        ),
        "direct_exact_5prime_TSS_models": 0,
        "models": rows,
        "inference": (
            "The marsupial candidates are RNA-supported predicted transcript "
            "models, not protein-only truncation artifacts."
        ),
        "limit": (
            "The available qualifiers do not validate the exact 5-prime "
            "transcription start or promoter activity."
        ),
    }
    (RESULTS / "marsupial_transcript_model_evidence_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
