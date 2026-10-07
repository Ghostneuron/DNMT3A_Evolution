#!/usr/bin/env python3
"""Test the tammar DNMT3A internal first exon in public ONT cDNA reads."""

from __future__ import annotations

import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import ahocorasick
from Bio import Align, SeqIO
from Bio.Seq import Seq

from dunnart_transcriptome_first_exon import junction_seeds, spaced_seeds


ROOT = Path(__file__).resolve().parents[1]
FASTQ = (
    ROOT / "data/raw/marsupial_transcriptomes/SRR30198375_1.fastq.gz"
)
GENBANK = (
    ROOT / "data/raw/promoter_tss/"
    "Notamacropus_eugenii_NC_092872.1_DNMT3A.gb"
)
RESULTS = ROOT / "results/06_isoform_evolution"
TARGET = "XM_072633891.1"


def target_parts() -> list[str]:
    record = SeqIO.read(GENBANK, "genbank")
    feature = next(
        item for item in record.features
        if item.type == "mRNA"
        and item.qualifiers.get("transcript_id", [""])[0] == TARGET
    )
    return [
        str(part.extract(record.seq)).upper()
        for part in feature.location.parts
    ]


def build_seeds() -> list[dict[str, str]]:
    parts = target_parts()
    rows = []
    for component, component_seeds in (
        ("first_exon", spaced_seeds(parts[0], count=20)),
        ("first_junction", junction_seeds(parts[0], parts[1])),
        ("shared_dnmt3a_exon", spaced_seeds(parts[1], count=10)),
    ):
        for _, seed in component_seeds:
            for orientation, pattern in (
                ("forward", seed),
                ("reverse_complement", str(Seq(seed).reverse_complement())),
            ):
                rows.append({
                    "component": component,
                    "orientation": orientation,
                    "seed": pattern,
                })
    return rows


def classify_read(
    sequence: str,
    seeds: list[dict[str, str]],
) -> dict[tuple[str, str], int]:
    sequence = sequence.upper()
    matches: dict[tuple[str, str], int] = defaultdict(int)
    for row in seeds:
        if row["seed"] in sequence:
            matches[(row["component"], row["orientation"])] += 1
    return dict(matches)


def main() -> None:
    if not FASTQ.exists():
        raise SystemExit(f"Missing FASTQ: {FASTQ}")
    seeds = build_seeds()
    seed_lookup: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in seeds:
        seed_lookup[row["seed"]].append(row)
    automaton = ahocorasick.Automaton()
    for seed in seed_lookup:
        automaton.add_word(seed, seed)
    automaton.make_automaton()
    evidence_rows = []
    matched_reads: list[tuple[str, str, str]] = []
    total_reads = 0
    with gzip.open(FASTQ, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip().upper()
            plus = handle.readline()
            quality = handle.readline().strip()
            if not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ near read {total_reads + 1}")
            total_reads += 1
            matched_seed_strings = {seed for _, seed in automaton.iter(sequence)}
            if not matched_seed_strings:
                continue
            read_id = header[1:].split()[0]
            matches: dict[tuple[str, str], int] = defaultdict(int)
            for seed in matched_seed_strings:
                for row in seed_lookup[seed]:
                    matches[(row["component"], row["orientation"])] += 1
            matched_reads.append((header.rstrip(), sequence, quality))
            for (component, orientation), count in matches.items():
                evidence_rows.append({
                    "run_accession": "SRR30198375",
                    "read_id": read_id,
                    "read_length": len(sequence),
                    "component": component,
                    "orientation": orientation,
                    "distinct_exact_35mer_hits": count,
                })

    evidence_rows.sort(key=lambda row: (
        row["read_id"], row["component"], row["orientation"]
    ))
    evidence_output = RESULTS / "tammar_ont_first_exon_evidence.tsv"
    with evidence_output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t",
            fieldnames=[
                "run_accession", "read_id", "read_length", "component",
                "orientation", "distinct_exact_35mer_hits",
            ],
        )
        writer.writeheader()
        writer.writerows(evidence_rows)
    with (RESULTS / "tammar_ont_dnmt3a_seed_reads.fastq").open("w") as handle:
        for header, sequence, quality in matched_reads:
            handle.write(f"{header}\n{sequence}\n+\n{quality}\n")

    components_by_read: dict[str, set[str]] = defaultdict(set)
    for row in evidence_rows:
        components_by_read[str(row["read_id"])].add(str(row["component"]))
    first_junction_reads = sorted(
        read_id for read_id, components in components_by_read.items()
        if "first_junction" in components
    )
    first_and_shared_reads = sorted(
        read_id for read_id, components in components_by_read.items()
        if {"first_exon", "shared_dnmt3a_exon"} <= components
    )

    read_sequences = {
        header[1:].split()[0]: sequence
        for header, sequence, _ in matched_reads
    }
    target_sequence = "".join(target_parts())
    aligner = Align.PairwiseAligner(
        mode="local",
        match_score=2,
        mismatch_score=-1,
        open_gap_score=-5,
        extend_gap_score=-1,
    )
    alignment_rows = []
    for read_id in first_junction_reads:
        read_sequence = read_sequences[read_id]
        alternatives = []
        for orientation, oriented_sequence in (
            ("forward", read_sequence),
            ("reverse_complement", str(Seq(read_sequence).reverse_complement())),
        ):
            alignment = aligner.align(target_sequence, oriented_sequence)[0]
            alternatives.append((alignment.score, orientation, alignment))
        score, orientation, alignment = max(alternatives, key=lambda item: item[0])
        counts = alignment.counts()
        coordinates = alignment.coordinates
        alignment_rows.append({
            "read_id": read_id,
            "read_length": len(read_sequence),
            "orientation_to_model": orientation,
            "alignment_score": f"{score:.1f}",
            "model_start_0based": int(coordinates[0, 0]),
            "model_end_0based_exclusive": int(coordinates[0, -1]),
            "read_start_0based": int(coordinates[1, 0]),
            "read_end_0based_exclusive": int(coordinates[1, -1]),
            "aligned_identical_bases": counts.identities,
            "aligned_nongap_bases": counts.aligned,
            "aligned_identity": f"{counts.identities / counts.aligned:.4f}",
            "reaches_annotated_transcript_base_1": str(
                int(coordinates[0, 0]) == 0
            ).lower(),
        })
    with (RESULTS / "tammar_ont_junction_read_alignments.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(alignment_rows[0])
        )
        writer.writeheader()
        writer.writerows(alignment_rows)

    summary = {
        "target_transcript": TARGET,
        "run_accession": "SRR30198375",
        "sample": "day-148 pouch-young testis",
        "platform": "Oxford Nanopore PromethION PCR-cDNA",
        "fastq": str(FASTQ.relative_to(ROOT)),
        "fastq_expected_md5": "ea114050c864de79d3634fbc2a46627a",
        "reads_scanned": total_reads,
        "reads_with_any_exact_35mer_seed": len(matched_reads),
        "reads_with_exact_first_junction_seed": first_junction_reads,
        "reads_with_first_and_shared_exon_seeds": first_and_shared_reads,
        "junction_reads_reaching_annotated_transcript_base_1": [
            row["read_id"] for row in alignment_rows
            if row["reaches_annotated_transcript_base_1"] == "true"
        ],
        "junction_read_alignment_table": (
            "results/06_isoform_evolution/"
            "tammar_ont_junction_read_alignments.tsv"
        ),
        "interpretation_limit": (
            "A junction-spanning cDNA read supports the internal transcript "
            "structure, but PCR-cDNA read starts do not establish an exact "
            "capped TSS."
        ),
    }
    (RESULTS / "tammar_ont_first_exon_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
