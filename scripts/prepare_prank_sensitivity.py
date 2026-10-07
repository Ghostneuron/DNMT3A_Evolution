#!/usr/bin/env python3
"""Prepare and clean a rooted-tree PRANK codon sensitivity alignment."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta  # noqa: E402


CURATED_CDS = ROOT / "results/01_curated/DNMT3A_representative_cds.fna"
TREE = ROOT / "results/03_tree/mammals/DNMT3A_mammals_tree_HyPhy_unique_rooted.nwk"
OUTDIR = ROOT / "results/02_alignment/sensitivity/prank"
INPUT = OUTDIR / "DNMT3A_mammals_unique_ungapped_no_stop.fasta"
PRANK_ALIGNMENT = OUTDIR / "DNMT3A_PRANK.best.fas"
CLEAN_ALIGNMENT = OUTDIR / "DNMT3A_PRANK_clean_mammals_HyPhy_unique.fasta"
CROSSWALK = OUTDIR / "PRANK_human_coordinate_crosswalk.tsv"
SEQUENCE_QC = OUTDIR / "PRANK_sequence_qc.tsv"
SUMMARY = OUTDIR / "PRANK_alignment_qc_summary.json"
HUMAN = "Homo_sapiens"
MIN_OCCUPANCY = 0.70


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start:start + width] + "\n")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def prepare() -> None:
    curated = {record.identifier: record.sequence.upper() for record in read_fasta(CURATED_CDS)}
    tree = Phylo.read(TREE, "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    if set(tips) - set(curated):
        raise RuntimeError(f"Missing curated CDS for {sorted(set(tips) - set(curated))[:5]}")
    records: list[tuple[str, str]] = []
    ambiguous_input_codons = 0
    for identifier in tips:
        sequence = curated[identifier]
        if len(sequence) % 3:
            raise RuntimeError(f"{identifier}: CDS length is not divisible by three")
        if CODON_TABLE.get(sequence[-3:]) != "*":
            raise RuntimeError(f"{identifier}: expected terminal stop codon")
        sequence = sequence[:-3]
        codons: list[str] = []
        for start in range(0, len(sequence), 3):
            codon = sequence[start:start + 3]
            if CODON_TABLE.get(codon) == "*":
                raise RuntimeError(f"{identifier}: internal stop codon")
            if CODON_TABLE.get(codon) is None:
                codon = "NNN"
                ambiguous_input_codons += 1
            codons.append(codon)
        sequence = "".join(codons)
        records.append((identifier, sequence))
    write_fasta(INPUT, records)
    print(
        f"Wrote {len(records)} ungapped CDS records without terminal stops to {INPUT}; "
        f"normalized {ambiguous_input_codons} ambiguous codons to NNN"
    )


def clean() -> None:
    records = read_fasta(PRANK_ALIGNMENT)
    if len(records) != 250:
        raise RuntimeError(f"Expected 250 PRANK sequences, observed {len(records)}")
    if len({record.identifier for record in records}) != len(records):
        raise RuntimeError("Duplicate PRANK identifiers")
    lengths = {len(record.sequence) for record in records}
    if len(lengths) != 1:
        raise RuntimeError("PRANK output sequences have unequal lengths")
    alignment_length = next(iter(lengths))
    if alignment_length % 3:
        raise RuntimeError("PRANK alignment length is not divisible by three")
    sequences: dict[str, str] = {}
    masked_ambiguous_codons = 0
    for record in records:
        raw = record.sequence.upper()
        normalized: list[str] = []
        for start in range(0, len(raw), 3):
            codon = raw[start:start + 3]
            if CODON_TABLE.get(codon) == "*":
                raise RuntimeError(f"{record.identifier}: PRANK output contains a stop codon")
            if codon != "---" and (
                len(codon) != 3
                or "-" in codon
                or any(base not in "ACGT" for base in codon)
                or CODON_TABLE.get(codon) is None
            ):
                codon = "---"
                masked_ambiguous_codons += 1
            normalized.append(codon)
        sequences[record.identifier] = "".join(normalized)
    expected_tips = {tip.name for tip in Phylo.read(TREE, "newick").get_terminals()}
    if set(sequences) != expected_tips:
        raise RuntimeError("PRANK alignment and guide-tree tip sets differ")

    codon_columns = alignment_length // 3
    human_position = 0
    retained_columns: list[int] = []
    crosswalk_rows: list[dict[str, object]] = []
    insertion_columns = 0
    for column in range(codon_columns):
        codons = [
            sequence[column * 3:(column + 1) * 3]
            for sequence in sequences.values()
        ]
        human_codon = sequences[HUMAN][column * 3:(column + 1) * 3]
        if human_codon != "---":
            human_position += 1
        occupancy = sum(codon != "---" for codon in codons) / len(codons)
        if occupancy >= MIN_OCCUPANCY and human_codon != "---":
            retained_columns.append(column)
            crosswalk_rows.append({
                "filtered_codon_site_1based": len(retained_columns),
                "original_prank_codon_1based": column + 1,
                "human_DNMT3A1_aa": human_position,
                "human_codon": human_codon,
                "human_residue": CODON_TABLE[human_codon],
                "occupancy": occupancy,
            })
        elif occupancy >= MIN_OCCUPANCY and human_codon == "---":
            insertion_columns += 1
    if human_position != 912:
        raise RuntimeError(f"Human sequence maps to {human_position} residues, expected 912")

    clean_records: list[tuple[str, str]] = []
    sequence_rows: list[dict[str, object]] = []
    for identifier, sequence in sequences.items():
        clean_sequence = "".join(
            sequence[column * 3:(column + 1) * 3] for column in retained_columns
        )
        clean_records.append((identifier, clean_sequence))
        non_gap = sum(
            clean_sequence[start:start + 3] != "---"
            for start in range(0, len(clean_sequence), 3)
        )
        sequence_rows.append({
            "safe_id": identifier,
            "retained_non_gap_codons": non_gap,
            "retained_fraction": non_gap / len(retained_columns),
        })
    write_fasta(CLEAN_ALIGNMENT, clean_records)
    write_tsv(
        CROSSWALK, crosswalk_rows,
        [
            "filtered_codon_site_1based", "original_prank_codon_1based",
            "human_DNMT3A1_aa", "human_codon", "human_residue", "occupancy",
        ],
    )
    write_tsv(
        SEQUENCE_QC, sequence_rows,
        ["safe_id", "retained_non_gap_codons", "retained_fraction"],
    )
    summary = {
        "sequences": len(records),
        "unaligned_human_residues": human_position,
        "prank_codon_columns": codon_columns,
        "retained_human_mapped_columns": len(retained_columns),
        "minimum_occupancy": MIN_OCCUPANCY,
        "high_occupancy_nonhuman_insertion_columns_excluded": insertion_columns,
        "ambiguous_codons_masked": masked_ambiguous_codons,
        "invalid_codons_after_masking": 0,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "clean"])
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    prepare() if args.stage == "prepare" else clean()
