#!/usr/bin/env python3
"""Prepare, project, and clean a COBALT protein-alignment sensitivity set."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import (  # noqa: E402
    CODON_TABLE,
    domain_for,
    load_domains,
    read_fasta,
    thread_cds_to_protein_alignment,
    write_fasta,
    write_tsv,
)


CURATED_PROTEIN = ROOT / "results/01_curated/DNMT3A_representative_proteins.faa"
CURATED_CDS = ROOT / "results/01_curated/DNMT3A_representative_cds.fna"
TREE = ROOT / "results/03_tree/mammals/DNMT3A_mammals_tree_HyPhy_unique_rooted.nwk"
OUT = ROOT / "results/02_alignment/sensitivity/cobalt"
INPUT = OUT / "DNMT3A_mammals_unique_proteins.faa"
COBALT_RAW_ALIGNMENT = OUT / "DNMT3A_COBALT_AA.fasta"
COBALT_ALIGNMENT = OUT / "DNMT3A_COBALT_AA_named.fasta"
PROJECTED = OUT / "DNMT3A_COBALT_projected_codons.fasta"
CLEAN = OUT / "DNMT3A_COBALT_clean_mammals_HyPhy_unique.fasta"
CROSSWALK = OUT / "COBALT_human_coordinate_crosswalk.tsv"
SEQUENCE_QC = OUT / "COBALT_sequence_qc.tsv"
SUMMARY = OUT / "COBALT_alignment_qc_summary.json"
HUMAN = "Homo_sapiens"
MIN_OCCUPANCY = 0.70


def prepare() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    proteins = {
        record.identifier: record.sequence.upper()
        for record in read_fasta(CURATED_PROTEIN)
    }
    tips = [tip.name for tip in Phylo.read(TREE, "newick").get_terminals()]
    missing = set(tips) - set(proteins)
    if missing:
        raise RuntimeError(f"Missing curated proteins: {sorted(missing)[:5]}")
    records = [(identifier, proteins[identifier]) for identifier in tips]
    write_fasta(records, INPUT)
    print(f"Wrote {len(records)} proteins to {INPUT}")


def clean() -> None:
    aligned = read_fasta(COBALT_RAW_ALIGNMENT)
    if len(aligned) != 250:
        raise RuntimeError(f"Expected 250 COBALT sequences, observed {len(aligned)}")
    if len({record.identifier for record in aligned}) != len(aligned):
        raise RuntimeError("Duplicate COBALT identifiers")
    lengths = {len(record.sequence) for record in aligned}
    if len(lengths) != 1:
        raise RuntimeError("COBALT output is not rectangular")

    input_records = read_fasta(INPUT)
    source = {record.identifier: record.sequence.upper() for record in input_records}
    normalized: list[tuple[str, str]] = []
    for index, record in enumerate(aligned, start=1):
        # COBALT 3.0 emits local numeric identifiers unless parseable NCBI
        # accessions are supplied. Its local IDs preserve FASTA input order.
        expected_local_id = f"lcl|{index}"
        if record.identifier != expected_local_id:
            raise RuntimeError(
                f"Expected COBALT identifier {expected_local_id}, observed "
                f"{record.identifier}"
            )
        identifier = input_records[index - 1].identifier
        if record.sequence.replace("-", "").upper() != source[identifier]:
            raise RuntimeError(f"COBALT changed residue content for {identifier}")
        normalized.append((identifier, record.sequence.upper()))
    write_fasta(normalized, COBALT_ALIGNMENT)

    thread_cds_to_protein_alignment(COBALT_ALIGNMENT, CURATED_CDS, PROJECTED)
    projected = read_fasta(PROJECTED)
    sequences: dict[str, str] = {}
    ambiguous_masked = 0
    for record in projected:
        normalized: list[str] = []
        for start in range(0, len(record.sequence), 3):
            value = record.sequence[start:start + 3].upper()
            if value != "---" and (
                len(value) != 3
                or "-" in value
                or any(base not in "ACGT" for base in value)
                or CODON_TABLE.get(value) in {None, "*"}
            ):
                value = "---"
                ambiguous_masked += 1
            normalized.append(value)
        sequences[record.identifier] = "".join(normalized)

    if set(sequences) != set(source):
        raise RuntimeError("COBALT projection tip set differs from input")
    codon_columns = next(iter(lengths))  # amino-acid columns
    human_position = 0
    retained: list[int] = []
    crosswalk_rows: list[dict[str, object]] = []
    domains = load_domains(ROOT)
    insertion_columns = 0
    for column in range(codon_columns):
        values = [
            sequence[column * 3:(column + 1) * 3]
            for sequence in sequences.values()
        ]
        human_codon = sequences[HUMAN][column * 3:(column + 1) * 3]
        position: int | None = None
        if human_codon != "---":
            human_position += 1
            position = human_position
        occupancy = sum(value != "---" for value in values) / len(values)
        if occupancy >= MIN_OCCUPANCY and human_codon != "---":
            retained.append(column)
            crosswalk_rows.append({
                "filtered_codon_site_1based": len(retained),
                "original_cobalt_codon_1based": column + 1,
                "human_DNMT3A1_aa": position,
                "human_codon": human_codon,
                "human_residue": CODON_TABLE[human_codon],
                "domain": domain_for(position, domains),
                "occupancy": occupancy,
            })
        elif occupancy >= MIN_OCCUPANCY:
            insertion_columns += 1
    if human_position != 912:
        raise RuntimeError(f"Human maps to {human_position} residues, expected 912")

    clean_records: list[tuple[str, str]] = []
    sequence_rows: list[dict[str, object]] = []
    for identifier, sequence in sequences.items():
        clean_sequence = "".join(
            sequence[column * 3:(column + 1) * 3] for column in retained
        )
        clean_records.append((identifier, clean_sequence))
        non_gap = sum(
            clean_sequence[start:start + 3] != "---"
            for start in range(0, len(clean_sequence), 3)
        )
        sequence_rows.append({
            "safe_id": identifier,
            "retained_non_gap_codons": non_gap,
            "retained_fraction": non_gap / len(retained),
        })
    write_fasta(clean_records, CLEAN)
    write_tsv(CROSSWALK, crosswalk_rows)
    write_tsv(SEQUENCE_QC, sequence_rows)
    summary = {
        "cobalt_version": "3.0.0",
        "constraint_mode": "norps_local_similarity_only",
        "sequences": len(aligned),
        "aligned_amino_acid_columns": codon_columns,
        "human_residues": human_position,
        "retained_human_mapped_columns": len(retained),
        "minimum_occupancy": MIN_OCCUPANCY,
        "high_occupancy_nonhuman_insertion_columns_excluded": insertion_columns,
        "ambiguous_codons_masked": ambiguous_masked,
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
