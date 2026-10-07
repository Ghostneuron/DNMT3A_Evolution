#!/usr/bin/env python3
"""Prepare the pre-MACSE MAFFT codon projection for FUBAR sensitivity analysis."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import (  # noqa: E402
    CODON_TABLE,
    domain_for,
    load_domains,
    read_fasta,
    safe_id,
    write_fasta,
    write_tsv,
)


SOURCE = ROOT / "results/02_alignment/DNMT3A_MAFFT_projected_codons.fasta"
TAXONOMY = ROOT / "results/01_curated/taxon_metadata.tsv"
DEDUPLICATION = ROOT / "results/04_selection/mammals/deduplication_audit.tsv"
OUT = ROOT / "results/02_alignment/sensitivity"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    records = read_fasta(SOURCE)
    lengths = {len(record.sequence) for record in records}
    if len(lengths) != 1 or next(iter(lengths)) % 3:
        raise SystemExit("ERROR: invalid MAFFT codon projection dimensions")
    codon_columns = next(iter(lengths)) // 3
    codons: dict[str, list[str]] = {}
    for record in records:
        row: list[str] = []
        for start in range(0, len(record.sequence), 3):
            codon = record.sequence[start:start + 3]
            if codon == "---":
                row.append(codon)
            elif "-" in codon or any(base not in "ACGT" for base in codon):
                row.append("---")
            elif CODON_TABLE.get(codon) == "*":
                row.append("---")
            else:
                row.append(codon)
        codons[record.identifier] = row

    threshold = 0.70
    retained = [
        column for column in range(codon_columns)
        if sum(row[column] != "---" for row in codons.values()) / len(codons) >= threshold
    ]
    filtered = {
        identifier: "".join(row[column] for column in retained)
        for identifier, row in codons.items()
    }
    OUT.mkdir(parents=True, exist_ok=True)
    write_fasta(filtered.items(), OUT / "DNMT3A_MAFFT_clean_all.fasta")

    human_id = safe_id("Homo sapiens")
    human_position = 0
    filtered_site = 0
    domains = load_domains(ROOT)
    retained_set = set(retained)
    crosswalk: list[dict[str, object]] = []
    for column in range(codon_columns):
        human_codon = codons[human_id][column]
        position: int | None = None
        if human_codon != "---":
            human_position += 1
            position = human_position
        if column in retained_set:
            filtered_site += 1
            occupancy = sum(row[column] != "---" for row in codons.values()) / len(codons)
            crosswalk.append({
                "filtered_codon_site_1based": filtered_site,
                "original_mafft_codon_1based": column + 1,
                "human_DNMT3A1_aa": position if position is not None else "",
                "human_residue": CODON_TABLE.get(human_codon, "-") if human_codon != "---" else "-",
                "domain": domain_for(position, domains),
                "occupancy": f"{occupancy:.6f}",
            })
    write_tsv(OUT / "MAFFT_human_coordinate_crosswalk.tsv", crosswalk)

    mammals = {
        row["safe_id"] for row in read_tsv(TAXONOMY) if row["is_mammal"] == "true"
    }
    removed = {
        row["removed_safe_id"] for row in read_tsv(DEDUPLICATION)
    }
    keep = mammals - removed
    mammal_rows = [
        (identifier, filtered[identifier]) for identifier in filtered if identifier in keep
    ]
    if len(mammal_rows) != 250:
        raise SystemExit(f"ERROR: expected 250 unique mammalian sequences, observed {len(mammal_rows)}")
    write_fasta(mammal_rows, OUT / "DNMT3A_MAFFT_clean_mammals_HyPhy_unique.fasta")

    summary = {
        "sequences_all": len(records),
        "input_codon_columns": codon_columns,
        "retained_codon_columns": len(retained),
        "human_mapped_residues": human_position,
        "mammal_unique_sequences": len(mammal_rows),
    }
    (OUT / "MAFFT_sensitivity_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
