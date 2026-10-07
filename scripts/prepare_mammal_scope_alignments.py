#!/usr/bin/env python3
"""Refilter raw MACSE/MAFFT codon alignments within the 250-mammal scope."""

from __future__ import annotations

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
    write_fasta,
    write_tsv,
)


TREE = ROOT / "results/03_tree/mammals/DNMT3A_mammals_tree_HyPhy_unique_rooted.nwk"
OUT = ROOT / "results/02_alignment/sensitivity/mammal_scope"
SOURCES = {
    "MACSE": ROOT / "results/02_alignment/DNMT3A_MACSE_NT.fasta",
    "MAFFT": ROOT / "results/02_alignment/DNMT3A_MAFFT_projected_codons.fasta",
}
HUMAN = "Homo_sapiens"
MIN_OCCUPANCY = 0.70


def valid(codon: str) -> str:
    codon = codon.upper().replace("U", "T")
    if (
        codon == "---"
        or len(codon) != 3
        or "!" in codon
        or "-" in codon
        or any(base not in "ACGT" for base in codon)
        or CODON_TABLE.get(codon) in {None, "*"}
    ):
        return "---"
    return codon


def prepare(label: str, source: Path, tips: list[str]) -> dict[str, object]:
    records = read_fasta(source)
    available = {record.identifier: record.sequence for record in records}
    missing = set(tips) - set(available)
    if missing:
        raise RuntimeError(f"{label} lacks mammal tips: {sorted(missing)[:5]}")
    lengths = {len(available[tip]) for tip in tips}
    if len(lengths) != 1 or next(iter(lengths)) % 3:
        raise RuntimeError(f"{label} source is not a rectangular codon alignment")
    codon_columns = next(iter(lengths)) // 3
    codons = {
        tip: [
            valid(available[tip][start:start + 3])
            for start in range(0, codon_columns * 3, 3)
        ]
        for tip in tips
    }

    human_position = 0
    retained: list[int] = []
    crosswalk: list[dict[str, object]] = []
    domains = load_domains(ROOT)
    insertion_columns = 0
    for column in range(codon_columns):
        human_codon = codons[HUMAN][column]
        position: int | None = None
        if human_codon != "---":
            human_position += 1
            position = human_position
        occupancy = sum(row[column] != "---" for row in codons.values()) / len(tips)
        if occupancy >= MIN_OCCUPANCY and human_codon != "---":
            retained.append(column)
            crosswalk.append({
                "filtered_codon_site_1based": len(retained),
                f"original_{label.lower()}_codon_1based": column + 1,
                "human_DNMT3A1_aa": position,
                "human_codon": human_codon,
                "human_residue": CODON_TABLE[human_codon],
                "domain": domain_for(position, domains),
                "mammal_occupancy": occupancy,
            })
        elif occupancy >= MIN_OCCUPANCY:
            insertion_columns += 1
    if human_position != 912:
        raise RuntimeError(f"{label} maps human to {human_position}, expected 912")

    directory = OUT / label.lower()
    alignment = directory / f"DNMT3A_{label}_mammal_scope_HyPhy_unique.fasta"
    write_fasta(
        [
            (
            tip,
            "".join(codons[tip][column] for column in retained),
            )
            for tip in tips
        ],
        alignment,
    )
    write_tsv(directory / f"{label}_mammal_scope_crosswalk.tsv", crosswalk)
    summary = {
        "method": label,
        "sequences": len(tips),
        "source_codon_columns": codon_columns,
        "retained_human_mapped_columns": len(retained),
        "minimum_mammal_occupancy": MIN_OCCUPANCY,
        "human_residues": human_position,
        "high_occupancy_nonhuman_insertion_columns_excluded": insertion_columns,
        "s97_mammal_occupancy": next(
            row["mammal_occupancy"]
            for row in crosswalk
            if row["human_DNMT3A1_aa"] == 97
        ),
    }
    (directory / f"{label}_mammal_scope_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    return summary


def main() -> None:
    tips = [tip.name for tip in Phylo.read(TREE, "newick").get_terminals()]
    summaries = [prepare(label, source, tips) for label, source in SOURCES.items()]
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
