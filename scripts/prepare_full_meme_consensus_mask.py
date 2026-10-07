#!/usr/bin/env python3
"""Mask discovery-site codons that disagree among all three alignments."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALIGNMENTS = {
    "primary": (
        ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
    ),
    "mafft": (
        ROOT
        / "results/02_alignment/sensitivity"
        / "DNMT3A_MAFFT_clean_mammals_HyPhy_unique.fasta"
    ),
    "prank": (
        ROOT
        / "results/02_alignment/sensitivity/prank"
        / "DNMT3A_PRANK_clean_mammals_HyPhy_unique.fasta"
    ),
}
CROSSWALKS = {
    "primary": ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv",
    "mafft": (
        ROOT
        / "results/02_alignment/sensitivity"
        / "MAFFT_human_coordinate_crosswalk.tsv"
    ),
    "prank": (
        ROOT
        / "results/02_alignment/sensitivity/prank"
        / "PRANK_human_coordinate_crosswalk.tsv"
    ),
}
OUTDIR = (
    ROOT
    / "results/04_selection/mammals/full_meme"
    / "consensus_masked"
)
SITES = {5: "P5", 6: "S6", 18: "E18"}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            identifier = line[1:].split()[0]
            records[identifier] = ""
        elif identifier is None:
            raise RuntimeError(f"Sequence before header in {path}")
        else:
            records[identifier] += line
    return records


def write_fasta(path: Path, records: dict[str, str], width: int = 80) -> None:
    with path.open("w") as handle:
        for identifier, sequence in records.items():
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")


def crosswalk(path: Path) -> dict[int, int]:
    with path.open() as handle:
        return {
            int(row["human_DNMT3A1_aa"]): int(
                row["filtered_codon_site_1based"]
            )
            for row in csv.DictReader(handle, delimiter="\t")
            if row["human_DNMT3A1_aa"]
        }


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    alignments = {
        name: read_fasta(path) for name, path in ALIGNMENTS.items()
    }
    maps = {name: crosswalk(path) for name, path in CROSSWALKS.items()}
    tip_sets = {frozenset(records) for records in alignments.values()}
    if len(tip_sets) != 1:
        raise RuntimeError("Alignment tip sets differ")

    audit: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    for human_position, label in SITES.items():
        primary_site = maps["primary"][human_position]
        ambiguous: list[str] = []
        for identifier in alignments["primary"]:
            codons = {
                name: records[identifier][
                    (maps[name][human_position] - 1) * 3 :
                    maps[name][human_position] * 3
                ]
                for name, records in alignments.items()
            }
            if len(set(codons.values())) > 1:
                ambiguous.append(identifier)
                audit.append({
                    "human_site": label,
                    "primary_filtered_site": primary_site,
                    "safe_id": identifier,
                    "primary_codon": codons["primary"],
                    "mafft_codon": codons["mafft"],
                    "prank_codon": codons["prank"],
                    "mask": "NNN",
                    "reason": "human_mapped_codon_disagrees_among_alignments",
                })
        if not ambiguous:
            raise RuntimeError(f"No ambiguous taxa found at {label}")

        masked = dict(alignments["primary"])
        start = (primary_site - 1) * 3
        for identifier in ambiguous:
            sequence = masked[identifier]
            masked[identifier] = sequence[:start] + "NNN" + sequence[start + 3 :]
        output = OUTDIR / f"DNMT3A.{label}.consensus_masked.fasta"
        write_fasta(output, masked)
        for identifier in ambiguous:
            if masked[identifier][start : start + 3] != "NNN":
                raise RuntimeError(f"Failed to mask {identifier} at {label}")
        summaries.append({
            "human_site": label,
            "primary_filtered_site": primary_site,
            "masked_taxa_count": len(ambiguous),
            "masked_taxa": sorted(ambiguous),
            "unmasked_taxa_count": len(masked) - len(ambiguous),
            "output": str(output.relative_to(ROOT)),
        })

    with (OUTDIR / "mask_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(audit[0]))
        writer.writeheader()
        writer.writerows(audit)
    (OUTDIR / "mask_summary.json").write_text(
        json.dumps(summaries, indent=2) + "\n"
    )
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
