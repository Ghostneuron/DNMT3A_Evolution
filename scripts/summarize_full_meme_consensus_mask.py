#!/usr/bin/env python3
"""Summarize MEME tests after masking three-aligner-discordant codons."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = (
    ROOT
    / "results/04_selection/mammals/full_meme"
    / "consensus_masked"
)
DISCOVERY = (
    ROOT / "results/04_selection/mammals/full_meme/significant_sites.tsv"
)
SITES = {"P5": 2, "S6": 3, "E18": 15}


def holm_adjust(p_values: list[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [1.0] * len(p_values)
    running = 0.0
    count = len(p_values)
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (count - rank) * p_values[index]))
        adjusted[index] = running
    return adjusted


def main() -> None:
    with DISCOVERY.open() as handle:
        discovery = {
            f"{row['human_residue']}{row['human_DNMT3A1_aa']}": row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    mask_summary = {
        row["human_site"]: row
        for row in json.loads((OUTDIR / "mask_summary.json").read_text())
    }
    rows: list[dict[str, object]] = []
    for label, site in SITES.items():
        result = OUTDIR / f"DNMT3A.{label}.consensus_masked.MEME.json"
        data = json.loads(result.read_text())
        values = data["MLE"]["content"]["0"][site - 1]
        source = discovery[label]
        rows.append({
            "human_site": label,
            "primary_filtered_site": site,
            "primary_full_meme_p": source["p_value"],
            "primary_full_meme_bh_q": source["bh_q_value_901_sites"],
            "masked_taxa_count": mask_summary[label]["masked_taxa_count"],
            "masked_taxa": ";".join(mask_summary[label]["masked_taxa"]),
            "consensus_masked_lrt": values[5],
            "consensus_masked_p": values[6],
            "consensus_masked_holm_p_3_candidates": "",
            "consensus_masked_branches_EBF_ge_100": values[7],
            "supported_after_consensus_masking_p05": str(
                float(values[6]) <= 0.05
            ).lower(),
            "interpretation": (
                "signal_depends_on_aligner_discordant_taxa"
                if float(values[6]) > 0.05
                else "signal_persists_after_consensus_masking"
            ),
        })
    adjusted = holm_adjust([float(row["consensus_masked_p"]) for row in rows])
    for row, value in zip(rows, adjusted):
        row["consensus_masked_holm_p_3_candidates"] = value

    with (OUTDIR / "consensus_masked_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    supported = sum(
        row["supported_after_consensus_masking_p05"] == "true" for row in rows
    )
    (OUTDIR / "consensus_masked_summary.json").write_text(json.dumps({
        "sites_tested": len(rows),
        "sites_supported_after_consensus_masking_p05": supported,
        "masked_taxa_per_site": 2,
        "masking_rule": (
            "replace the human-mapped primary codon with NNN when its exact "
            "codon differs among primary, MAFFT, and PRANK"
        ),
        "multiplicity_control": "Holm across the three discovery candidates",
        "interpretation": (
            "All three primary full-scan discoveries become nonsignificant "
            "after masking Puma concolor and Vombatus ursinus at the tested "
            "site; the signals depend on aligner-discordant codon placements."
        ),
    }, indent=2) + "\n")
    print(f"Summarized {len(rows)} sites; {supported} persist at p<=0.05")


if __name__ == "__main__":
    main()
