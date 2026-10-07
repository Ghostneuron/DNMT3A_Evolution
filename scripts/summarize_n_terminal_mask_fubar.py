#!/usr/bin/env python3
"""Summarize FUBAR results after the taxon-wide N-terminal reliability mask."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MASK_DIR = ROOT / "results/02_alignment/reliability_mask"
SELECTION = ROOT / "results/04_selection/mammals"
FUBAR = {
    "primary": SELECTION / "DNMT3A.FUBAR.json",
    "masked": SELECTION / "DNMT3A.n_terminal_reliability_mask.FUBAR.json",
    "mafft": SELECTION / "DNMT3A.mafft_sensitivity.FUBAR.json",
    "prank": SELECTION / "DNMT3A.prank_sensitivity.FUBAR.json",
}
CROSSWALKS = {
    "primary": ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv",
    "mafft": (
        ROOT / "results/02_alignment/sensitivity/MAFFT_human_coordinate_crosswalk.tsv"
    ),
    "prank": (
        ROOT
        / "results/02_alignment/sensitivity/prank/PRANK_human_coordinate_crosswalk.tsv"
    ),
}
AUDITED = {5, 6, 12, 16, 18, 34, 114}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def crosswalk(path: Path) -> dict[int, dict[str, str]]:
    return {
        int(row["human_DNMT3A1_aa"]): row
        for row in read_tsv(path)
        if row["human_DNMT3A1_aa"]
    }


def posterior(path: Path) -> dict[int, float]:
    data = json.loads(path.read_text())
    headers = [
        item[0] if isinstance(item, list) else str(item)
        for item in data["MLE"]["headers"]
    ]
    index = headers.index("Prob[alpha<beta]")
    return {
        site: float(values[index])
        for site, values in enumerate(data["MLE"]["content"]["0"], start=1)
    }


def status(primary: float, masked: float) -> str:
    if primary >= 0.9 and masked >= 0.9:
        return "retained_above_0.90_after_mask"
    if primary < 0.9 <= masked:
        return "mask_sensitivity_emergent"
    if primary >= 0.9 > masked:
        return "lost_after_mask"
    return "below_0.90_in_primary_and_masked"


def main() -> None:
    maps = {name: crosswalk(path) for name, path in CROSSWALKS.items()}
    fits = {name: posterior(path) for name, path in FUBAR.items()}
    reliability = {
        int(row["human_DNMT3A1_aa"]): row
        for row in read_tsv(MASK_DIR / "site_reliability.tsv")
    }

    rows: list[dict[str, object]] = []
    for human_position in sorted(AUDITED):
        primary_row = maps["primary"][human_position]
        primary_site = int(primary_row["filtered_codon_site_1based"])
        mafft_site = int(maps["mafft"][human_position]["filtered_codon_site_1based"])
        prank_site = int(maps["prank"][human_position]["filtered_codon_site_1based"])
        primary_p = fits["primary"][primary_site]
        masked_p = fits["masked"][primary_site]
        reliability_row = reliability[human_position]
        rows.append({
            "human_site": f"{primary_row['human_residue']}{human_position}",
            "human_DNMT3A1_aa": human_position,
            "primary_filtered_site": primary_site,
            "three_alignment_agreement_fraction": reliability_row[
                "all_three_equal_fraction"
            ],
            "discordant_taxa_count": reliability_row["discordant_count"],
            "newly_masked_taxa_count": reliability_row[
                "newly_masked_informative_count"
            ],
            "primary_fubar_posterior": primary_p,
            "masked_fubar_posterior": masked_p,
            "posterior_change_masked_minus_primary": masked_p - primary_p,
            "mafft_fubar_posterior": fits["mafft"][mafft_site],
            "prank_fubar_posterior": fits["prank"][prank_site],
            "mask_status": status(primary_p, masked_p),
        })

    with (MASK_DIR / "candidate_fubar_comparison.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    masked_positive = []
    primary_map_by_site = {
        int(row["filtered_codon_site_1based"]): row
        for row in maps["primary"].values()
    }
    for site, value in fits["masked"].items():
        if value < 0.9:
            continue
        mapped = primary_map_by_site[site]
        masked_positive.append({
            "primary_filtered_site": site,
            "human_site": f"{mapped['human_residue']}{mapped['human_DNMT3A1_aa']}",
            "human_DNMT3A1_aa": int(mapped["human_DNMT3A1_aa"]),
            "masked_fubar_posterior": value,
            "primary_fubar_posterior": fits["primary"][site],
            "mask_status": status(fits["primary"][site], value),
        })

    summary = {
        "masked_fubar_sites_posterior_ge_0.90": masked_positive,
        "retained_primary_candidates": [
            row["human_site"]
            for row in masked_positive
            if row["mask_status"] == "retained_above_0.90_after_mask"
        ],
        "mask_sensitivity_emergent": [
            row["human_site"]
            for row in masked_positive
            if row["mask_status"] == "mask_sensitivity_emergent"
        ],
        "candidate_comparison": rows,
        "interpretation": (
            "T12 and G34 retain FUBAR posterior >=0.90 after the taxon-wide "
            "alignment-reliability mask. A16 crosses 0.90 only after masking and "
            "is therefore sensitivity-emergent, not an independent discovery."
        ),
    }
    (MASK_DIR / "fubar_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
