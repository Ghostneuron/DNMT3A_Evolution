#!/usr/bin/env python3
"""Summarize an alignment-wide DNMT3A MEME scan with site-wide correction."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/04_selection/mammals/full_meme"
RESULT = OUTDIR / "DNMT3A.full.MEME.json"
CROSSWALK = ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv"
ALL_SITES = OUTDIR / "all_sites.tsv"
SIGNIFICANT = OUTDIR / "significant_sites.tsv"
SUMMARY = OUTDIR / "summary.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bh_adjust(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values, monotone in ranked order."""
    count = len(p_values)
    order = sorted(range(count), key=p_values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for reverse_rank, index in enumerate(reversed(order), start=1):
        rank = count - reverse_rank + 1
        running = min(running, p_values[index] * count / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def find_column(headers: list[str], token: str) -> int:
    for index, header in enumerate(headers):
        if token.casefold() in header.casefold():
            return index
    raise RuntimeError(f"Missing MEME column containing {token!r}")


def main() -> None:
    data = json.loads(RESULT.read_text())
    site_count = int(data["input"]["number of sites"])
    rows = data["MLE"]["content"]["0"]
    if len(rows) != site_count:
        raise RuntimeError(f"MEME row count {len(rows)} != input sites {site_count}")
    if site_count != 901:
        raise RuntimeError(f"Expected the production 901-codon alignment, found {site_count}")

    headers = [
        item[0] if isinstance(item, list) else str(item)
        for item in data["MLE"]["headers"]
    ]
    p_index = find_column(headers, "p-value")
    lrt_index = find_column(headers, "LRT")
    alpha_index = find_column(headers, "&alpha;")
    beta_index = find_column(headers, "&beta;<sup>+</sup>")
    weight_index = find_column(headers, "p<sup>+</sup>")
    branch_index = find_column(headers, "# branches")
    meme_logl_index = find_column(headers, "MEME LogL")

    mapping = {
        int(row["filtered_codon_site_1based"]): row
        for row in read_tsv(CROSSWALK)
    }
    p_values = [float(row[p_index]) for row in rows]
    bh_values = bh_adjust(p_values)
    output: list[dict[str, object]] = []
    for site, (values, p_value, q_value) in enumerate(
        zip(rows, p_values, bh_values), start=1
    ):
        mapped = mapping[site]
        output.append({
            "filtered_codon_site_1based": site,
            "human_DNMT3A1_aa": mapped["human_DNMT3A1_aa"],
            "human_residue": mapped["human_residue"],
            "domain": mapped["domain"],
            "occupancy": mapped["occupancy"],
            "alpha": values[alpha_index],
            "beta_positive": values[beta_index],
            "positive_class_weight": values[weight_index],
            "LRT": values[lrt_index],
            "p_value": p_value,
            "bh_q_value_901_sites": q_value,
            "bonferroni_p_901_sites": min(1.0, p_value * site_count),
            "branches_EBF_ge_100": values[branch_index],
            "meme_log_likelihood": values[meme_logl_index],
            "bh_significant_0.05": str(q_value <= 0.05).lower(),
            "bonferroni_significant_0.05": str(
                p_value * site_count <= 0.05
            ).lower(),
        })

    fields = list(output[0])
    for path, selected in (
        (ALL_SITES, output),
        (
            SIGNIFICANT,
            [row for row in output if float(row["bh_q_value_901_sites"]) <= 0.05],
        ),
    ):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
            writer.writeheader()
            writer.writerows(selected)

    targeted = {9: "T12", 30: "G34", 107: "A114"}
    SUMMARY.write_text(json.dumps({
        "sequences": int(data["input"]["number of sequences"]),
        "codon_sites_tested": site_count,
        "sites_p_le_0.05": sum(p <= 0.05 for p in p_values),
        "sites_bh_q_le_0.05": sum(q <= 0.05 for q in bh_values),
        "sites_bonferroni_le_0.05": sum(
            p * site_count <= 0.05 for p in p_values
        ),
        "targeted_candidate_full_scan_results": {
            label: {
                "filtered_site": site,
                "p_value": p_values[site - 1],
                "bh_q_value": bh_values[site - 1],
                "bonferroni_p": min(1.0, p_values[site - 1] * site_count),
            }
            for site, label in targeted.items()
        },
        "interpretation": (
            "Discovery statistics use all 901 codons as the correction family. "
            "The earlier three-site MEME jobs remain targeted follow-up analyses."
        ),
    }, indent=2) + "\n")
    print(
        f"Summarized {site_count} sites; "
        f"{sum(q <= 0.05 for q in bh_values)} pass BH q<=0.05"
    )


if __name__ == "__main__":
    main()
