#!/usr/bin/env python3
"""Summarize the corrected 909-site mammal-scope MACSE MEME scan."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/04_selection/mammals/mammal_scope/full_meme"
RESULT = OUT / "DNMT3A.MACSE.mammal_scope.full.MEME.json"
CROSSWALK = ROOT / "results/02_alignment/sensitivity/mammal_scope/macse/MACSE_mammal_scope_crosswalk.tsv"
OLD_ALL = ROOT / "results/04_selection/mammals/full_meme/all_sites.tsv"
CANDIDATES = ("G34", "T12", "S97")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(rows[0]) if rows else []
        )
        if rows:
            writer.writeheader()
            writer.writerows(rows)


def bh_adjust(values: list[float]) -> list[float]:
    count = len(values)
    order = sorted(range(count), key=values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank, index in reversed(list(enumerate(order, start=1))):
        running = min(running, values[index] * count / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def index(headers: list[str], token: str) -> int:
    return next(
        position for position, header in enumerate(headers)
        if token.casefold() in header.casefold()
    )


def main() -> None:
    data = json.loads(RESULT.read_text())
    values = data["MLE"]["content"]["0"]
    if len(values) != 909:
        raise RuntimeError(f"Expected 909 MEME rows, observed {len(values)}")
    headers = [
        item[0] if isinstance(item, list) else str(item)
        for item in data["MLE"]["headers"]
    ]
    p_index = index(headers, "p-value")
    lrt_index = index(headers, "LRT")
    branch_index = index(headers, "# branches")
    mapping = {
        int(row["filtered_codon_site_1based"]): row
        for row in read_tsv(CROSSWALK)
    }
    p_values = [float(row[p_index]) for row in values]
    q_values = bh_adjust(p_values)
    rows: list[dict[str, object]] = []
    for site, (fit, p_value, q_value) in enumerate(
        zip(values, p_values, q_values), start=1
    ):
        mapped = mapping[site]
        rows.append({
            "filtered_codon_site_1based": site,
            "human_site": mapped["human_residue"] + mapped["human_DNMT3A1_aa"],
            "human_DNMT3A1_aa": mapped["human_DNMT3A1_aa"],
            "human_residue": mapped["human_residue"],
            "domain": mapped["domain"],
            "mammal_occupancy": mapped["mammal_occupancy"],
            "LRT": fit[lrt_index],
            "p_value": p_value,
            "bh_q_value_909_sites": q_value,
            "bonferroni_p_909_sites": min(1.0, p_value * 909),
            "branches_EBF_ge_100": fit[branch_index],
            "bh_significant_0.05": str(q_value <= 0.05).lower(),
        })
    write_tsv(OUT / "all_sites.tsv", rows)
    write_tsv(
        OUT / "significant_sites.tsv",
        [row for row in rows if float(row["bh_q_value_909_sites"]) <= 0.05],
    )

    by_label = {str(row["human_site"]): row for row in rows}
    old = {
        row["human_residue"] + row["human_DNMT3A1_aa"]: row
        for row in read_tsv(OLD_ALL)
    }
    comparison = []
    for label in sorted(set(old) & set(by_label)):
        comparison.append({
            "human_site": label,
            "old_901_p": old[label]["p_value"],
            "mammal_scope_909_p": by_label[label]["p_value"],
            "old_901_bh_q": old[label]["bh_q_value_901_sites"],
            "mammal_scope_909_bh_q": by_label[label]["bh_q_value_909_sites"],
        })
    write_tsv(OUT / "common_site_comparison.tsv", comparison)

    restored = {"D22", "S97", "P106", "A107", "A116", "E117", "Q231", "G232"}
    write_tsv(
        OUT / "restored_site_results.tsv",
        [row for row in rows if row["human_site"] in restored],
    )
    summary = {
        "sequences": int(data["input"]["number of sequences"]),
        "codon_sites_tested": len(rows),
        "sites_p_le_0.05": sum(value <= 0.05 for value in p_values),
        "sites_bh_q_le_0.05": sum(value <= 0.05 for value in q_values),
        "sites_bonferroni_le_0.05": sum(
            value * len(rows) <= 0.05 for value in p_values
        ),
        "retained_candidate_results": {
            label: {
                "filtered_site": by_label[label]["filtered_codon_site_1based"],
                "p_value": by_label[label]["p_value"],
                "bh_q_value_909_sites": by_label[label]["bh_q_value_909_sites"],
                "bonferroni_p_909_sites": by_label[label][
                    "bonferroni_p_909_sites"
                ],
            }
            for label in CANDIDATES
        },
        "restored_site_results": {
            label: {
                "p_value": by_label[label]["p_value"],
                "bh_q_value_909_sites": by_label[label]["bh_q_value_909_sites"],
            }
            for label in sorted(restored)
        },
        "interpretation": (
            "Discovery statistics are corrected across all 909 codons after "
            "occupancy filtering within the 250-mammal analysis scope."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
