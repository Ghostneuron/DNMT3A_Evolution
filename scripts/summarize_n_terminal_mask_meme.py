#!/usr/bin/env python3
"""Summarize exact-site MEME tests on the N-terminal reliability mask."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/04_selection/mammals/reliability_mask"
MASK_DIR = ROOT / "results/02_alignment/reliability_mask"
SITES = {"T12": 9, "A16": 13, "G34": 30}


def holm_adjust(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [1.0] * len(values)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (len(values) - rank) * values[index]))
        adjusted[index] = running
    return adjusted


def result(label: str, site: int) -> dict[str, object]:
    data = json.loads((OUT / f"DNMT3A.{label}.MEME.json").read_text())
    headers = [
        item[0] if isinstance(item, list) else str(item)
        for item in data["MLE"]["headers"]
    ]
    values = data["MLE"]["content"]["0"][site - 1]
    record = dict(zip(headers, values))
    return {
        "human_site": label,
        "primary_filtered_site": site,
        "lrt": float(record["LRT"]),
        "meme_p": float(record["p-value"]),
        "branches_ebf_ge_100": int(record["# branches under selection"]),
    }


def main() -> None:
    rows = [result(label, site) for label, site in SITES.items()]
    adjusted = holm_adjust([float(row["meme_p"]) for row in rows])
    for row, value in zip(rows, adjusted):
        row["holm_p_3_mask_followups"] = value
        row["classification"] = (
            "nominal_only_not_holm_significant"
            if float(row["meme_p"]) <= 0.05 and value > 0.05
            else "unsupported"
        )

    with (MASK_DIR / "candidate_meme_comparison.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "tests": rows,
        "nominal_p_le_0.05": [
            row["human_site"] for row in rows if float(row["meme_p"]) <= 0.05
        ],
        "holm_significant_3_tests": [
            row["human_site"]
            for row in rows
            if float(row["holm_p_3_mask_followups"]) <= 0.05
        ],
        "interpretation": (
            "T12 and G34 retain nominal exact-site MEME support after the "
            "taxon-wide reliability mask, but neither passes Holm correction "
            "across T12, A16, and G34. A16 is unsupported."
        ),
    }
    (MASK_DIR / "meme_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
