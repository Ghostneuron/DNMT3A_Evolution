#!/usr/bin/env python3
"""Summarize the targeted G34 aBSREL run with transparent multiplicity control."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/04_selection/mammals/g34_absrel"
MANIFEST = OUTDIR / "foreground_manifest.tsv"
RESULT = OUTDIR / "DNMT3A.G34_foreground.aBSREL.json"
SUMMARY = OUTDIR / "branch_summary.tsv"
SUMMARY_JSON = OUTDIR / "summary.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def holm_adjust(p_values: list[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [0.0] * len(p_values)
    running = 0.0
    count = len(p_values)
    for rank, index in enumerate(order):
        value = min(1.0, (count - rank) * p_values[index])
        running = max(running, value)
        adjusted[index] = running
    return adjusted


def main() -> None:
    foreground = read_tsv(MANIFEST)
    data = json.loads(RESULT.read_text())
    attributes = data["branch attributes"]["0"]
    rows: list[dict[str, object]] = []
    for selected in foreground:
        branch = selected["branch"]
        if branch not in attributes:
            raise RuntimeError(f"aBSREL output is missing foreground branch {branch}")
        values = attributes[branch]
        tested = str(values.get("original name", "")).endswith("{Foreground}")
        if not tested:
            # HyPhy versions differ in whether they retain the annotation in
            # "original name"; presence of an LRT is the definitive check.
            tested = "LRT" in values and "Uncorrected P-value" in values
        if not tested:
            raise RuntimeError(f"aBSREL did not test {branch}")
        rows.append({
            "branch": branch,
            "descendant_tips": selected["descendant_tips"],
            "reconstructed_G34_change": selected["reconstructed_change"],
            "LRT": values.get("LRT", ""),
            "uncorrected_p_value": values.get("Uncorrected P-value", ""),
            "hyphy_corrected_p_value": values.get("Corrected P-value", ""),
            "holm_p_value_two_foregrounds": "",
            "significant_holm_0.05": "",
            "rate_distributions": json.dumps(
                values.get("Rate Distributions", []), separators=(",", ":")
            ),
            "interpretation_scope": (
                "post_hoc_branch_wide_concordance;not_independent_G34_validation;"
                "not_brain_phenotype_evidence"
            ),
        })
    p_values = [float(row["uncorrected_p_value"]) for row in rows]
    for row, adjusted in zip(rows, holm_adjust(p_values)):
        row["holm_p_value_two_foregrounds"] = adjusted
        row["significant_holm_0.05"] = str(adjusted <= 0.05).lower()

    with SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    significant = sum(row["significant_holm_0.05"] == "true" for row in rows)
    SUMMARY_JSON.write_text(json.dumps({
        "foreground_branches_tested": len(rows),
        "holm_significant_branches_0.05": significant,
        "multiplicity_control": "Holm across the two post-hoc G34 foreground branches",
        "interpretation": (
            "aBSREL is a whole-branch episodic-selection test. Because foregrounds "
            "were chosen after the G34 MEME result, this is a concordance analysis, "
            "not an independent confirmation of the site or a brain-development claim."
        ),
    }, indent=2) + "\n")
    print(f"Summarized {len(rows)} branches; {significant} Holm-significant at 0.05")


if __name__ == "__main__":
    main()
