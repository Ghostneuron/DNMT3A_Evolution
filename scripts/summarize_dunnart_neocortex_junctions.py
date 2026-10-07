#!/usr/bin/env python3
"""Summarize replicated dunnart developmental-neocortex junction scans."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/06_isoform_evolution"
INPUTS = (
    RESULTS / "dunnart_P20_neocortex_rep1_junction_summary.json",
    RESULTS / "dunnart_P20_neocortex_rep2_junction_summary.json",
    RESULTS / "dunnart_P20_neocortex_rep3_junction_summary.json",
)
COMPONENTS = (
    "full_length_first_junction",
    "internal_067_first_junction",
    "internal_068_first_junction",
    "common_downstream_junction",
)


def summarize(records: list[dict]) -> tuple[list[dict], dict]:
    rows = []
    pooled_counts = Counter()
    pooled_pairs = 0
    for record in records:
        pair_count = record["paired_read_count"]
        counts = record["junction_supporting_pair_counts"]
        diagnostic_total = (
            counts["full_length_first_junction"]
            + counts["internal_067_first_junction"]
            + counts["internal_068_first_junction"]
        )
        row = {
            "run_accession": record["run_accession"],
            "sample": record["sample"],
            "paired_read_count": pair_count,
        }
        for component in COMPONENTS:
            count = counts[component]
            row[f"{component}_pairs"] = count
            row[f"{component}_pairs_per_million"] = round(
                count * 1_000_000 / pair_count, 4
            )
            pooled_counts[component] += count
        row["internal_067_fraction_of_diagnostic_first_junctions"] = round(
            counts["internal_067_first_junction"] / diagnostic_total, 4
        ) if diagnostic_total else None
        rows.append(row)
        pooled_pairs += pair_count

    pooled_diagnostic = (
        pooled_counts["full_length_first_junction"]
        + pooled_counts["internal_067_first_junction"]
        + pooled_counts["internal_068_first_junction"]
    )
    pooled = {
        "replicate_count": len(records),
        "pooled_paired_read_count": pooled_pairs,
        "pooled_junction_supporting_pair_counts": {
            component: pooled_counts[component]
            for component in COMPONENTS
        },
        "pooled_pairs_per_million": {
            component: round(
                pooled_counts[component] * 1_000_000 / pooled_pairs, 4
            )
            for component in COMPONENTS
        },
        "internal_067_fraction_of_diagnostic_first_junctions": round(
            pooled_counts["internal_067_first_junction"] / pooled_diagnostic, 4
        ) if pooled_diagnostic else None,
        "replication_statement": (
            "The internal XM_074300066.1/XM_074300067.1 first-junction "
            "architecture is supported by exact junction-spanning read pairs "
            f"in {len(records)} independent P20 neocortex animals."
        ),
        "limits": [
            "The reads validate a splice architecture, not the exact capped TSS.",
            "The shared junction does not distinguish XM_074300066.1 from "
            "XM_074300067.1.",
            "P20 replication alone does not establish developmental-stage "
            "regulation; use the matched P12/P20 stage summary for that test.",
            "Sparse diagnostic-junction counts should not be interpreted as "
            "whole-transcript abundance estimates.",
        ],
    }
    return rows, pooled


def main() -> None:
    missing = [str(path) for path in INPUTS if not path.exists()]
    if missing:
        raise SystemExit(f"Missing replicate summaries: {missing}")
    records = [json.loads(path.read_text()) for path in INPUTS]
    rows, pooled = summarize(records)

    output_tsv = RESULTS / "dunnart_P20_neocortex_junction_replicates.tsv"
    with output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    (
        RESULTS / "dunnart_P20_neocortex_junction_replicates_summary.json"
    ).write_text(json.dumps(pooled, indent=2) + "\n")
    print(json.dumps(pooled, indent=2))


if __name__ == "__main__":
    main()
