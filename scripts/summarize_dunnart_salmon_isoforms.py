#!/usr/bin/env python3
"""Summarize competitive Salmon estimates for dunnart DNMT3A models."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from summarize_dunnart_neocortex_stages import exact_permutation_p_greater


ROOT = Path(__file__).resolve().parents[1]
SALMON_DIR = ROOT / "results/06_isoform_evolution/salmon"
OUTPUT_DIR = ROOT / "results/06_isoform_evolution"
RUNS = {
    "SRR13036043": ("P12", "436"),
    "SRR13036044": ("P12", "451"),
    "SRR13036045": ("P12", "458"),
    "SRR13036046": ("P20", "666"),
    "SRR13036047": ("P20", "746"),
    "SRR13036048": ("P20", "758"),
}
FAMILIES = {
    "full_length_first_junction_family": {
        "XM_074300059.1", "XM_074300062.1", "XM_074300064.1",
    },
    "internal_066_067_family": {
        "XM_074300066.1", "XM_074300067.1",
    },
    "internal_068_family": {"XM_074300068.1"},
    "other_annotated_dnmt3a_models": {
        "XM_074300060.1", "XM_074300061.1",
        "XM_074300063.1", "XM_074300065.1",
    },
}


def family_for_model(model: str) -> str:
    matches = [family for family, models in FAMILIES.items() if model in models]
    if len(matches) != 1:
        raise ValueError(f"Model must belong to exactly one family: {model}")
    return matches[0]


def aggregate_families(rows: list[dict]) -> dict[str, dict[str, float]]:
    families = defaultdict(lambda: {"NumReads": 0.0, "TPM": 0.0})
    for row in rows:
        family = family_for_model(row["Name"])
        families[family]["NumReads"] += float(row["NumReads"])
        families[family]["TPM"] += float(row["TPM"])
    return dict(families)


def main() -> None:
    model_rows = []
    family_rows = []
    replicate_metrics = []
    for run, (stage, specimen) in RUNS.items():
        quant_path = SALMON_DIR / run / "quant.sf"
        meta_path = SALMON_DIR / run / "aux_info/meta_info.json"
        if not quant_path.exists() or not meta_path.exists():
            raise SystemExit(f"Missing Salmon output for {run}")
        with quant_path.open() as handle:
            quant_rows = list(csv.DictReader(handle, delimiter="\t"))
        meta = json.loads(meta_path.read_text())
        families = aggregate_families(quant_rows)
        total_assigned = sum(float(row["NumReads"]) for row in quant_rows)
        internal = families["internal_066_067_family"]
        full = families["full_length_first_junction_family"]

        for row in quant_rows:
            model_rows.append({
                "run_accession": run,
                "stage": stage,
                "specimen": specimen,
                "family": family_for_model(row["Name"]),
                "transcript_id": row["Name"],
                "length": row["Length"],
                "effective_length": row["EffectiveLength"],
                "tpm": row["TPM"],
                "estimated_fragments": row["NumReads"],
            })
        for family, values in families.items():
            family_rows.append({
                "run_accession": run,
                "stage": stage,
                "specimen": specimen,
                "family": family,
                "tpm": round(values["TPM"], 6),
                "estimated_fragments": round(values["NumReads"], 6),
                "fraction_of_assigned_dnmt3a_fragments": round(
                    values["NumReads"] / total_assigned, 6
                ),
            })
        replicate_metrics.append({
            "run_accession": run,
            "stage": stage,
            "specimen": specimen,
            "input_fragments": meta["num_processed"],
            "retained_dnmt3a_fragments": meta["num_mapped"],
            "decoy_fragments": meta["num_decoy_fragments"],
            "library_type": ",".join(meta["library_types"]),
            "internal_066_067_estimated_fragments": internal["NumReads"],
            "internal_066_067_fragments_per_million_input": (
                internal["NumReads"] * 1_000_000 / meta["num_processed"]
            ),
            "internal_066_067_tpm": internal["TPM"],
            "internal_066_067_fraction_of_assigned_dnmt3a": (
                internal["NumReads"] / total_assigned
            ),
            "internal_066_067_to_full_length_count_ratio": (
                internal["NumReads"] / full["NumReads"]
                if full["NumReads"] else None
            ),
        })

    for path, rows in (
        (OUTPUT_DIR / "dunnart_salmon_model_estimates.tsv", model_rows),
        (OUTPUT_DIR / "dunnart_salmon_family_estimates.tsv", family_rows),
        (OUTPUT_DIR / "dunnart_salmon_replicate_metrics.tsv", replicate_metrics),
    ):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    stages = {}
    for stage in ("P12", "P20"):
        rows = [row for row in replicate_metrics if row["stage"] == stage]
        total_input = sum(row["input_fragments"] for row in rows)
        total_target = sum(row["retained_dnmt3a_fragments"] for row in rows)
        total_internal = sum(
            row["internal_066_067_estimated_fragments"] for row in rows
        )
        stages[stage] = {
            "replicate_count": len(rows),
            "pooled_input_fragments": total_input,
            "pooled_retained_dnmt3a_fragments": total_target,
            "pooled_internal_066_067_estimated_fragments": round(
                total_internal, 6
            ),
            "pooled_internal_fragments_per_million_input": round(
                total_internal * 1_000_000 / total_input, 4
            ),
            "pooled_internal_fraction_of_assigned_dnmt3a": round(
                total_internal / total_target, 4
            ),
            "replicate_internal_fragments_per_million_range": [
                round(min(
                    row["internal_066_067_fragments_per_million_input"]
                    for row in rows
                ), 4),
                round(max(
                    row["internal_066_067_fragments_per_million_input"]
                    for row in rows
                ), 4),
            ],
            "replicate_internal_fraction_range": [
                round(min(
                    row["internal_066_067_fraction_of_assigned_dnmt3a"]
                    for row in rows
                ), 4),
                round(max(
                    row["internal_066_067_fraction_of_assigned_dnmt3a"]
                    for row in rows
                ), 4),
            ],
            "replicate_internal_to_full_count_ratio_range": [
                round(min(
                    row["internal_066_067_to_full_length_count_ratio"]
                    for row in rows
                ), 4),
                round(max(
                    row["internal_066_067_to_full_length_count_ratio"]
                    for row in rows
                ), 4),
            ],
        }

    contrast = {
        "P12_over_P20_internal_fragments_per_million_fold": round(
            stages["P12"]["pooled_internal_fragments_per_million_input"]
            / stages["P20"]["pooled_internal_fragments_per_million_input"],
            4,
        ),
        "P12_over_P20_internal_fraction_fold": round(
            stages["P12"]["pooled_internal_fraction_of_assigned_dnmt3a"]
            / stages["P20"]["pooled_internal_fraction_of_assigned_dnmt3a"],
            4,
        ),
        "all_P12_internal_rates_exceed_all_P20_rates": (
            stages["P12"]["replicate_internal_fragments_per_million_range"][0]
            > stages["P20"]["replicate_internal_fragments_per_million_range"][1]
        ),
        "all_P12_internal_fractions_exceed_all_P20_fractions": (
            stages["P12"]["replicate_internal_fraction_range"][0]
            > stages["P20"]["replicate_internal_fraction_range"][1]
        ),
        "internal_to_full_ratio_ranges_overlap": not (
            stages["P12"]["replicate_internal_to_full_count_ratio_range"][0]
            > stages["P20"]["replicate_internal_to_full_count_ratio_range"][1]
            or stages["P20"]["replicate_internal_to_full_count_ratio_range"][0]
            > stages["P12"]["replicate_internal_to_full_count_ratio_range"][1]
        ),
        "animal_level_exact_permutation_p_greater": {
            "internal_fragments_per_million_input": exact_permutation_p_greater(
                [
                    row["internal_066_067_fragments_per_million_input"]
                    for row in replicate_metrics if row["stage"] == "P12"
                ],
                [
                    row["internal_066_067_fragments_per_million_input"]
                    for row in replicate_metrics if row["stage"] == "P20"
                ],
            ),
            "internal_fraction_of_assigned_dnmt3a": exact_permutation_p_greater(
                [
                    row["internal_066_067_fraction_of_assigned_dnmt3a"]
                    for row in replicate_metrics if row["stage"] == "P12"
                ],
                [
                    row["internal_066_067_fraction_of_assigned_dnmt3a"]
                    for row in replicate_metrics if row["stage"] == "P20"
                ],
            ),
            "internal_to_full_length_count_ratio": exact_permutation_p_greater(
                [
                    row["internal_066_067_to_full_length_count_ratio"]
                    for row in replicate_metrics if row["stage"] == "P12"
                ],
                [
                    row["internal_066_067_to_full_length_count_ratio"]
                    for row in replicate_metrics if row["stage"] == "P20"
                ],
            ),
        },
    }
    summary = {
        "method": (
            "Salmon selective alignment against ten RefSeq DNMT3A targets "
            "with 2,093,960 non-DNMT3A Trinity transcripts as decoys; ISR "
            "paired libraries, concordant mappings only, sequence and GC bias "
            "correction."
        ),
        "stages": stages,
        "descriptive_contrast": contrast,
        "interpretation": (
            "Competitive transcript quantification supports higher pooled "
            "absolute and proportional internal-066/067 signal at P12. The "
            "absolute rate is higher in all three P12 animals, but the "
            "internal fraction and internal-to-full-length model ratio overlap "
            "between stages. The exact-junction result therefore remains the "
            "stronger evidence for a relative splice-architecture shift."
        ),
        "limits": [
            (
                f"Only {min(stages['P12']['replicate_count'], stages['P20']['replicate_count'])} "
                "biological replicates per stage are available."
            ),
            "Closely related transcript models create abundance uncertainty.",
            "Salmon model allocations do not replace direct junction counts.",
            "Decoy-filtered mapping rates are expectedly low for a "
            "single-locus target analysis.",
        ],
    }
    (
        OUTPUT_DIR / "dunnart_salmon_isoform_summary.json"
    ).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
