#!/usr/bin/env python3
"""Design a lineage-supported DNMT3A1 functional substitution panel."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MASKED = ROOT / "results/04_selection/mammals/reliability_mask/lineages/prioritized_change_coupled_branches.tsv"
S97 = ROOT / "results/04_selection/mammals/s97_homology_audit/lineages/prioritized_change_coupled_branches.tsv"
OUT = ROOT / "results/07_functional_prioritization"
ATLAS_SUMMARY = (
    ROOT
    / "results/05_brain_integration/human_cell_atlas_import_summary.json"
)
HUMAN = {"T12": "T", "G34": "G", "S97": "S"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    events = read_tsv(MASKED) + read_tsv(S97)
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in events:
        site = row["human_site"]
        derived = row["derived_aa"]
        if site in HUMAN and derived != HUMAN[site]:
            grouped[(site, derived)].append(row)

    rows: list[dict[str, object]] = []
    for (site, derived), members in grouped.items():
        internal = sum(row["branch_type"] == "internal" for row in members)
        terminal = len(members) - internal
        if internal >= 1 and len(members) >= 2:
            tier = 1
            rationale = "recurrent_with_internal_lineage_support"
        elif internal >= 1 or len(members) >= 2:
            tier = 2
            rationale = "internal_or_recurrent_lineage_support"
        else:
            tier = 3
            rationale = "single_terminal_change"
        rows.append({
            "priority_tier": tier,
            "human_site": site,
            "human_reference_residue": HUMAN[site],
            "experimental_substitution": f"{site}{derived}",
            "derived_residue": derived,
            "change_event_count": len(members),
            "internal_branch_events": internal,
            "terminal_branch_events": terminal,
            "maximum_MEME_positive_class_posterior": max(
                float(row["positive_class_posterior"]) for row in members
            ),
            "supporting_branches": ";".join(row["branch"] for row in members),
            "supporting_lineages": ";".join(row["branch_label"] for row in members),
            "evolutionary_rationale": rationale,
            "construct_background": "human_DNMT3A1",
            "primary_assay_class": (
                "protein_abundance;chromatin_localization;"
                "tail_interaction_profile;methylation_activity"
            ),
            "brain_claim_allowed": "false",
        })
    rows.sort(key=lambda row: (
        int(row["priority_tier"]),
        {"G34": 1, "T12": 2, "S97": 3}[str(row["human_site"])],
        -int(row["change_event_count"]),
        str(row["experimental_substitution"]),
    ))
    with (OUT / "DNMT3A1_variant_panel.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    controls = [
        {
            "control": "human_DNMT3A1_WT",
            "purpose": "reference for every substitution",
        },
        {
            "control": "human_DNMT3A2",
            "purpose": "control lacking the DNMT3A1-specific N-terminal tail",
        },
        {
            "control": "matched_empty_vector",
            "purpose": "background assay and localization control",
        },
    ]
    with (OUT / "control_panel.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(controls[0]))
        writer.writeheader()
        writer.writerows(controls)

    atlas = json.loads(ATLAS_SUMMARY.read_text())
    neuronal_mch = float(atlas["motor_cortex_neuronal_mCH_median"])
    non_neuronal_mch = float(atlas["motor_cortex_non_neuronal_mCH_median"])
    donor_ratios = atlas[
        "motor_cortex_donor_neuronal_to_non_neuronal_mCH_median_ratio"
    ]
    assay_context = [
        {
            "context": "adult_human_primary_motor_cortex_neuronal",
            "metric": "median_mCH_fraction",
            "value": neuronal_mch,
            "unit": "fraction",
            "biological_replicates_in_reference": len(donor_ratios),
            "experimental_role": "later_stage_neuronal_benchmark",
            "interpretive_limit": (
                "reference distribution; not a candidate-variant effect"
            ),
        },
        {
            "context": "adult_human_primary_motor_cortex_non_neuronal",
            "metric": "median_mCH_fraction",
            "value": non_neuronal_mch,
            "unit": "fraction",
            "biological_replicates_in_reference": len(donor_ratios),
            "experimental_role": "matched_lineage_context_control",
            "interpretive_limit": (
                "reference distribution; not a candidate-variant effect"
            ),
        },
        {
            "context": "neuronal_vs_non_neuronal",
            "metric": "median_mCH_ratio",
            "value": neuronal_mch / non_neuronal_mch,
            "unit": "fold",
            "biological_replicates_in_reference": len(donor_ratios),
            "experimental_role": "assay_dynamic_range_reference",
            "interpretive_limit": (
                "single cells are not independent biological replicates"
            ),
        },
    ]
    with (OUT / "human_motor_cortex_assay_context.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=list(assay_context[0])
        )
        writer.writeheader()
        writer.writerows(assay_context)

    summary = {
        "tier_1_substitutions": [
            row["experimental_substitution"]
            for row in rows if row["priority_tier"] == 1
        ],
        "tier_2_substitutions": [
            row["experimental_substitution"]
            for row in rows if row["priority_tier"] == 2
        ],
        "tier_3_substitutions": [
            row["experimental_substitution"]
            for row in rows if row["priority_tier"] == 3
        ],
        "recommended_first_pass": (
            "Test all tier-1 substitutions individually in human DNMT3A1 "
            "with WT, DNMT3A2, and empty-vector controls. Add tier 2 only "
            "after a reproducible molecular phenotype."
        ),
        "interpretive_boundary": (
            "These assays can test molecular consequences of evolutionary "
            "substitutions; neuronal or developmental experiments are "
            "required for brain-development inference."
        ),
        "human_motor_cortex_assay_context": {
            "neuronal_mCH_median": neuronal_mch,
            "non_neuronal_mCH_median": non_neuronal_mch,
            "neuronal_to_non_neuronal_median_ratio": (
                neuronal_mch / non_neuronal_mch
            ),
            "reference_donors": len(donor_ratios),
            "required_design_boundary": (
                "Use biological replicates in matched neuronal and "
                "non-neuronal contexts; do not use cells as replicate "
                "organisms."
            ),
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
