#!/usr/bin/env python3
"""Profile evolution of functional subregions within the DNMT3A1 N terminus."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/04_selection/mammals/site_domain_metrics.tsv"
CANDIDATES = ROOT / "results/04_selection/mammals/final_candidate_evidence.tsv"
REGIONS = ROOT / "config/human_regulatory_subregions.tsv"
OUT = ROOT / "results/04_selection/mammals/regulatory_tail"
NEGATIVE_POSTERIOR_THRESHOLD = 0.90


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def hypergeometric_lower_tail(
    population: int, successes: int, draws: int, observed: int
) -> float:
    """P(X <= observed) for a hypergeometric random variable."""
    denominator = math.comb(population, draws)
    lower = max(0, draws - (population - successes))
    upper = min(observed, successes, draws)
    return sum(
        math.comb(successes, value)
        * math.comb(population - successes, draws - value)
        for value in range(lower, upper + 1)
    ) / denominator


def odds_ratio(
    focal_positive: int,
    focal_negative: int,
    reference_positive: int,
    reference_negative: int,
) -> float:
    if focal_negative == 0 or reference_positive == 0:
        return float("inf")
    return (focal_positive * reference_negative) / (
        focal_negative * reference_positive
    )


def main() -> None:
    sites = read_tsv(INPUT)
    regions = read_tsv(REGIONS)
    if not sites:
        raise RuntimeError(f"No site metrics found in {INPUT}")

    site_positions = {int(row["human_DNMT3A1_aa"]) for row in sites}
    region_by_position: dict[int, dict[str, str]] = {}
    for region in regions:
        start = int(region["start_aa"])
        end = int(region["end_aa"])
        for position in range(start, end + 1):
            if position in region_by_position:
                raise RuntimeError(f"Overlapping subregion definition at {position}")
            region_by_position[position] = region
    if not site_positions <= set(region_by_position):
        missing = sorted(site_positions - set(region_by_position))
        raise RuntimeError(f"Retained sites lack a subregion: {missing[:10]}")

    grouped: dict[str, list[dict[str, str]]] = {
        row["region"]: [] for row in regions
    }
    for site in sites:
        region = region_by_position[int(site["human_DNMT3A1_aa"])]["region"]
        grouped[region].append(site)

    summary_rows: list[dict[str, object]] = []
    for region in regions:
        name = region["region"]
        rows = grouped[name]
        if not rows:
            continue
        polymorphic = sum(row["polymorphic"] == "True" for row in rows)
        recurrent = sum(row["recurrently_variable"] == "True" for row in rows)
        constrained = sum(
            float(row["primary_negative_posterior"])
            >= NEGATIVE_POSTERIOR_THRESHOLD
            for row in rows
        )
        summary_rows.append({
            "region": name,
            "start_aa": int(region["start_aa"]),
            "end_aa": int(region["end_aa"]),
            "evidence_class": region["evidence_class"],
            "retained_sites": len(rows),
            "mean_mammal_unique_occupancy": statistics.fmean(
                float(row["mammal_unique_occupancy"]) for row in rows
            ),
            "polymorphic_sites": polymorphic,
            "polymorphic_fraction": polymorphic / len(rows),
            "recurrently_variable_sites": recurrent,
            "recurrently_variable_fraction": recurrent / len(rows),
            "mean_normalized_amino_acid_entropy": statistics.fmean(
                float(row["normalized_amino_acid_entropy"]) for row in rows
            ),
            "negative_selection_sites_posterior_ge_0.90": constrained,
            "negative_selection_fraction_posterior_ge_0.90": (
                constrained / len(rows)
            ),
            "primary_fubar_candidates": sum(
                row["primary_candidate"] == "True" for row in rows
            ),
            "note": region["note"],
        })

    upstream = grouped["upstream_disordered_tail"]
    engagement = grouped["nucleosome_H2AK119ub_engagement"]
    comparison_rows: list[dict[str, object]] = []
    combined = upstream + engagement
    for label, field in (
        ("polymorphic", "polymorphic"),
        ("recurrently_variable", "recurrently_variable"),
    ):
        upstream_positive = sum(row[field] == "True" for row in upstream)
        engagement_positive = sum(row[field] == "True" for row in engagement)
        total_positive = upstream_positive + engagement_positive
        comparison_rows.append({
            "comparison": (
                "nucleosome_H2AK119ub_engagement_vs_upstream_disordered_tail"
            ),
            "variability_definition": label,
            "engagement_positive_sites": engagement_positive,
            "engagement_retained_sites": len(engagement),
            "engagement_positive_fraction": engagement_positive / len(engagement),
            "upstream_positive_sites": upstream_positive,
            "upstream_retained_sites": len(upstream),
            "upstream_positive_fraction": upstream_positive / len(upstream),
            "engagement_vs_upstream_odds_ratio": odds_ratio(
                engagement_positive,
                len(engagement) - engagement_positive,
                upstream_positive,
                len(upstream) - upstream_positive,
            ),
            "one_sided_hypergeometric_depletion_p": hypergeometric_lower_tail(
                len(combined),
                total_positive,
                len(engagement),
                engagement_positive,
            ),
        })

    candidate_rows: list[dict[str, object]] = []
    for candidate in read_tsv(CANDIDATES):
        if candidate["human_site"] not in {"T12", "G34", "S97"}:
            continue
        position = int(candidate["human_DNMT3A1_aa"])
        region = region_by_position[position]
        candidate_rows.append({
            "priority_rank": int(candidate["priority_rank"]),
            "human_site": candidate["human_site"],
            "human_DNMT3A1_aa": position,
            "regulatory_subregion": region["region"],
            "inside_supported_164_219_engagement_region": (
                region["region"] == "nucleosome_H2AK119ub_engagement"
            ),
            "mammal_scope_MEME_BH_q_909": candidate[
                "mammal_scope_full_meme_bh_q_909"
            ],
            "evidence_class": candidate["evidence_class"],
            "functional_status": (
                "evolutionary_follow-up hypothesis; no demonstrated "
                "molecular or brain phenotype"
            ),
            "recommended_role": (
                "secondary experimental candidate, not a central "
                "mechanistic claim"
            ),
        })
    candidate_rows.sort(key=lambda row: int(row["priority_rank"]))

    known_region_rows = []
    for row in engagement:
        known_region_rows.append({
            "human_DNMT3A1_aa": int(row["human_DNMT3A1_aa"]),
            "human_residue": row["human_residue"],
            "mammal_unique_occupancy": row["mammal_unique_occupancy"],
            "amino_acid_states": int(row["amino_acid_states"]),
            "nonmajor_residue_count": int(row["nonmajor_residue_count"]),
            "polymorphic": row["polymorphic"],
            "recurrently_variable": row["recurrently_variable"],
            "normalized_amino_acid_entropy": row[
                "normalized_amino_acid_entropy"
            ],
            "primary_negative_posterior": row["primary_negative_posterior"],
            "primary_positive_posterior": row["primary_positive_posterior"],
        })
    known_region_rows.sort(
        key=lambda row: float(row["normalized_amino_acid_entropy"]),
        reverse=True,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT / "regulatory_subregion_summary.tsv", summary_rows)
    write_tsv(OUT / "regulatory_subregion_comparisons.tsv", comparison_rows)
    write_tsv(OUT / "candidate_region_context.tsv", candidate_rows)
    write_tsv(OUT / "engagement_region_site_metrics.tsv", known_region_rows)

    by_name = {str(row["region"]): row for row in summary_rows}
    summary = {
        "coordinate_system": "human_DNMT3A1_912aa",
        "mammal_unique_sequences": 250,
        "retained_sites": len(sites),
        "engagement_region": {
            "start_aa": 164,
            "end_aa": 219,
            "retained_sites": by_name[
                "nucleosome_H2AK119ub_engagement"
            ]["retained_sites"],
            "polymorphic_sites": by_name[
                "nucleosome_H2AK119ub_engagement"
            ]["polymorphic_sites"],
            "recurrently_variable_sites": by_name[
                "nucleosome_H2AK119ub_engagement"
            ]["recurrently_variable_sites"],
            "mean_normalized_amino_acid_entropy": by_name[
                "nucleosome_H2AK119ub_engagement"
            ]["mean_normalized_amino_acid_entropy"],
        },
        "candidate_sites": [row["human_site"] for row in candidate_rows],
        "candidates_inside_engagement_region": [
            row["human_site"]
            for row in candidate_rows
            if row["inside_supported_164_219_engagement_region"]
        ],
        "interpretation": (
            "The experimentally supported 164-219 engagement region is "
            "strongly conserved relative to residues 1-163. T12, G34, and "
            "S97 lie in the variable upstream tail and remain secondary "
            "functional hypotheses."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
