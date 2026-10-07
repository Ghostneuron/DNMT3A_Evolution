#!/usr/bin/env python3
"""Overlay retained candidate states on available expression trajectories."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from dnmt3a_pipeline import CODON_TABLE, read_fasta


ROOT = Path(__file__).resolve().parents[1]
ALIGNMENT = (
    ROOT
    / "results/02_alignment/sensitivity/mammal_scope/macse"
    / "DNMT3A_MACSE_mammal_scope_HyPhy_unique.fasta"
)
PRANK = ROOT / "results/02_alignment/sensitivity/prank/DNMT3A_PRANK_clean_mammals_HyPhy_unique.fasta"
PRIMARY_CROSSWALK = (
    ROOT
    / "results/02_alignment/sensitivity/mammal_scope/macse"
    / "MACSE_mammal_scope_crosswalk.tsv"
)
PRANK_CROSSWALK = ROOT / "results/02_alignment/sensitivity/prank/PRANK_human_coordinate_crosswalk.tsv"
TAXA = ROOT / "results/01_curated/taxon_metadata.tsv"
TRAJECTORIES = ROOT / "results/05_brain_integration/expression_trajectory_monotonicity.tsv"
RESULTS = ROOT / "results/05_brain_integration"
HUMAN_POSITIONS = {"T12": 12, "G34": 34, "S97": 97}


def alignment_sites(path: Path) -> dict[str, int]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_human = {
        int(row["human_DNMT3A1_aa"]): int(row["filtered_codon_site_1based"])
        for row in rows if row["human_DNMT3A1_aa"]
    }
    return {label: by_human[position] for label, position in HUMAN_POSITIONS.items()}


def residues(path: Path, sites: dict[str, int]) -> dict[str, dict[str, str]]:
    output = {}
    for record in read_fasta(path):
        output[record.identifier] = {
            label: CODON_TABLE.get(
                record.sequence[(site - 1) * 3:site * 3], "-"
            )
            for label, site in sites.items()
        }
    return output


def main() -> None:
    with TAXA.open(newline="") as handle:
        safe_by_species = {
            row["species"]: row["safe_id"]
            for row in csv.DictReader(handle, delimiter="\t")
        }
    with TRAJECTORIES.open(newline="") as handle:
        trajectories = list(csv.DictReader(handle, delimiter="\t"))
    primary = residues(ALIGNMENT, alignment_sites(PRIMARY_CROSSWALK))
    prank = residues(PRANK, alignment_sites(PRANK_CROSSWALK))

    by_species: dict[str, dict[str, str]] = {}
    for row in trajectories:
        species = row["species"]
        region_key = (
            "cerebrum_rho" if row["brain_region"] == "brain/cerebrum series"
            else "cerebellum_rho"
        )
        by_species.setdefault(species, {})[region_key] = (
            row["spearman_rho_stage_vs_median_log2_CPM"]
        )

    output = []
    for species in sorted(by_species):
        safe_id = safe_by_species[species]
        output.append({
            "species": species,
            "safe_id": safe_id,
            **{
                f"mammal_scope_MACSE_{label}": primary[safe_id][label]
                for label in HUMAN_POSITIONS
            },
            **{
                f"PRANK_{label}": prank[safe_id][label]
                for label in HUMAN_POSITIONS
            },
            "cerebrum_trajectory_rho": by_species[species].get("cerebrum_rho", ""),
            "cerebellum_trajectory_rho": by_species[species].get("cerebellum_rho", ""),
            "candidate_association_status": (
                "descriptive_only; insufficient independent residue-state replication"
            ),
        })

    fields = [
        "species", "safe_id",
        *[
            field
            for label in HUMAN_POSITIONS
            for field in (f"mammal_scope_MACSE_{label}", f"PRANK_{label}")
        ],
        "cerebrum_trajectory_rho",
        "cerebellum_trajectory_rho", "candidate_association_status",
    ]
    with (RESULTS / "candidate_expression_overlay.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)

    state_counts = {
        label: dict(sorted(Counter(
            row[f"mammal_scope_MACSE_{label}"] for row in output
        ).items()))
        for label in HUMAN_POSITIONS
    }
    summary = {
        "expression_species": len(output),
        **{
            f"{label}_state_counts": counts
            for label, counts in state_counts.items()
        },
        "mammal_scope_MACSE_PRANK_state_disagreements": sum(
            any(
                row[f"mammal_scope_MACSE_{label}"] != row[f"PRANK_{label}"]
                for label in HUMAN_POSITIONS
            )
            for row in output
        ),
        "candidate_expression_association_testable": False,
        "reason": (
            "G34 and S97 are invariant and the non-reference T12 state occurs "
            "only in opossum among the six expression species"
        ),
    }
    (RESULTS / "candidate_expression_overlay_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
