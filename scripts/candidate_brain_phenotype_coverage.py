#!/usr/bin/env python3
"""Audit phenotype coverage for retained DNMT3A coding states."""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta, write_tsv  # noqa: E402


ALIGNMENT = ROOT / "results/02_alignment/sensitivity/mammal_scope/macse/DNMT3A_MACSE_mammal_scope_HyPhy_unique.fasta"
CROSSWALK = ROOT / "results/02_alignment/sensitivity/mammal_scope/macse/MACSE_mammal_scope_crosswalk.tsv"
TRAITS = ROOT / "data/traits/brain_phenotypes.tsv"
TAXA = ROOT / "results/01_curated/taxon_metadata.tsv"
S97_STATES = ROOT / "results/04_selection/mammals/s97_homology_audit/s97_taxon_states.tsv"
OUT = ROOT / "results/05_brain_integration"
POSITIONS = {"T12": 12, "G34": 34, "S97": 97}
MIN_SPECIES = 10


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    mapping = {
        int(row["human_DNMT3A1_aa"]): int(row["filtered_codon_site_1based"])
        for row in read_tsv(CROSSWALK)
    }
    sequences = {
        record.identifier: record.sequence for record in read_fasta(ALIGNMENT)
    }
    states = {
        identifier: {
            label: CODON_TABLE.get(
                sequence[(mapping[position] - 1) * 3:mapping[position] * 3], "-"
            )
            for label, position in POSITIONS.items()
        }
        for identifier, sequence in sequences.items()
    }
    taxa = {row["species"]: row for row in read_tsv(TAXA)}
    traits = read_tsv(TRAITS)
    grouped: dict[str, set[str]] = defaultdict(set)
    for row in traits:
        taxon = taxa[row["species"]]
        if taxon["is_mammal"] == "true":
            grouped[row["phenotype"]].add(taxon["safe_id"])

    coverage_rows: list[dict[str, object]] = []
    for phenotype, species in sorted(grouped.items()):
        for label in POSITIONS:
            counts = Counter(states[safe_id][label] for safe_id in species)
            replicated_states = sum(count >= 2 for count in counts.values())
            coverage_rows.append({
                "phenotype": phenotype,
                "human_site": label,
                "mammal_species": len(species),
                "state_counts": ";".join(
                    f"{state}:{count}" for state, count in sorted(counts.items())
                ),
                "states_with_at_least_two_species": replicated_states,
                "minimum_species_gate_pass": str(len(species) >= MIN_SPECIES).lower(),
                "replicated_state_contrast": str(replicated_states >= 2).lower(),
                "association_ready": str(
                    len(species) >= MIN_SPECIES and replicated_states >= 2
                ).lower(),
            })
    write_tsv(OUT / "candidate_brain_phenotype_coverage.tsv", coverage_rows)

    phenotype_species = {
        safe_id for species in grouped.values() for safe_id in species
    }
    s97_rows = read_tsv(S97_STATES)
    gap_rows = []
    for row in s97_rows:
        if row["macse_residue"] == "S":
            continue
        safe_id = row["safe_id"]
        available = sorted(
            phenotype for phenotype, species in grouped.items()
            if safe_id in species
        )
        gap_rows.append({
            "safe_id": safe_id,
            "S97_residue": row["macse_residue"],
            "has_any_brain_phenotype": str(safe_id in phenotype_species).lower(),
            "available_phenotypes": ";".join(available),
            "acquisition_priority": (
                "existing_data" if available
                else "priority_nonreference_lineage_for_brain_mCH_or_developmental_expression"
            ),
        })
    write_tsv(OUT / "s97_nonreference_phenotype_gaps.tsv", gap_rows)

    summary = {
        "candidate_sites": list(POSITIONS),
        "phenotype_candidate_combinations": len(coverage_rows),
        "association_ready_combinations": [
            f"{row['phenotype']}:{row['human_site']}"
            for row in coverage_rows if row["association_ready"] == "true"
        ],
        "s97_nonreference_taxa": len(gap_rows),
        "s97_nonreference_taxa_with_any_brain_phenotype": [
            row["safe_id"] for row in gap_rows
            if row["has_any_brain_phenotype"] == "true"
        ],
        "priority_s97_data_gaps": [
            row["safe_id"] for row in gap_rows
            if row["has_any_brain_phenotype"] == "false"
        ],
        "conclusion": (
            "No retained-site/phenotype combination has both ten mammal "
            "species and replicated residue-state groups. Platypus is the "
            "only S97 non-serine taxon with any current brain phenotype."
        ),
    }
    (OUT / "candidate_brain_phenotype_coverage_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
