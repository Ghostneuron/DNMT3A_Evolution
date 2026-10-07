#!/usr/bin/env python3
"""Map retained DNMT3A candidates to UniProt and AlphaFold annotations."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from Bio.PDB import PDBParser


ROOT = Path(__file__).resolve().parents[1]
UNIPROT = ROOT / "data/external/uniprot/Q9Y6K1.json"
ALPHAFOLD = ROOT / "data/external/alphafold/AF-Q9Y6K1-F1-model_v6.pdb"
FINAL = ROOT / "results/04_selection/mammals/final_candidate_evidence.tsv"
OUT = ROOT / "results/04_selection/mammals/candidate_biological_context"
RETAINED = ("G34", "T12", "S97")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    entry = json.loads(UNIPROT.read_text())
    sequence = entry["sequence"]["value"]
    features = entry["features"]
    structure = PDBParser(QUIET=True).get_structure("DNMT3A", ALPHAFOLD)
    chain = next(structure.get_chains())
    ranked = {row["human_site"]: row for row in read_tsv(FINAL)}
    modified = [
        feature for feature in features if feature["type"] == "Modified residue"
    ]

    rows: list[dict[str, object]] = []
    for label in RETAINED:
        position = int(ranked[label]["human_DNMT3A1_aa"])
        overlaps = [
            feature for feature in features
            if feature["location"]["start"]["value"] <= position
            <= feature["location"]["end"]["value"]
        ]
        nearest = min(
            modified,
            key=lambda feature: abs(
                feature["location"]["start"]["value"] - position
            ),
        )
        nearest_position = nearest["location"]["start"]["value"]
        residue = chain[position]
        plddt = sum(atom.bfactor for atom in residue) / len(residue)
        exact_modified = [
            feature.get("description", "")
            for feature in overlaps if feature["type"] == "Modified residue"
        ]
        exact_variants = [
            feature.get("description", "")
            for feature in overlaps if feature["type"] == "Natural variant"
        ]
        rows.append({
            "priority_rank": ranked[label]["priority_rank"],
            "human_site": label,
            "human_DNMT3A1_aa": position,
            "human_residue": sequence[position - 1],
            "local_sequence_21aa": sequence[
                max(0, position - 11):min(len(sequence), position + 10)
            ],
            "uniprot_disordered_region": str(any(
                feature["type"] == "Region"
                and feature.get("description") == "Disordered"
                for feature in overlaps
            )).lower(),
            "dnmt3a1_specific_absent_from_isoform2": str(position <= 213).lower(),
            "alphafold_plddt": plddt,
            "alphafold_confidence_class": (
                "very_low" if plddt < 50 else
                "low" if plddt < 70 else
                "confident" if plddt < 90 else "very_high"
            ),
            "exact_uniprot_modified_residue": ";".join(exact_modified),
            "exact_uniprot_natural_variant": ";".join(exact_variants),
            "nearest_uniprot_modified_residue": (
                f"{nearest.get('description', '')}@{nearest_position}"
            ),
            "distance_to_nearest_modified_residue_aa": abs(
                nearest_position - position
            ),
            "overlapping_uniprot_features": ";".join(
                f"{feature['type']}:{feature.get('description', '')}"
                for feature in overlaps
            ),
            "experimental_priority": (
                "regulatory_tail_substitution_or_deletion_assay;"
                "avoid_structure_contact_claim"
            ),
        })
    rows.sort(key=lambda row: int(row["priority_rank"]))
    with (OUT / "retained_candidate_context.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "uniprot_accession": entry["primaryAccession"],
        "uniprot_entry": entry["uniProtkbId"],
        "alphafold_model": ALPHAFOLD.name,
        "retained_candidates": rows,
        "shared_context": (
            "G34, T12, and S97 are in the UniProt-annotated disordered "
            "DNMT3A1-specific N terminus and have very-low AlphaFold "
            "confidence. None has an exact UniProt modified-residue or "
            "natural-variant annotation."
        ),
        "s97_context": (
            "S97 is eight residues from annotated phosphoserine S105. "
            "This proximity motivates a regulatory-tail/PTM hypothesis but "
            "does not establish that S97 itself is phosphorylated."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
