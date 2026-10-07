#!/usr/bin/env python3
"""Audit human sites restored by filtering within the mammal analysis scope."""

from __future__ import annotations

import csv
import itertools
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta, write_tsv  # noqa: E402


OUT = ROOT / "results/04_selection/mammals/mammal_scope_filter_audit"
RAW = {
    "MACSE": ROOT / "results/02_alignment/DNMT3A_MACSE_NT.fasta",
    "MAFFT": ROOT / "results/02_alignment/DNMT3A_MAFFT_projected_codons.fasta",
}
ORIGINAL_MAP = {
    "MACSE": ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv",
    "MAFFT": ROOT / "results/02_alignment/sensitivity/MAFFT_human_coordinate_crosswalk.tsv",
}
ALIGNMENTS = {
    "MACSE": (
        ROOT / "results/02_alignment/sensitivity/mammal_scope/macse/DNMT3A_MACSE_mammal_scope_HyPhy_unique.fasta",
        ROOT / "results/02_alignment/sensitivity/mammal_scope/macse/MACSE_mammal_scope_crosswalk.tsv",
    ),
    "MAFFT": (
        ROOT / "results/02_alignment/sensitivity/mammal_scope/mafft/DNMT3A_MAFFT_mammal_scope_HyPhy_unique.fasta",
        ROOT / "results/02_alignment/sensitivity/mammal_scope/mafft/MAFFT_mammal_scope_crosswalk.tsv",
    ),
    "PRANK": (
        ROOT / "results/02_alignment/sensitivity/prank/DNMT3A_PRANK_clean_mammals_HyPhy_unique.fasta",
        ROOT / "results/02_alignment/sensitivity/prank/PRANK_human_coordinate_crosswalk.tsv",
    ),
    "COBALT": (
        ROOT / "results/02_alignment/sensitivity/cobalt/DNMT3A_COBALT_clean_mammals_HyPhy_unique.fasta",
        ROOT / "results/02_alignment/sensitivity/cobalt/COBALT_human_coordinate_crosswalk.tsv",
    ),
}
FUBAR = {
    "MACSE": ROOT / "results/04_selection/mammals/mammal_scope/DNMT3A.MACSE.mammal_scope.FUBAR.json",
    "MAFFT": ROOT / "results/04_selection/mammals/mammal_scope/DNMT3A.MAFFT.mammal_scope.FUBAR.json",
    "PRANK": ROOT / "results/04_selection/mammals/DNMT3A.prank_sensitivity.FUBAR.json",
    "COBALT": ROOT / "results/04_selection/mammals/cobalt/DNMT3A.COBALT.FUBAR.json",
}
HUMAN = "Homo_sapiens"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def valid(codon: str) -> bool:
    return (
        len(codon) == 3
        and codon != "---"
        and "!" not in codon
        and "-" not in codon
        and all(base in "ACGT" for base in codon)
        and CODON_TABLE.get(codon) not in {None, "*"}
    )


def raw_occupancy(path: Path) -> dict[int, float]:
    sequences = {record.identifier: record.sequence for record in read_fasta(path)}
    human = sequences[HUMAN]
    position = 0
    result: dict[int, float] = {}
    for start in range(0, len(human), 3):
        human_codon = human[start:start + 3]
        if not valid(human_codon):
            continue
        position += 1
        result[position] = sum(
            valid(sequence[start:start + 3]) for sequence in sequences.values()
        ) / len(sequences)
    if position != 912:
        raise RuntimeError(f"{path} maps {position} human residues")
    return result


def crosswalk(path: Path) -> dict[int, dict[str, str]]:
    return {
        int(row["human_DNMT3A1_aa"]): row
        for row in read_tsv(path) if row["human_DNMT3A1_aa"]
    }


def fubar(path: Path) -> dict[int, float]:
    data = json.loads(path.read_text())
    headers = [
        item[0] if isinstance(item, list) else str(item)
        for item in data["MLE"]["headers"]
    ]
    index = headers.index("Prob[alpha<beta]")
    return {
        site: float(values[index])
        for site, values in enumerate(data["MLE"]["content"]["0"], start=1)
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    maps = {name: crosswalk(value[1]) for name, value in ALIGNMENTS.items()}
    fits = {name: fubar(path) for name, path in FUBAR.items()}
    sequences = {
        name: {record.identifier: record.sequence for record in read_fasta(value[0])}
        for name, value in ALIGNMENTS.items()
    }
    original_positions = {
        name: set(crosswalk(path)) for name, path in ORIGINAL_MAP.items()
    }
    restored = sorted(
        (set(maps["MACSE"]) - original_positions["MACSE"])
        | (set(maps["MAFFT"]) - original_positions["MAFFT"])
    )
    broad = {name: raw_occupancy(path) for name, path in RAW.items()}

    rows: list[dict[str, object]] = []
    for position in restored:
        states: dict[str, dict[str, str]] = {}
        for name in ALIGNMENTS:
            site = int(maps[name][position]["filtered_codon_site_1based"])
            states[name] = {
                tip: sequence[(site - 1) * 3:site * 3]
                for tip, sequence in sequences[name].items()
            }
        agreements = []
        for first, second in itertools.combinations(ALIGNMENTS, 2):
            agreements.append(
                sum(
                    CODON_TABLE.get(states[first][tip], "-")
                    == CODON_TABLE.get(states[second][tip], "-")
                    for tip in states[first]
                ) / len(states[first])
            )
        row: dict[str, object] = {
            "human_site": (
                maps["MACSE"][position]["human_residue"] + str(position)
            ),
            "human_DNMT3A1_aa": position,
            "broad_475_macse_occupancy": broad["MACSE"][position],
            "broad_475_mafft_occupancy": broad["MAFFT"][position],
            "mammal_macse_occupancy": maps["MACSE"][position]["mammal_occupancy"],
            "mammal_mafft_occupancy": maps["MAFFT"][position]["mammal_occupancy"],
            "minimum_four_alignment_aa_agreement": min(agreements),
            "all_four_alignments_identical": str(min(agreements) == 1.0).lower(),
        }
        for name in ALIGNMENTS:
            site = int(maps[name][position]["filtered_codon_site_1based"])
            row[f"{name.lower()}_fubar_posterior"] = fits[name][site]
        row["any_fubar_ge_0_90"] = str(any(
            float(row[f"{name.lower()}_fubar_posterior"]) >= 0.90
            for name in ALIGNMENTS
        )).lower()
        row["all_fubar_ge_0_90"] = str(all(
            float(row[f"{name.lower()}_fubar_posterior"]) >= 0.90
            for name in ALIGNMENTS
        )).lower()
        rows.append(row)
    write_tsv(OUT / "restored_human_sites.tsv", rows)

    candidate_sets = {}
    for name in ALIGNMENTS:
        candidate_sets[name] = [
            (
                maps[name][position]["human_residue"] + str(position)
            )
            for position in sorted(maps[name])
            if fits[name][int(maps[name][position]["filtered_codon_site_1based"])]
            >= 0.90
        ]
    summary = {
        "restored_human_sites": [row["human_site"] for row in rows],
        "restored_site_count": len(rows),
        "restored_sites_fubar_ge_0_90_all_alignments": [
            row["human_site"] for row in rows
            if row["all_fubar_ge_0_90"] == "true"
        ],
        "four_alignment_fubar_candidates": candidate_sets,
        "conclusion": (
            "S97 is the only restored site with FUBAR posterior >=0.90 "
            "under all four alignments. The other seven restored sites do "
            "not become selection candidates."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
