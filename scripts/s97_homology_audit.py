#!/usr/bin/env python3
"""Audit S97 occupancy, residue homology, and selection across four alignments."""

from __future__ import annotations

import csv
import itertools
import json
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta, write_tsv  # noqa: E402


OUT = ROOT / "results/04_selection/mammals/s97_homology_audit"
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
MEME = {
    "MACSE": ROOT / "results/04_selection/mammals/mammal_scope/DNMT3A.MACSE.mammal_scope.S97.MEME.json",
    "MAFFT": ROOT / "results/04_selection/mammals/mammal_scope/DNMT3A.MAFFT.mammal_scope.S97.MEME.json",
    "PRANK": ROOT / "results/04_selection/mammals/DNMT3A.PRANK.site_94_human_S97.MEME.json",
    "COBALT": ROOT / "results/04_selection/mammals/cobalt/DNMT3A.COBALT.S97.MEME.json",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def site(path: Path, position: int = 97) -> tuple[int, float]:
    row = next(
        row for row in read_tsv(path)
        if int(row["human_DNMT3A1_aa"]) == position
    )
    occupancy = row.get("mammal_occupancy", row.get("occupancy", ""))
    return int(row["filtered_codon_site_1based"]), float(occupancy)


def posterior(path: Path, alignment_site: int) -> float:
    data = json.loads(path.read_text())
    headers = [
        item[0] if isinstance(item, list) else str(item)
        for item in data["MLE"]["headers"]
    ]
    return float(
        data["MLE"]["content"]["0"][alignment_site - 1][
            headers.index("Prob[alpha<beta]")
        ]
    )


def meme(path: Path, alignment_site: int) -> tuple[float, float, int]:
    values = json.loads(path.read_text())["MLE"]["content"]["0"][alignment_site - 1]
    return float(values[5]), float(values[6]), int(values[7])


def amino_acid(codon: str) -> str:
    return "-" if codon == "---" else CODON_TABLE.get(codon, "X")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    maps = {name: site(crosswalk) for name, (_, crosswalk) in ALIGNMENTS.items()}
    records = {
        name: {record.identifier: record.sequence for record in read_fasta(alignment)}
        for name, (alignment, _) in ALIGNMENTS.items()
    }
    tip_sets = {frozenset(value) for value in records.values()}
    if len(tip_sets) != 1 or len(next(iter(tip_sets))) != 250:
        raise RuntimeError("Four-alignment taxon sets differ")

    states: dict[str, dict[str, str]] = {}
    codons: dict[str, dict[str, str]] = {}
    for name, rows in records.items():
        alignment_site = maps[name][0]
        codons[name] = {
            identifier: sequence[(alignment_site - 1) * 3:alignment_site * 3]
            for identifier, sequence in rows.items()
        }
        states[name] = {
            identifier: amino_acid(value)
            for identifier, value in codons[name].items()
        }

    agreement_rows: list[dict[str, object]] = []
    for first, second in itertools.combinations(ALIGNMENTS, 2):
        matches = sum(
            states[first][tip] == states[second][tip]
            for tip in states[first]
        )
        agreement_rows.append({
            "alignment_1": first,
            "alignment_2": second,
            "matching_taxa": matches,
            "taxa": len(states[first]),
            "amino_acid_agreement": matches / len(states[first]),
        })
    write_tsv(OUT / "pairwise_s97_alignment_agreement.tsv", agreement_rows)

    taxon_rows: list[dict[str, object]] = []
    for tip in sorted(states["MACSE"]):
        row: dict[str, object] = {"safe_id": tip}
        for name in ALIGNMENTS:
            row[f"{name.lower()}_codon"] = codons[name][tip]
            row[f"{name.lower()}_residue"] = states[name][tip]
        row["four_alignment_consensus"] = (
            len({states[name][tip] for name in ALIGNMENTS}) == 1
        )
        row["nonhuman_reference_state"] = states["MACSE"][tip] != "S"
        taxon_rows.append(row)
    write_tsv(OUT / "s97_taxon_states.tsv", taxon_rows)

    selection_rows: list[dict[str, object]] = []
    for name in ALIGNMENTS:
        alignment_site, occupancy = maps[name]
        lrt, p_value, branches = meme(MEME[name], alignment_site)
        selection_rows.append({
            "alignment": name,
            "alignment_site": alignment_site,
            "mammal_occupancy": occupancy,
            "fubar_positive_posterior": posterior(FUBAR[name], alignment_site),
            "meme_lrt": lrt,
            "meme_p": p_value,
            "meme_branches_EBF_ge_100": branches,
        })
    write_tsv(OUT / "s97_selection_sensitivity.tsv", selection_rows)

    primary_states = Counter(states["MACSE"].values())
    summary = {
        "human_site": "S97",
        "mammal_taxa": len(states["MACSE"]),
        "mammal_occupancy_all_alignments": {
            name: maps[name][1] for name in ALIGNMENTS
        },
        "all_pairwise_amino_acid_agreement": all(
            row["amino_acid_agreement"] == 1.0 for row in agreement_rows
        ),
        "residue_counts": dict(sorted(primary_states.items())),
        "non_serine_taxa": [
            row["safe_id"] for row in taxon_rows
            if row["nonhuman_reference_state"]
        ],
        "selection": selection_rows,
        "filtering_diagnosis": (
            "S97 was omitted from the original MACSE and MAFFT mammal "
            "subsets because occupancy filtering was performed across 475 "
            "vertebrates before mammal subsetting, not because mammalian "
            "residue homology was uncertain."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
