#!/usr/bin/env python3
"""Summarize COBALT alignment and selection sensitivity in human coordinates."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta, write_tsv  # noqa: E402


ALIGNMENT = ROOT / "results/02_alignment/sensitivity/cobalt/DNMT3A_COBALT_clean_mammals_HyPhy_unique.fasta"
CROSSWALK = ROOT / "results/02_alignment/sensitivity/cobalt/COBALT_human_coordinate_crosswalk.tsv"
PRIMARY_ALIGNMENT = ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
PRIMARY_CROSSWALK = ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv"
OUT = ROOT / "results/04_selection/mammals/cobalt"
FUBAR = OUT / "DNMT3A.COBALT.FUBAR.json"
CANDIDATES = (5, 6, 12, 16, 18, 34, 97, 114)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def mapped(path: Path) -> dict[int, dict[str, str]]:
    return {
        int(row["human_DNMT3A1_aa"]): row
        for row in read_tsv(path)
        if row["human_DNMT3A1_aa"]
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


def amino_acid(codon: str) -> str:
    return "-" if codon == "---" else CODON_TABLE.get(codon, "X")


def holm(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [1.0] * len(values)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (len(values) - rank) * values[index]))
        result[index] = running
    return result


def main() -> None:
    primary_map, cobalt_map = mapped(PRIMARY_CROSSWALK), mapped(CROSSWALK)
    primary = {r.identifier: r.sequence for r in read_fasta(PRIMARY_ALIGNMENT)}
    cobalt = {r.identifier: r.sequence for r in read_fasta(ALIGNMENT)}
    if set(primary) != set(cobalt):
        raise RuntimeError("Primary and COBALT taxon sets differ")
    fit = fubar(FUBAR)
    if len(fit) != len(cobalt_map):
        raise RuntimeError("COBALT FUBAR and crosswalk site sets differ")

    site_rows: list[dict[str, object]] = []
    for position in sorted(cobalt_map):
        cobalt_site = int(cobalt_map[position]["filtered_codon_site_1based"])
        primary_site = (
            int(primary_map[position]["filtered_codon_site_1based"])
            if position in primary_map else None
        )
        aa_matches = codon_matches = 0
        differences: list[str] = []
        if primary_site is not None:
            for identifier in primary:
                first = primary[identifier][(primary_site - 1) * 3:primary_site * 3]
                second = cobalt[identifier][(cobalt_site - 1) * 3:cobalt_site * 3]
                aa_matches += amino_acid(first) == amino_acid(second)
                codon_matches += first == second
                if amino_acid(first) != amino_acid(second):
                    differences.append(
                        f"{identifier}:{amino_acid(first)}/{first}>"
                        f"{amino_acid(second)}/{second}"
                    )
        site_rows.append({
            "human_DNMT3A1_aa": position,
            "human_residue": cobalt_map[position]["human_residue"],
            "domain": cobalt_map[position]["domain"],
            "primary_filtered_site": primary_site,
            "cobalt_filtered_site": cobalt_site,
            "cobalt_occupancy": cobalt_map[position]["occupancy"],
            "primary_cobalt_amino_acid_agreement": (
                aa_matches / len(primary) if primary_site is not None else ""
            ),
            "primary_cobalt_codon_agreement": (
                codon_matches / len(primary) if primary_site is not None else ""
            ),
            "amino_acid_disagreement_count": len(differences),
            "amino_acid_disagreements": ";".join(differences),
            "cobalt_fubar_positive_posterior": fit[cobalt_site],
        })
    write_tsv(OUT / "cobalt_alignment_site_agreement.tsv", site_rows)

    by_position = {int(row["human_DNMT3A1_aa"]): row for row in site_rows}
    candidate_rows: list[dict[str, object]] = []
    for position in CANDIDATES:
        row = dict(by_position[position])
        label = f"{row['human_residue']}{position}"
        data = json.loads((OUT / f"DNMT3A.COBALT.{label}.MEME.json").read_text())
        site = int(row["cobalt_filtered_site"])
        values = data["MLE"]["content"]["0"][site - 1]
        row.update({
            "human_site": label,
            "cobalt_meme_lrt": float(values[5]),
            "cobalt_meme_p": float(values[6]),
            "cobalt_meme_holm_p_8_candidates": "",
            "cobalt_meme_branches_EBF_ge_100": int(values[7]),
        })
        candidate_rows.append(row)
    adjusted = holm([float(row["cobalt_meme_p"]) for row in candidate_rows])
    for row, value in zip(candidate_rows, adjusted):
        row["cobalt_meme_holm_p_8_candidates"] = value
    write_tsv(OUT / "cobalt_candidate_sensitivity.tsv", candidate_rows)

    fubar_candidates: list[dict[str, object]] = []
    for position, row in by_position.items():
        posterior = float(row["cobalt_fubar_positive_posterior"])
        if posterior >= 0.90:
            fubar_candidates.append({
                "human_site": f"{row['human_residue']}{position}",
                "human_DNMT3A1_aa": position,
                "cobalt_filtered_site": row["cobalt_filtered_site"],
                "cobalt_fubar_positive_posterior": posterior,
            })
    write_tsv(OUT / "cobalt_fubar_candidates.tsv", fubar_candidates)
    summary = {
        "cobalt_human_positions": len(site_rows),
        "shared_primary_cobalt_human_positions": sum(
            row["primary_filtered_site"] is not None for row in site_rows
        ),
        "mean_primary_cobalt_amino_acid_agreement": sum(
            float(row["primary_cobalt_amino_acid_agreement"])
            for row in site_rows
            if row["primary_cobalt_amino_acid_agreement"] != ""
        ) / sum(
            row["primary_cobalt_amino_acid_agreement"] != "" for row in site_rows
        ),
        "cobalt_fubar_candidates_ge_0.90": [
            row["human_site"] for row in fubar_candidates
        ],
        "exact_meme_candidates": [
            {
                "human_site": row["human_site"],
                "p": row["cobalt_meme_p"],
                "holm_p_8": row["cobalt_meme_holm_p_8_candidates"],
            }
            for row in candidate_rows
        ],
        "constraint_mode": "norps_local_similarity_only",
    }
    (OUT / "cobalt_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
