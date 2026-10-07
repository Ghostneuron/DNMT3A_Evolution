#!/usr/bin/env python3
"""Compare PRANK candidate coordinates and selection results with primary analyses."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta  # noqa: E402


OUT = ROOT / "results/04_selection/mammals"
PRIMARY_ALIGNMENT = ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
PRIMARY_CROSSWALK = ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv"
PRANK_DIR = ROOT / "results/02_alignment/sensitivity/prank"
PRANK_ALIGNMENT = PRANK_DIR / "DNMT3A_PRANK_clean_mammals_HyPhy_unique.fasta"
PRANK_CROSSWALK = PRANK_DIR / "PRANK_human_coordinate_crosswalk.tsv"
PRIMARY_FUBAR = OUT / "DNMT3A.FUBAR.json"
COMPOSITION_FUBAR = OUT / "DNMT3A.composition_pass.FUBAR.json"
MAFFT_FUBAR = OUT / "DNMT3A.mafft_sensitivity.FUBAR.json"
MAFFT_CROSSWALK = ROOT / "results/02_alignment/sensitivity/MAFFT_human_coordinate_crosswalk.tsv"
PRANK_FUBAR = OUT / "DNMT3A.prank_sensitivity.FUBAR.json"
CANDIDATES = {
    12: OUT / "DNMT3A.PRANK.site_9_human_T12.MEME.json",
    34: OUT / "DNMT3A.PRANK.site_31_human_G34.MEME.json",
    114: OUT / "DNMT3A.PRANK.site_111_human_A114.MEME.json",
}
PRANK_MEME_FILES = {
    **CANDIDATES,
    16: OUT / "DNMT3A.PRANK.site_13_human_A16.MEME.json",
    97: OUT / "DNMT3A.PRANK.site_94_human_S97.MEME.json",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def crosswalk(path: Path) -> dict[int, dict[str, str]]:
    return {
        int(row["human_DNMT3A1_aa"]): row
        for row in read_tsv(path) if row["human_DNMT3A1_aa"]
    }


def fubar(path: Path) -> dict[int, dict[str, float]]:
    data = json.loads(path.read_text())
    headers = [item[0] if isinstance(item, list) else str(item) for item in data["MLE"]["headers"]]
    return {
        site: {header: values[index] for index, header in enumerate(headers)}
        for site, values in enumerate(data["MLE"]["content"]["0"], start=1)
    }


def aa(codon: str) -> str:
    return "-" if codon == "---" else CODON_TABLE.get(codon, "X")


def main() -> None:
    primary_map = crosswalk(PRIMARY_CROSSWALK)
    prank_map = crosswalk(PRANK_CROSSWALK)
    primary_records = {
        record.identifier: record.sequence for record in read_fasta(PRIMARY_ALIGNMENT)
    }
    prank_records = {
        record.identifier: record.sequence for record in read_fasta(PRANK_ALIGNMENT)
    }
    if set(primary_records) != set(prank_records):
        raise RuntimeError("Primary and PRANK taxon sets differ")
    primary_fit = fubar(PRIMARY_FUBAR)
    composition_fit = fubar(COMPOSITION_FUBAR)
    mafft_fit = fubar(MAFFT_FUBAR)
    mafft_map = crosswalk(MAFFT_CROSSWALK)
    prank_fit = fubar(PRANK_FUBAR)
    if len(prank_fit) != len(prank_map):
        raise RuntimeError("PRANK FUBAR and crosswalk site sets differ")

    site_rows: list[dict[str, object]] = []
    for position in sorted(set(primary_map) & set(prank_map)):
        primary_site = int(primary_map[position]["filtered_codon_site_1based"])
        prank_site = int(prank_map[position]["filtered_codon_site_1based"])
        aa_agree = 0
        codon_agree = 0
        disagreements: list[str] = []
        for identifier in primary_records:
            primary_codon = primary_records[identifier][
                (primary_site - 1) * 3:primary_site * 3
            ]
            prank_codon = prank_records[identifier][
                (prank_site - 1) * 3:prank_site * 3
            ]
            primary_aa, prank_aa = aa(primary_codon), aa(prank_codon)
            aa_agree += primary_aa == prank_aa
            codon_agree += primary_codon == prank_codon
            if primary_aa != prank_aa:
                disagreements.append(
                    f"{identifier}:{primary_aa}/{primary_codon}>{prank_aa}/{prank_codon}"
                )
        site_rows.append({
            "human_DNMT3A1_aa": position,
            "human_residue": primary_map[position]["human_residue"],
            "domain": primary_map[position]["domain"],
            "primary_filtered_site": primary_site,
            "prank_filtered_site": prank_site,
            "primary_occupancy": primary_map[position]["occupancy"],
            "prank_occupancy": prank_map[position]["occupancy"],
            "amino_acid_agreement": aa_agree / len(primary_records),
            "codon_agreement": codon_agree / len(primary_records),
            "amino_acid_disagreement_count": len(disagreements),
            "amino_acid_disagreements": ";".join(disagreements),
            "primary_fubar_positive_posterior": primary_fit[primary_site]["Prob[alpha<beta]"],
            "prank_fubar_positive_posterior": prank_fit[prank_site]["Prob[alpha<beta]"],
        })

    by_position = {int(row["human_DNMT3A1_aa"]): row for row in site_rows}
    candidate_rows: list[dict[str, object]] = []
    for position, meme_path in CANDIDATES.items():
        row = dict(by_position[position])
        prank_site = int(row["prank_filtered_site"])
        data = json.loads(meme_path.read_text())
        values = data["MLE"]["content"]["0"][prank_site - 1]
        row.update({
            "prank_meme_lrt": values[5],
            "prank_meme_p_value": values[6],
            "prank_meme_approximate_branches": values[7],
        })
        candidate_rows.append(row)

    prank_candidate_rows: list[dict[str, object]] = []
    for position, mapped in prank_map.items():
        site = int(mapped["filtered_codon_site_1based"])
        posterior = float(prank_fit[site]["Prob[alpha<beta]"])
        if posterior >= 0.90:
            primary_posterior = (
                float(primary_fit[int(primary_map[position]["filtered_codon_site_1based"])]["Prob[alpha<beta]"])
                if position in primary_map else None
            )
            composition_posterior = (
                float(composition_fit[int(primary_map[position]["filtered_codon_site_1based"])]["Prob[alpha<beta]"])
                if position in primary_map else None
            )
            mafft_posterior = (
                float(mafft_fit[int(mafft_map[position]["filtered_codon_site_1based"])]["Prob[alpha<beta]"])
                if position in mafft_map else None
            )
            meme_data = json.loads(PRANK_MEME_FILES[position].read_text())
            meme_values = meme_data["MLE"]["content"]["0"][site - 1]
            support = (
                "primary_and_prank_fubar"
                if primary_posterior is not None and primary_posterior >= 0.90
                else "prank_fubar_only"
            )
            support += (
                "_and_prank_meme_p05"
                if float(meme_values[6]) <= 0.05 else "_without_prank_meme_p05"
            )
            prank_candidate_rows.append({
                "prank_filtered_site": site,
                "human_DNMT3A1_aa": position,
                "human_residue": mapped["human_residue"],
                "prank_occupancy": mapped["occupancy"],
                "prank_fubar_positive_posterior": posterior,
                "primary_fubar_positive_posterior": primary_posterior if primary_posterior is not None else "",
                "composition_fubar_positive_posterior": composition_posterior if composition_posterior is not None else "",
                "mafft_fubar_positive_posterior": mafft_posterior if mafft_posterior is not None else "",
                "prank_meme_lrt": meme_values[5],
                "prank_meme_p_value": meme_values[6],
                "support_classification": support,
            })

    fields = list(site_rows[0])
    write_tsv(OUT / "prank_alignment_site_agreement.tsv", site_rows, fields)
    write_tsv(
        OUT / "prank_candidate_sensitivity.tsv",
        candidate_rows,
        list(candidate_rows[0]),
    )
    write_tsv(
        OUT / "prank_fubar_candidates.tsv",
        prank_candidate_rows,
        list(prank_candidate_rows[0]) if prank_candidate_rows else [
            "prank_filtered_site", "human_DNMT3A1_aa", "human_residue",
            "prank_occupancy", "prank_fubar_positive_posterior",
            "primary_fubar_positive_posterior",
            "composition_fubar_positive_posterior",
            "mafft_fubar_positive_posterior", "prank_meme_lrt",
            "prank_meme_p_value", "support_classification",
        ],
    )
    summary = {
        "shared_human_positions": len(site_rows),
        "mean_amino_acid_agreement": sum(
            float(row["amino_acid_agreement"]) for row in site_rows
        ) / len(site_rows),
        "prank_fubar_candidates_ge_0.90": len(prank_candidate_rows),
        "original_candidate_positions": sorted(CANDIDATES),
        "prank_fubar_candidate_positions": sorted(
            int(row["human_DNMT3A1_aa"]) for row in prank_candidate_rows
        ),
    }
    (OUT / "prank_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
