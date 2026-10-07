#!/usr/bin/env python3
"""Combine primary and sensitivity FUBAR results in human DNMT3A coordinates."""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta  # noqa: E402


CROSSWALK = ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv"
PRIMARY = ROOT / "results/04_selection/mammals/DNMT3A.FUBAR.json"
SENSITIVITY = ROOT / "results/04_selection/mammals/DNMT3A.composition_pass.FUBAR.json"
MAFFT_SENSITIVITY = ROOT / "results/04_selection/mammals/DNMT3A.mafft_sensitivity.FUBAR.json"
MEME_CANDIDATE_TEMPLATE = ROOT / "results/04_selection/mammals/DNMT3A.site_{site}.MEME.json"
MAFFT_CROSSWALK = ROOT / "results/02_alignment/sensitivity/MAFFT_human_coordinate_crosswalk.tsv"
ALIGNMENT = ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
MAFFT_ALIGNMENT = (
    ROOT / "results/02_alignment/sensitivity/DNMT3A_MAFFT_clean_mammals_HyPhy_unique.fasta"
)
TAXONOMY = ROOT / "results/01_curated/taxon_metadata.tsv"
OUTPUT_DIR = ROOT / "results/04_selection/mammals"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fubar_rows(path: Path) -> list[dict[str, float]]:
    data = json.loads(path.read_text())
    headers = [entry[0] if isinstance(entry, list) else str(entry) for entry in data["MLE"]["headers"]]
    values = data["MLE"]["content"]["0"]
    rows: list[dict[str, float]] = []
    for site, row in enumerate(values, start=1):
        rows.append({"site": site, **{header: row[index] for index, header in enumerate(headers)}})
    return rows


def meme_row(path: Path, site: int) -> dict[str, float]:
    data = json.loads(path.read_text())
    headers = [entry[0] if isinstance(entry, list) else str(entry) for entry in data["MLE"]["headers"]]
    values = data["MLE"]["content"]["0"][site - 1]
    return {header: values[index] for index, header in enumerate(headers)}


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    crosswalk = {
        int(row["filtered_codon_site_1based"]): row for row in read_tsv(CROSSWALK)
    }
    primary = {int(row["site"]): row for row in fubar_rows(PRIMARY)}
    sensitivity = {int(row["site"]): row for row in fubar_rows(SENSITIVITY)}
    if set(primary) != set(crosswalk) or set(sensitivity) != set(crosswalk):
        raise SystemExit("ERROR: FUBAR and human-coordinate site sets differ")

    mafft_crosswalk = {
        int(row["filtered_codon_site_1based"]): row for row in read_tsv(MAFFT_CROSSWALK)
    }
    mafft = {int(row["site"]): row for row in fubar_rows(MAFFT_SENSITIVITY)}
    mafft_by_human = {
        int(mafft_crosswalk[site]["human_DNMT3A1_aa"]): row
        for site, row in mafft.items()
        if mafft_crosswalk[site]["human_DNMT3A1_aa"]
    }
    primary_site_by_human = {
        int(row["human_DNMT3A1_aa"]): site
        for site, row in crosswalk.items()
        if row["human_DNMT3A1_aa"]
    }
    threshold = 0.9
    candidate_positions = {
        int(crosswalk[site]["human_DNMT3A1_aa"])
        for site in crosswalk
        if crosswalk[site]["human_DNMT3A1_aa"] and (
            float(primary[site]["Prob[alpha<beta]"]) >= threshold
            or float(sensitivity[site]["Prob[alpha<beta]"]) >= threshold
        )
    } | {
        position for position, row in mafft_by_human.items()
        if float(row["Prob[alpha<beta]"]) >= threshold
    }
    candidate_sites = sorted(primary_site_by_human[position] for position in candidate_positions)
    records = read_fasta(ALIGNMENT)
    mafft_records = {record.identifier: record for record in read_fasta(MAFFT_ALIGNMENT)}
    taxon_order = {
        row["safe_id"]: row["order"] for row in read_tsv(TAXONOMY)
        if row["is_mammal"] == "true"
    }

    candidate_rows: list[dict[str, object]] = []
    state_rows: list[dict[str, object]] = []
    order_states: dict[tuple[int, str], Counter[str]] = defaultdict(Counter)
    for site in candidate_sites:
        p = float(primary[site]["Prob[alpha<beta]"])
        s = float(sensitivity[site]["Prob[alpha<beta]"])
        human_position = int(crosswalk[site]["human_DNMT3A1_aa"])
        mafft_row = mafft_by_human[human_position]
        m = float(mafft_row["Prob[alpha<beta]"])
        supported = [
            label for label, value in (
                ("primary", p), ("composition_pass", s), ("mafft_alignment", m)
            ) if value >= threshold
        ]
        classification = "supported_" + "_and_".join(supported) + "_only"
        if len(supported) == 3:
            classification = "supported_in_all_three"
        states: Counter[str] = Counter()
        mafft_site = next(
            alt_site for alt_site, row in mafft_crosswalk.items()
            if row["human_DNMT3A1_aa"] == str(human_position)
        )
        amino_acid_agreements = 0
        codon_agreements = 0
        for record in records:
            codon = record.sequence[(site - 1) * 3:site * 3]
            residue = "-" if codon == "---" else CODON_TABLE.get(codon, "X")
            states[residue] += 1
            order_states[(site, taxon_order[record.identifier])][residue] += 1
            mafft_codon = mafft_records[record.identifier].sequence[
                (mafft_site - 1) * 3:mafft_site * 3
            ]
            mafft_residue = "-" if mafft_codon == "---" else CODON_TABLE.get(mafft_codon, "X")
            amino_acid_agreements += residue == mafft_residue
            codon_agreements += codon == mafft_codon
        mapped = crosswalk[site]
        meme_path = Path(str(MEME_CANDIDATE_TEMPLATE).format(site=site))
        meme = meme_row(meme_path, site) if meme_path.exists() else {}
        candidate_rows.append({
            "filtered_codon_site_1based": site,
            "human_DNMT3A1_aa": mapped["human_DNMT3A1_aa"],
            "human_residue": mapped["human_residue"],
            "domain": mapped["domain"],
            "full_panel_occupancy": mapped["occupancy"],
            "mammal_unique_occupancy": f"{1 - states['-'] / len(records):.6f}",
            "macse_mafft_amino_acid_agreement": f"{amino_acid_agreements / len(records):.6f}",
            "macse_mafft_codon_agreement": f"{codon_agreements / len(records):.6f}",
            "primary_alpha": primary[site]["alpha"],
            "primary_beta": primary[site]["beta"],
            "primary_positive_posterior": p,
            "primary_bayes_factor": primary[site]["BayesFactor[alpha<beta]"],
            "composition_pass_positive_posterior": s,
            "composition_pass_bayes_factor": sensitivity[site]["BayesFactor[alpha<beta]"],
            "mafft_alignment_positive_posterior": m,
            "mafft_alignment_bayes_factor": mafft_row["BayesFactor[alpha<beta]"],
            "meme_lrt": meme.get("LRT", ""),
            "meme_p_value": meme.get("p-value", ""),
            "meme_branches_ebf_ge_100": meme.get("# branches under selection", ""),
            "sensitivity_classification": classification,
        })
        for residue, count in sorted(states.items(), key=lambda item: (-item[1], item[0])):
            state_rows.append({
                "filtered_codon_site_1based": site,
                "human_DNMT3A1_aa": mapped["human_DNMT3A1_aa"],
                "residue": residue,
                "count": count,
                "fraction": f"{count / len(records):.6f}",
            })

    order_rows: list[dict[str, object]] = []
    for (site, order), states in sorted(order_states.items()):
        order_rows.append({
            "filtered_codon_site_1based": site,
            "human_DNMT3A1_aa": crosswalk[site]["human_DNMT3A1_aa"],
            "order": order,
            "taxa": sum(states.values()),
            "residue_counts": ";".join(f"{aa}:{count}" for aa, count in states.most_common()),
        })

    write_tsv(
        OUTPUT_DIR / "fubar_candidate_sensitivity.tsv",
        candidate_rows,
        [
            "filtered_codon_site_1based", "human_DNMT3A1_aa", "human_residue",
            "domain", "full_panel_occupancy", "mammal_unique_occupancy",
            "macse_mafft_amino_acid_agreement", "macse_mafft_codon_agreement",
            "primary_alpha", "primary_beta", "primary_positive_posterior",
            "primary_bayes_factor", "composition_pass_positive_posterior",
            "composition_pass_bayes_factor", "mafft_alignment_positive_posterior",
            "mafft_alignment_bayes_factor", "meme_lrt", "meme_p_value",
            "meme_branches_ebf_ge_100", "sensitivity_classification",
        ],
    )
    write_tsv(
        OUTPUT_DIR / "fubar_candidate_residue_states.tsv",
        state_rows,
        ["filtered_codon_site_1based", "human_DNMT3A1_aa", "residue", "count", "fraction"],
    )
    write_tsv(
        OUTPUT_DIR / "fubar_candidate_order_states.tsv",
        order_rows,
        ["filtered_codon_site_1based", "human_DNMT3A1_aa", "order", "taxa", "residue_counts"],
    )
    print(json.dumps(candidate_rows, indent=2))


if __name__ == "__main__":
    main()
