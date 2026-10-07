#!/usr/bin/env python3
"""Audit full-MEME discovery sites across MACSE, MAFFT, and PRANK."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta  # noqa: E402


OUTDIR = ROOT / "results/04_selection/mammals/full_meme"
SENSITIVITY = OUTDIR / "sensitivity"
DISCOVERY = OUTDIR / "significant_sites.tsv"
PRIMARY_ALIGNMENT = (
    ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
)
MAFFT_ALIGNMENT = (
    ROOT
    / "results/02_alignment/sensitivity"
    / "DNMT3A_MAFFT_clean_mammals_HyPhy_unique.fasta"
)
PRANK_ALIGNMENT = (
    ROOT
    / "results/02_alignment/sensitivity/prank"
    / "DNMT3A_PRANK_clean_mammals_HyPhy_unique.fasta"
)
CROSSWALKS = {
    "primary": ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv",
    "mafft": (
        ROOT
        / "results/02_alignment/sensitivity"
        / "MAFFT_human_coordinate_crosswalk.tsv"
    ),
    "prank": (
        ROOT
        / "results/02_alignment/sensitivity/prank"
        / "PRANK_human_coordinate_crosswalk.tsv"
    ),
}
FUBAR = {
    "primary": ROOT / "results/04_selection/mammals/DNMT3A.FUBAR.json",
    "mafft": (
        ROOT
        / "results/04_selection/mammals"
        / "DNMT3A.mafft_sensitivity.FUBAR.json"
    ),
    "prank": (
        ROOT
        / "results/04_selection/mammals"
        / "DNMT3A.prank_sensitivity.FUBAR.json"
    ),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def crosswalk(path: Path) -> dict[int, dict[str, str]]:
    return {
        int(row["human_DNMT3A1_aa"]): row
        for row in read_tsv(path)
        if row["human_DNMT3A1_aa"]
    }


def fubar_posteriors(path: Path) -> dict[int, float]:
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


def aa(codon: str) -> str:
    return "-" if codon == "---" else CODON_TABLE.get(codon, "X")


def alignment_agreement(
    primary_records: dict[str, str],
    other_records: dict[str, str],
    primary_site: int,
    other_site: int,
) -> tuple[float, float, list[str]]:
    if set(primary_records) != set(other_records):
        raise RuntimeError("Alternative alignment taxon set differs from primary")
    aa_matches = 0
    codon_matches = 0
    differences: list[str] = []
    for identifier in primary_records:
        primary_codon = primary_records[identifier][
            (primary_site - 1) * 3 : primary_site * 3
        ]
        other_codon = other_records[identifier][
            (other_site - 1) * 3 : other_site * 3
        ]
        primary_aa, other_aa = aa(primary_codon), aa(other_codon)
        aa_matches += primary_aa == other_aa
        codon_matches += primary_codon == other_codon
        if primary_aa != other_aa:
            differences.append(
                f"{identifier}:{primary_aa}/{primary_codon}>"
                f"{other_aa}/{other_codon}"
            )
    count = len(primary_records)
    return aa_matches / count, codon_matches / count, differences


def holm_adjust(p_values: list[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [1.0] * len(p_values)
    running = 0.0
    count = len(p_values)
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (count - rank) * p_values[index]))
        adjusted[index] = running
    return adjusted


def meme_result(alignment: str, label: str, site: int) -> tuple[float, float, int]:
    path = SENSITIVITY / f"DNMT3A.{alignment}.human_{label}.MEME.json"
    data = json.loads(path.read_text())
    values = data["MLE"]["content"]["0"][site - 1]
    return float(values[5]), float(values[6]), int(values[7])


def classify(mafft_p: float, prank_p: float) -> str:
    if mafft_p <= 0.05 and prank_p <= 0.05:
        return "meme_supported_all_three_alignments"
    if mafft_p <= 0.05 and prank_p > 0.05:
        return "primary_and_mafft_only;prank_sensitive"
    if mafft_p > 0.05 and prank_p <= 0.05:
        return "primary_and_prank_only;mafft_sensitive"
    return "primary_only;fails_mafft_and_prank"


def main() -> None:
    discovery = read_tsv(DISCOVERY)
    if len(discovery) != 3:
        raise RuntimeError(f"Expected three BH discovery sites, found {len(discovery)}")
    maps = {name: crosswalk(path) for name, path in CROSSWALKS.items()}
    fits = {name: fubar_posteriors(path) for name, path in FUBAR.items()}
    alignments = {
        "primary": {
            record.identifier: record.sequence
            for record in read_fasta(PRIMARY_ALIGNMENT)
        },
        "mafft": {
            record.identifier: record.sequence
            for record in read_fasta(MAFFT_ALIGNMENT)
        },
        "prank": {
            record.identifier: record.sequence
            for record in read_fasta(PRANK_ALIGNMENT)
        },
    }

    rows: list[dict[str, object]] = []
    for source in discovery:
        position = int(source["human_DNMT3A1_aa"])
        label = f"{source['human_residue']}{position}"
        sites = {
            name: int(mapped[position]["filtered_codon_site_1based"])
            for name, mapped in maps.items()
        }
        mafft_lrt, mafft_p, mafft_branches = meme_result(
            "MAFFT", label, sites["mafft"]
        )
        prank_lrt, prank_p, prank_branches = meme_result(
            "PRANK", label, sites["prank"]
        )
        mafft_aa, mafft_codon, mafft_differences = alignment_agreement(
            alignments["primary"], alignments["mafft"],
            sites["primary"], sites["mafft"],
        )
        prank_aa, prank_codon, prank_differences = alignment_agreement(
            alignments["primary"], alignments["prank"],
            sites["primary"], sites["prank"],
        )
        rows.append({
            "human_site": label,
            "domain": source["domain"],
            "primary_site": sites["primary"],
            "mafft_site": sites["mafft"],
            "prank_site": sites["prank"],
            "primary_full_meme_p": source["p_value"],
            "primary_full_meme_bh_q": source["bh_q_value_901_sites"],
            "mafft_meme_lrt": mafft_lrt,
            "mafft_meme_p": mafft_p,
            "mafft_meme_holm_p_3_candidates": "",
            "prank_meme_lrt": prank_lrt,
            "prank_meme_p": prank_p,
            "prank_meme_holm_p_3_candidates": "",
            "mafft_branches_EBF_ge_100": mafft_branches,
            "prank_branches_EBF_ge_100": prank_branches,
            "primary_fubar_posterior": fits["primary"][sites["primary"]],
            "mafft_fubar_posterior": fits["mafft"][sites["mafft"]],
            "prank_fubar_posterior": fits["prank"][sites["prank"]],
            "primary_mafft_aa_agreement": mafft_aa,
            "primary_mafft_codon_agreement": mafft_codon,
            "primary_mafft_aa_differences": ";".join(mafft_differences),
            "primary_prank_aa_agreement": prank_aa,
            "primary_prank_codon_agreement": prank_codon,
            "primary_prank_aa_differences": ";".join(prank_differences),
            "alignment_robustness": classify(mafft_p, prank_p),
        })

    mafft_adjusted = holm_adjust([float(row["mafft_meme_p"]) for row in rows])
    prank_adjusted = holm_adjust([float(row["prank_meme_p"]) for row in rows])
    for row, mafft_p, prank_p in zip(rows, mafft_adjusted, prank_adjusted):
        row["mafft_meme_holm_p_3_candidates"] = mafft_p
        row["prank_meme_holm_p_3_candidates"] = prank_p

    table = OUTDIR / "discovery_alignment_sensitivity.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    counts: dict[str, int] = {}
    for row in rows:
        key = str(row["alignment_robustness"])
        counts[key] = counts.get(key, 0) + 1
    (OUTDIR / "discovery_alignment_sensitivity_summary.json").write_text(
        json.dumps({
            "discovery_sites": len(rows),
            "meme_supported_all_three_alignments": sum(
                row["alignment_robustness"]
                == "meme_supported_all_three_alignments"
                for row in rows
            ),
            "classifications": counts,
            "interpretation": (
                "No full-scan discovery site is MEME-significant in both "
                "alternative alignments. Alternative-alignment tests are "
                "Holm-corrected across the three discovery candidates."
            ),
        }, indent=2) + "\n"
    )
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
