#!/usr/bin/env python3
"""Extract descriptive branch evidence for targeted DNMT3A MEME candidates."""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from io import StringIO
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta  # noqa: E402


RESULTS = ROOT / "results/04_selection/mammals"
CROSSWALK = ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv"
TAXONOMY = ROOT / "results/01_curated/taxon_metadata.tsv"
COMPOSITION_FAILURES = ROOT / "results/03_tree/mammals/composition_chi2_failures.tsv"
ALIGNMENT = ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
CANDIDATE_SITES = (9, 30, 107)
HIGH_POSTERIOR_THRESHOLD = 0.90


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def translate(codon: str) -> str:
    if codon == "---":
        return "-"
    return CODON_TABLE.get(codon.upper(), "X")


def change_class(ancestral: str, derived: str) -> str:
    if "---" in (ancestral, derived):
        return "gap_or_missing"
    ancestral_aa, derived_aa = translate(ancestral), translate(derived)
    if "X" in (ancestral_aa, derived_aa):
        return "ambiguous"
    return "synonymous" if ancestral_aa == derived_aa else "nonsynonymous"


def deepest_shared_lineage(descendants: list[str], taxonomy: dict[str, dict[str, str]]) -> str:
    paths = [taxonomy[tip]["lineage"].split(";") for tip in descendants]
    shared: list[str] = []
    for values in zip(*paths):
        if len(set(values)) != 1:
            break
        shared.append(values[0])
    return shared[-1] if shared else "Mammalia"


def main() -> None:
    crosswalk = {
        int(row["filtered_codon_site_1based"]): row for row in read_tsv(CROSSWALK)
    }
    taxonomy = {row["safe_id"]: row for row in read_tsv(TAXONOMY)}
    composition_failures = {
        row["safe_id"] for row in read_tsv(COMPOSITION_FAILURES)
    }
    observed_sequences = {
        record.identifier: record.sequence for record in read_fasta(ALIGNMENT)
    }

    branch_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []

    for site in CANDIDATE_SITES:
        path = RESULTS / f"DNMT3A.site_{site}.MEME.json"
        data = json.loads(path.read_text())
        mle = data["MLE"]["content"]["0"][site - 1]
        p_value = float(mle[6])
        site_significant = p_value <= 0.05
        mapped = crosswalk[site]
        tree = Phylo.read(StringIO(data["input"]["trees"]["0"]), "newick")
        attributes = data["branch attributes"]["0"]
        substitutions = data["substitutions"]["0"][str(site - 1)]
        if not isinstance(substitutions, dict) or "root" not in substitutions:
            raise RuntimeError(f"Missing ancestral reconstruction for site {site}")

        parents = {child: parent for parent in tree.find_clades() for child in parent.clades}
        states = {tree.root: substitutions["root"]}
        site_branches: list[dict[str, object]] = []
        site_events: list[dict[str, object]] = []

        for clade in tree.find_clades(order="preorder"):
            if clade is tree.root:
                continue
            branch = clade.name
            if not branch:
                raise RuntimeError(f"Unnamed non-root branch at site {site}")
            parent_state = states[parents[clade]]
            derived_state = substitutions.get(branch, parent_state)
            states[clade] = derived_state
            posterior = float(
                attributes[branch]["Posterior prob omega class by site"][1][site - 1]
            )
            descendants = [tip.name for tip in clade.get_terminals()]
            orders = sorted({taxonomy[tip]["order"] for tip in descendants})
            shared_lineage = deepest_shared_lineage(descendants, taxonomy)
            failed_descendants = sorted(set(descendants) & composition_failures)
            branch_type = "terminal" if clade.is_terminal() else "internal"
            label = (
                taxonomy[branch]["species"]
                if branch_type == "terminal"
                else (
                    f"{branch} MRCA ({len(descendants)} tips; "
                    f"{shared_lineage}; {', '.join(orders)})"
                )
            )
            reconstructed_change = branch in substitutions
            row = {
                "filtered_codon_site_1based": site,
                "human_DNMT3A1_aa": mapped["human_DNMT3A1_aa"],
                "human_residue": mapped["human_residue"],
                "site_meme_p_value": p_value,
                "site_meme_significant_p05": str(site_significant).lower(),
                "branch": branch,
                "branch_type": branch_type,
                "branch_label": label,
                "descendant_tip_count": len(descendants),
                "deepest_shared_lineage": shared_lineage,
                "descendant_orders": ";".join(orders),
                "positive_class_posterior": posterior,
                "posterior_ge_0_90": str(posterior >= HIGH_POSTERIOR_THRESHOLD).lower(),
                "reconstructed_change": str(reconstructed_change).lower(),
                "ancestral_codon": parent_state,
                "derived_codon": derived_state,
                "ancestral_aa": translate(parent_state),
                "derived_aa": translate(derived_state),
                "change_class": (
                    change_class(parent_state, derived_state)
                    if reconstructed_change else "none"
                ),
                "composition_failure_descendant_count": len(failed_descendants),
                "composition_failure_descendants": ";".join(failed_descendants),
                "descendant_tips": ";".join(descendants),
            }
            branch_rows.append(row)
            site_branches.append(row)
            if reconstructed_change:
                event_rows.append(row)
                site_events.append(row)

        reconstruction_mismatches = []
        for tip in tree.get_terminals():
            observed = observed_sequences[tip.name][(site - 1) * 3:site * 3]
            if states[tip] != observed:
                reconstruction_mismatches.append((tip.name, states[tip], observed))
        if reconstruction_mismatches:
            raise RuntimeError(
                f"Site {site} reconstruction disagrees with terminal alignment states: "
                f"{reconstruction_mismatches[:5]}"
            )

        posterior_ranks = sorted(
            site_branches, key=lambda row: float(row["positive_class_posterior"]), reverse=True
        )
        for rank, row in enumerate(posterior_ranks, start=1):
            row["posterior_rank"] = rank

        change_counts = Counter(str(row["change_class"]) for row in site_events)
        high = [
            row for row in posterior_ranks
            if float(row["positive_class_posterior"]) >= HIGH_POSTERIOR_THRESHOLD
        ]
        summaries.append({
            "filtered_codon_site_1based": site,
            "human_DNMT3A1_aa": int(mapped["human_DNMT3A1_aa"]),
            "human_residue": mapped["human_residue"],
            "meme_p_value": p_value,
            "meme_significant_p05": site_significant,
            "reconstructed_change_events": len(site_events),
            "nonsynonymous_events": change_counts["nonsynonymous"],
            "synonymous_events": change_counts["synonymous"],
            "gap_or_missing_events": change_counts["gap_or_missing"],
            "branches_posterior_ge_0_90": len(high),
            "terminal_reconstruction_mismatches": 0,
            "high_posterior_branches": [str(row["branch"]) for row in high],
        })

    branch_fields = [
        "filtered_codon_site_1based", "human_DNMT3A1_aa", "human_residue",
        "site_meme_p_value", "site_meme_significant_p05", "branch",
        "branch_type", "branch_label", "descendant_tip_count", "descendant_orders",
        "deepest_shared_lineage",
        "positive_class_posterior", "posterior_rank", "posterior_ge_0_90",
        "reconstructed_change", "ancestral_codon", "derived_codon",
        "ancestral_aa", "derived_aa", "change_class",
        "composition_failure_descendant_count", "composition_failure_descendants",
        "descendant_tips",
    ]
    write_tsv(RESULTS / "candidate_branch_evidence.tsv", branch_rows, branch_fields)
    write_tsv(RESULTS / "candidate_substitution_events.tsv", event_rows, branch_fields)
    write_tsv(
        RESULTS / "candidate_high_posterior_branches.tsv",
        [
            row for row in branch_rows
            if float(row["positive_class_posterior"]) >= HIGH_POSTERIOR_THRESHOLD
        ],
        branch_fields,
    )
    (RESULTS / "candidate_lineage_summary.json").write_text(
        json.dumps(summaries, indent=2) + "\n"
    )
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
