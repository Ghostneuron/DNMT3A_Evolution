#!/usr/bin/env python3
"""Extract T12/G34 lineage evidence from reliability-masked MEME fits."""

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
MEME_DIR = RESULTS / "reliability_mask"
OUT = MEME_DIR / "lineages"
ALIGNMENT = (
    ROOT
    / "results/02_alignment/reliability_mask"
    / "DNMT3A_mammals_HyPhy_unique_N_terminal_consensus_masked.fasta"
)
CROSSWALK = ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv"
TAXONOMY = ROOT / "results/01_curated/taxon_metadata.tsv"
COMPOSITION_FAILURES = (
    ROOT / "results/03_tree/mammals/composition_chi2_failures.tsv"
)
ORIGINAL_BRANCHES = RESULTS / "candidate_branch_evidence.tsv"
SITES = {"T12": 9, "G34": 30}
HIGH_POSTERIOR = 0.90


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def translate(codon: str) -> str:
    if codon == "---":
        return "-"
    return CODON_TABLE.get(codon.upper(), "X")


def informative(codon: str) -> bool:
    return len(codon) == 3 and all(base in "ACGT" for base in codon)


def change_class(ancestral: str, derived: str) -> str:
    if "---" in (ancestral, derived):
        return "gap_or_missing"
    ancestral_aa, derived_aa = translate(ancestral), translate(derived)
    if "X" in (ancestral_aa, derived_aa):
        return "ambiguous"
    return "synonymous" if ancestral_aa == derived_aa else "nonsynonymous"


def deepest_shared_lineage(
    descendants: list[str], taxonomy: dict[str, dict[str, str]]
) -> str:
    paths = [taxonomy[tip]["lineage"].split(";") for tip in descendants]
    shared: list[str] = []
    for values in zip(*paths):
        if len(set(values)) != 1:
            break
        shared.append(values[0])
    return shared[-1] if shared else "Mammalia"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    crosswalk = {
        int(row["filtered_codon_site_1based"]): row
        for row in read_tsv(CROSSWALK)
    }
    taxonomy = {row["safe_id"]: row for row in read_tsv(TAXONOMY)}
    composition_failures = {
        row["safe_id"] for row in read_tsv(COMPOSITION_FAILURES)
    }
    sequences = {
        record.identifier: record.sequence for record in read_fasta(ALIGNMENT)
    }
    original = {
        (int(row["filtered_codon_site_1based"]), row["branch"]): row
        for row in read_tsv(ORIGINAL_BRANCHES)
        if int(row["filtered_codon_site_1based"]) in set(SITES.values())
    }

    branch_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    comparison_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []

    for label, site in SITES.items():
        data = json.loads((MEME_DIR / f"DNMT3A.{label}.MEME.json").read_text())
        mle = data["MLE"]["content"]["0"][site - 1]
        meme_p = float(mle[6])
        tree = Phylo.read(StringIO(data["input"]["trees"]["0"]), "newick")
        attributes = data["branch attributes"]["0"]
        substitutions = data["substitutions"]["0"][str(site - 1)]
        if not isinstance(substitutions, dict) or "root" not in substitutions:
            raise RuntimeError(f"Missing ancestral reconstruction for {label}")

        parents = {
            child: parent for parent in tree.find_clades() for child in parent.clades
        }
        states = {tree.root: substitutions["root"]}
        site_rows: list[dict[str, object]] = []
        for clade in tree.find_clades(order="preorder"):
            if clade is tree.root:
                continue
            branch = clade.name
            if not branch:
                raise RuntimeError(f"Unnamed branch at {label}")
            parent_state = states[parents[clade]]
            derived_state = substitutions.get(branch, parent_state)
            states[clade] = derived_state
            posterior = float(
                attributes[branch]["Posterior prob omega class by site"][1][site - 1]
            )
            descendants = [tip.name for tip in clade.get_terminals()]
            orders = sorted({taxonomy[tip]["order"] for tip in descendants})
            shared = deepest_shared_lineage(descendants, taxonomy)
            failures = sorted(set(descendants) & composition_failures)
            branch_type = "terminal" if clade.is_terminal() else "internal"
            observed_terminal_codon = (
                sequences[branch][(site - 1) * 3:site * 3]
                if branch_type == "terminal" else ""
            )
            terminal_informative = (
                informative(observed_terminal_codon)
                if branch_type == "terminal" else True
            )
            reconstructed_change = branch in substitutions
            event_class = (
                change_class(parent_state, derived_state)
                if reconstructed_change else "none"
            )
            lineage_priority = (
                posterior >= HIGH_POSTERIOR
                and reconstructed_change
                and event_class == "nonsynonymous"
                and terminal_informative
            )
            branch_label = (
                taxonomy[branch]["species"]
                if branch_type == "terminal"
                else f"{branch} MRCA ({len(descendants)} tips; {shared})"
            )
            row = {
                "human_site": label,
                "filtered_codon_site_1based": site,
                "site_meme_p": meme_p,
                "branch": branch,
                "branch_type": branch_type,
                "branch_label": branch_label,
                "descendant_tip_count": len(descendants),
                "deepest_shared_lineage": shared,
                "descendant_orders": ";".join(orders),
                "positive_class_posterior": posterior,
                "posterior_ge_0_90": str(posterior >= HIGH_POSTERIOR).lower(),
                "observed_terminal_codon": observed_terminal_codon,
                "terminal_state_informative": (
                    str(terminal_informative).lower()
                    if branch_type == "terminal" else "not_applicable"
                ),
                "reconstructed_change": str(reconstructed_change).lower(),
                "ancestral_codon": parent_state,
                "derived_codon": derived_state,
                "ancestral_aa": translate(parent_state),
                "derived_aa": translate(derived_state),
                "change_class": event_class,
                "lineage_priority_change_coupled": str(lineage_priority).lower(),
                "composition_failure_descendant_count": len(failures),
                "composition_failure_descendants": ";".join(failures),
                "descendant_tips": ";".join(descendants),
            }
            branch_rows.append(row)
            site_rows.append(row)
            if reconstructed_change:
                event_rows.append(row)

        mismatches = []
        skipped_missing = 0
        for tip in tree.get_terminals():
            observed = sequences[tip.name][(site - 1) * 3:site * 3]
            if not informative(observed):
                skipped_missing += 1
                continue
            if states[tip] != observed:
                mismatches.append((tip.name, states[tip], observed))
        if mismatches:
            raise RuntimeError(f"{label} terminal reconstruction mismatches: {mismatches[:5]}")

        raw_masked_high = {
            str(row["branch"]): row
            for row in site_rows
            if float(row["positive_class_posterior"]) >= HIGH_POSTERIOR
        }
        masked_priority = {
            str(row["branch"]): row
            for row in site_rows
            if row["lineage_priority_change_coupled"] == "true"
        }
        original_priority = {
            branch: row
            for (original_site, branch), row in original.items()
            if original_site == site
            and float(row["positive_class_posterior"]) >= HIGH_POSTERIOR
            and row["reconstructed_change"] == "true"
            and row["change_class"] == "nonsynonymous"
        }
        for branch in sorted(set(masked_priority) | set(original_priority)):
            masked_row = masked_priority.get(branch)
            original_row = original_priority.get(branch)
            if masked_row and original_row:
                status = "retained"
            elif masked_row:
                status = "new_after_mask"
            else:
                status = "lost_after_mask"
            source = masked_row or original_row
            assert source is not None
            comparison_rows.append({
                "human_site": label,
                "branch": branch,
                "branch_type": source["branch_type"],
                "branch_label": source["branch_label"],
                "original_positive_class_posterior": (
                    original_row["positive_class_posterior"] if original_row else ""
                ),
                "masked_positive_class_posterior": (
                    masked_row["positive_class_posterior"] if masked_row else ""
                ),
                "status": status,
                "masked_reconstructed_change": (
                    masked_row["reconstructed_change"] if masked_row else ""
                ),
                "masked_change": (
                    f"{masked_row['ancestral_aa']}>{masked_row['derived_aa']}"
                    if masked_row and masked_row["reconstructed_change"] == "true"
                    else ""
                ),
            })

        changes = [
            row for row in site_rows if row["reconstructed_change"] == "true"
        ]
        change_counts = Counter(str(row["change_class"]) for row in changes)
        summaries.append({
            "human_site": label,
            "filtered_codon_site_1based": site,
            "masked_meme_p": meme_p,
            "informative_terminal_states": len(sequences) - skipped_missing,
            "masked_terminal_states_skipped": skipped_missing,
            "terminal_reconstruction_mismatches": 0,
            "reconstructed_change_events": len(changes),
            "nonsynonymous_events": change_counts["nonsynonymous"],
            "synonymous_events": change_counts["synonymous"],
            "raw_high_posterior_branch_count": len(raw_masked_high),
            "prioritized_nonsynonymous_change_branches": sorted(masked_priority),
            "original_prioritized_nonsynonymous_change_branches": sorted(
                original_priority
            ),
            "retained_prioritized_branches": sorted(
                set(masked_priority) & set(original_priority)
            ),
            "new_prioritized_branches": sorted(
                set(masked_priority) - set(original_priority)
            ),
            "lost_prioritized_branches": sorted(
                set(original_priority) - set(masked_priority)
            ),
        })

    write_tsv(OUT / "branch_evidence.tsv", branch_rows)
    write_tsv(OUT / "substitution_events.tsv", event_rows)
    write_tsv(
        OUT / "high_posterior_branches.tsv",
        [
            row for row in branch_rows
            if float(row["positive_class_posterior"]) >= HIGH_POSTERIOR
        ],
    )
    write_tsv(
        OUT / "prioritized_change_coupled_branches.tsv",
        [
            row for row in branch_rows
            if row["lineage_priority_change_coupled"] == "true"
        ],
    )
    write_tsv(OUT / "high_posterior_comparison.tsv", comparison_rows)
    summary = {
        "sites": summaries,
        "interpretation": (
            "Raw branch posteriors are descriptive statistics, not branch-wise "
            "significance tests. Because the T12 positive-class posterior saturates "
            "most branches after masking, lineage priority additionally requires a "
            "reconstructed nonsynonymous change and an informative terminal state "
            "(or an internal branch). Masked terminal states are excluded."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
