#!/usr/bin/env python3
"""Extract descriptive S97 substitution lineages from mammal-scope MEME."""

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
from dnmt3a_pipeline import CODON_TABLE, read_fasta, write_tsv  # noqa: E402


SITE = 94
HUMAN_SITE = "S97"
HIGH_POSTERIOR = 0.90
RESULTS = ROOT / "results/04_selection/mammals"
MEME = RESULTS / "mammal_scope/DNMT3A.MACSE.mammal_scope.S97.MEME.json"
ALIGNMENT = ROOT / "results/02_alignment/sensitivity/mammal_scope/macse/DNMT3A_MACSE_mammal_scope_HyPhy_unique.fasta"
TAXONOMY = ROOT / "results/01_curated/taxon_metadata.tsv"
OUT = RESULTS / "s97_homology_audit/lineages"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def aa(codon: str) -> str:
    return "-" if codon == "---" else CODON_TABLE.get(codon, "X")


def shared_lineage(descendants: list[str], taxonomy: dict[str, dict[str, str]]) -> str:
    paths = [taxonomy[tip]["lineage"].split(";") for tip in descendants]
    shared: list[str] = []
    for values in zip(*paths):
        if len(set(values)) != 1:
            break
        shared.append(values[0])
    return shared[-1] if shared else "Mammalia"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = json.loads(MEME.read_text())
    tree = Phylo.read(StringIO(data["input"]["trees"]["0"]), "newick")
    attributes = data["branch attributes"]["0"]
    substitutions = data["substitutions"]["0"][str(SITE - 1)]
    taxonomy = {row["safe_id"]: row for row in read_tsv(TAXONOMY)}
    sequences = {
        record.identifier: record.sequence for record in read_fasta(ALIGNMENT)
    }
    parents = {
        child: parent for parent in tree.find_clades() for child in parent.clades
    }
    states = {tree.root: substitutions["root"]}
    rows: list[dict[str, object]] = []
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        branch = clade.name
        if not branch:
            raise RuntimeError("Unnamed S97 branch")
        ancestral = states[parents[clade]]
        derived = substitutions.get(branch, ancestral)
        states[clade] = derived
        changed = branch in substitutions
        change_class = (
            "none" if not changed
            else "synonymous" if aa(ancestral) == aa(derived)
            else "nonsynonymous"
        )
        posterior = float(
            attributes[branch]["Posterior prob omega class by site"][1][SITE - 1]
        )
        descendants = [tip.name for tip in clade.get_terminals()]
        branch_type = "terminal" if clade.is_terminal() else "internal"
        lineage = shared_lineage(descendants, taxonomy)
        rows.append({
            "human_site": HUMAN_SITE,
            "alignment_site": SITE,
            "branch": branch,
            "branch_type": branch_type,
            "branch_label": (
                taxonomy[branch]["species"] if branch_type == "terminal"
                else f"{branch} MRCA ({len(descendants)} tips; {lineage})"
            ),
            "descendant_tip_count": len(descendants),
            "deepest_shared_lineage": lineage,
            "positive_class_posterior": posterior,
            "reconstructed_change": str(changed).lower(),
            "ancestral_codon": ancestral,
            "derived_codon": derived,
            "ancestral_aa": aa(ancestral),
            "derived_aa": aa(derived),
            "change_class": change_class,
            "priority_change_coupled": str(
                posterior >= HIGH_POSTERIOR and change_class == "nonsynonymous"
            ).lower(),
            "descendant_tips": ";".join(descendants),
        })

    mismatches = []
    for tip in tree.get_terminals():
        observed = sequences[tip.name][(SITE - 1) * 3:SITE * 3]
        if observed != states[tip]:
            mismatches.append((tip.name, states[tip], observed))
    if mismatches:
        raise RuntimeError(f"S97 reconstruction mismatches: {mismatches[:5]}")

    events = [row for row in rows if row["reconstructed_change"] == "true"]
    priority = [row for row in rows if row["priority_change_coupled"] == "true"]
    write_tsv(OUT / "branch_evidence.tsv", rows)
    write_tsv(OUT / "substitution_events.tsv", events)
    write_tsv(OUT / "prioritized_change_coupled_branches.tsv", priority)
    counts = Counter(str(row["change_class"]) for row in events)
    summary = {
        "human_site": HUMAN_SITE,
        "meme_p": float(data["MLE"]["content"]["0"][SITE - 1][6]),
        "root_codon": substitutions["root"],
        "root_residue": aa(substitutions["root"]),
        "terminal_reconstruction_mismatches": 0,
        "reconstructed_events": len(events),
        "nonsynonymous_events": counts["nonsynonymous"],
        "synonymous_events": counts["synonymous"],
        "prioritized_change_coupled_branches": [
            {
                "branch": row["branch"],
                "branch_label": row["branch_label"],
                "change": f"{row['ancestral_aa']}>{row['derived_aa']}",
                "posterior": row["positive_class_posterior"],
            }
            for row in priority
        ],
        "interpretation": (
            "Branches are descriptive MEME posterior/change-coupled leads, "
            "not branch-wide significance tests."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
