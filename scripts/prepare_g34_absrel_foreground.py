#!/usr/bin/env python3
"""Prepare an audited aBSREL foreground tree for robust internal G34 branches."""

from __future__ import annotations

import csv
import json
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Seq import Seq


ROOT = Path(__file__).resolve().parents[1]
SELECTION = ROOT / "results/04_selection/mammals"
MEME_JSON = SELECTION / "DNMT3A.site_30.MEME.json"
CANDIDATES = SELECTION / "candidate_high_posterior_branches.tsv"
IQTREE = ROOT / "results/03_tree/mammals/DNMT3A_mammals_tree.treefile"
ALIGNMENT = ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
PRANK_ALIGNMENT = (
    ROOT
    / "results/02_alignment/sensitivity/prank"
    / "DNMT3A_PRANK_clean_mammals_HyPhy_unique.fasta"
)
PRANK_CROSSWALK = (
    ROOT
    / "results/02_alignment/sensitivity/prank"
    / "PRANK_human_coordinate_crosswalk.tsv"
)
OUTDIR = SELECTION / "g34_absrel"
MANIFEST = OUTDIR / "foreground_manifest.tsv"
TREE_OUT = OUTDIR / "DNMT3A.G34_foreground.nwk"
SITE = 30
HUMAN_POSITION = 34
POSTERIOR_THRESHOLD = 0.90


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            if name is None:
                raise ValueError(f"Sequence before header in {path}")
            chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks)
    return records


def exact_clade(tree, descendants: set[str]):
    clade = tree.common_ancestor(sorted(descendants))
    observed = {tip.name for tip in clade.get_terminals()}
    if observed != descendants:
        raise RuntimeError(
            f"Descendants do not define an exact clade: expected={sorted(descendants)}, "
            f"observed={sorted(observed)}"
        )
    return clade


def parse_support(name: str | None) -> tuple[float, float]:
    if not name:
        raise RuntimeError("Candidate edge lacks IQ-TREE support labels")
    parts = name.split("/")
    if len(parts) < 2:
        raise RuntimeError(f"Cannot parse SH-aLRT/UFBoot label {name!r}")
    return float(parts[-2]), float(parts[-1])


def aa_at(sequence: str, site: int) -> str:
    codon = sequence[(site - 1) * 3 : site * 3]
    if codon == "---":
        return "-"
    if len(codon) != 3 or any(base not in "ACGT" for base in codon):
        return "X"
    return str(Seq(codon).translate())


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    primary = read_fasta(ALIGNMENT)
    prank = read_fasta(PRANK_ALIGNMENT)
    if len(primary) != len(set(primary.values())):
        raise RuntimeError("Primary HyPhy alignment still contains exact duplicate sequences")
    if set(primary) != set(prank):
        raise RuntimeError("Primary and PRANK alignment tip sets differ")

    prank_site_rows = [
        row for row in read_tsv(PRANK_CROSSWALK)
        if int(row["human_DNMT3A1_aa"]) == HUMAN_POSITION
    ]
    if len(prank_site_rows) != 1:
        raise RuntimeError("Could not uniquely map human G34 into the PRANK alignment")
    prank_site = int(prank_site_rows[0]["filtered_codon_site_1based"])

    selected = []
    for row in read_tsv(CANDIDATES):
        if (
            int(row["human_DNMT3A1_aa"]) == HUMAN_POSITION
            and row["branch_type"] == "internal"
            and row["reconstructed_change"] == "true"
            and row["change_class"] == "nonsynonymous"
            and float(row["positive_class_posterior"]) >= POSTERIOR_THRESHOLD
            and int(row["composition_failure_descendant_count"]) == 0
        ):
            selected.append(row)
    if not selected:
        raise RuntimeError("No G34 foreground branches passed the declared criteria")

    meme_data = json.loads(MEME_JSON.read_text())
    meme_newick = meme_data["input"]["trees"]["0"]
    meme_tree = Phylo.read(StringIO(meme_newick), "newick")
    iqtree = Phylo.read(IQTREE, "newick")
    all_tips = {tip.name for tip in iqtree.get_terminals()}

    manifest_rows: list[dict[str, object]] = []
    for row in selected:
        branch = row["branch"]
        descendants = set(row["descendant_tips"].split(";"))
        meme_clade = exact_clade(meme_tree, descendants)
        iq_clade = exact_clade(iqtree, descendants)
        if meme_clade.name != branch:
            raise RuntimeError(
                f"Candidate label mismatch for {sorted(descendants)}: "
                f"{branch} versus {meme_clade.name}"
            )
        if meme_clade.branch_length is None or meme_clade.branch_length <= 1e-6:
            raise RuntimeError(f"{branch} is a zero or near-zero branch")
        sh_alrt, ufboot = parse_support(iq_clade.name)
        entire_order_constraint = len(descendants) > 1 and (
            len({tip.name for tip in iq_clade.get_terminals()}) == len(all_tips)
        )
        # The production constraint contains only full mammalian order clades.
        # A candidate nested within an order is therefore not constraint-enforced.
        orders = set(row["descendant_orders"].split(";"))
        same_order = len(orders) == 1
        order_tip_count = sum(
            1
            for candidate in read_tsv(
                ROOT / "results/01_curated/taxon_metadata.tsv"
            )
            if candidate["is_mammal"] == "true"
            and candidate["order"] in orders
            and candidate["safe_id"] in all_tips
        )
        constraint_enforced = same_order and len(descendants) == order_tip_count
        if entire_order_constraint or constraint_enforced:
            raise RuntimeError(f"{branch} is a constraint-enforced order edge")

        primary_states = sorted({aa_at(primary[tip], SITE) for tip in descendants})
        prank_states = sorted({aa_at(prank[tip], prank_site) for tip in descendants})
        prank_concordant = primary_states == prank_states and "-" not in prank_states
        if not prank_concordant:
            raise RuntimeError(f"{branch} G34 descendant state is not PRANK-concordant")

        token = f"{branch}:"
        replacement = f"{branch}{{Foreground}}:"
        if meme_newick.count(token) != 1:
            raise RuntimeError(f"Expected exactly one Newick token for {branch}")
        meme_newick = meme_newick.replace(token, replacement)
        manifest_rows.append({
            "branch": branch,
            "human_site": "G34",
            "descendant_tip_count": len(descendants),
            "descendant_tips": ";".join(sorted(descendants)),
            "reconstructed_change": (
                f"{row['ancestral_aa']}>{row['derived_aa']} "
                f"({row['ancestral_codon']}>{row['derived_codon']})"
            ),
            "meme_positive_class_posterior": row["positive_class_posterior"],
            "branch_length": meme_clade.branch_length,
            "sh_alrt": sh_alrt,
            "ufboot": ufboot,
            "composition_failure_descendants": 0,
            "unique_alignment_exact_duplicates": 0,
            "constraint_enforced_order_edge": "false",
            "primary_descendant_states": ",".join(primary_states),
            "prank_descendant_states": ",".join(prank_states),
            "prank_state_concordant": "true",
            "selection_status": "foreground",
            "selection_rationale": (
                "post_hoc_internal_G34_nonsynonymous_change;"
                "MEME_positive_class_posterior>=0.90;"
                "composition_pass;positive_branch_length;"
                "not_order_constraint_edge;PRANK_state_concordant"
            ),
        })

    TREE_OUT.write_text(meme_newick.rstrip() + "\n")
    annotated = Phylo.read(TREE_OUT, "newick")
    annotations = [
        clade.name for clade in annotated.find_clades()
        if clade.name and clade.name.endswith("{Foreground}")
    ]
    if len(annotations) != len(manifest_rows):
        raise RuntimeError("Foreground annotations were not preserved in output Newick")
    if {tip.name for tip in annotated.get_terminals()} != set(primary):
        raise RuntimeError("Foreground tree and alignment tip sets differ")

    fields = list(manifest_rows[0])
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Prepared {len(manifest_rows)} G34 foreground branches in {TREE_OUT}")


if __name__ == "__main__":
    main()
