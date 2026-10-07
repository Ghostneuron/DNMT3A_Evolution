#!/usr/bin/env python3
"""Audit the production DNMT3A mammalian IQ-TREE result."""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
from collections import defaultdict
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PREFIX = ROOT / "results/03_tree/mammals/DNMT3A_mammals_tree"
TAXONOMY = ROOT / "results/01_curated/taxon_metadata.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, default=DEFAULT_PREFIX)
    args = parser.parse_args()
    prefix = args.prefix if args.prefix.is_absolute() else ROOT / args.prefix
    tree_path = Path(str(prefix) + ".treefile")
    log_path = Path(str(prefix) + ".log")
    report_path = Path(str(prefix) + ".iqtree")
    if not tree_path.exists() or not log_path.exists() or not report_path.exists():
        raise SystemExit(f"ERROR: incomplete IQ-TREE output for prefix {prefix}")

    tree = Phylo.read(tree_path, "newick")
    terminals = tree.get_terminals()
    terminal_by_name = {tip.name: tip for tip in terminals}
    if len(terminal_by_name) != len(terminals):
        raise SystemExit("ERROR: duplicate or unnamed terminal labels in tree")

    metadata = [row for row in read_tsv(TAXONOMY) if row["is_mammal"] == "true"]
    expected = {row["safe_id"] for row in metadata}
    observed = set(terminal_by_name)
    if expected != observed:
        raise SystemExit(
            f"ERROR: tree/taxonomy tip mismatch; missing={sorted(expected-observed)[:10]}, "
            f"extra={sorted(observed-expected)[:10]}"
        )

    orders: dict[str, list[str]] = defaultdict(list)
    for row in metadata:
        orders[row["order"]].append(row["safe_id"])
    all_tips = set(terminal_by_name)
    edge_splits: list[set[str]] = []
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        edge_splits.append({tip.name for tip in clade.get_terminals()})
    order_rows: list[dict[str, object]] = []
    for order, names in sorted(orders.items()):
        target = set(names)
        monophyletic = (
            len(target) == 1
            or target in edge_splits
            or (all_tips - target) in edge_splits
        )
        order_rows.append({
            "order": order,
            "taxa": len(names),
            "monophyletic": str(monophyletic).lower(),
        })

    terminal_rows = sorted(
        (
            {
                "safe_id": tip.name,
                "terminal_branch_length": tip.branch_length or 0.0,
                "order": next(row["order"] for row in metadata if row["safe_id"] == tip.name),
            }
            for tip in terminals
        ),
        key=lambda row: float(row["terminal_branch_length"]),
        reverse=True,
    )

    composition_rows: list[dict[str, object]] = []
    for line in log_path.read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 5 and fields[0].isdigit() and "failed" in fields:
            composition_rows.append({
                "safe_id": fields[1],
                "gap_or_ambiguity": fields[2],
                "p_value": fields[-1],
            })

    constraint_labels = set(orders) | {"Mammalia"}
    support_rows: list[dict[str, object]] = []
    for clade in tree.get_nonterminals():
        if not clade.name:
            continue
        parts = clade.name.split("/")
        if len(parts) < 2:
            continue
        try:
            sh_alrt = float(parts[-2])
            ufboot = float(parts[-1])
        except ValueError:
            continue
        label = "/".join(parts[:-2])
        support_rows.append({
            "node_label": label,
            "constraint_node": str(label in constraint_labels).lower(),
            "descendant_tips": len(clade.get_terminals()),
            "branch_length": clade.branch_length or 0.0,
            "sh_alrt": sh_alrt,
            "ufboot": ufboot,
        })
    unconstrained_support = [
        row for row in support_rows if row["constraint_node"] != "true"
    ]
    report_text = report_path.read_text(errors="replace")
    near_zero_match = re.search(r"WARNING: (\d+) near-zero internal branches", report_text)

    branch_lengths = [
        float(clade.branch_length)
        for clade in tree.find_clades()
        if clade.branch_length is not None
    ]
    summary = {
        "tree": str(tree_path.relative_to(ROOT)),
        "tips": len(terminals),
        "total_branch_length": tree.total_branch_length(),
        "minimum_length_terminal_branches_le_1e-6": sum(
            (tip.branch_length or 0.0) <= 1.0000001e-6 for tip in terminals
        ),
        "maximum_terminal_branch_length": terminal_rows[0]["terminal_branch_length"],
        "median_terminal_branch_length": statistics.median(
            float(row["terminal_branch_length"]) for row in terminal_rows
        ),
        "orders": len(order_rows),
        "nonmonophyletic_orders": [
            row["order"] for row in order_rows if row["monophyletic"] != "true"
        ],
        "composition_chi2_failures": len(composition_rows),
        "iqtree_reported_near_zero_internal_branches_lt_0.0004": (
            int(near_zero_match.group(1)) if near_zero_match else None
        ),
        "unconstrained_internal_branches_with_ufboot_lt_70": sum(
            float(row["ufboot"]) < 70 for row in unconstrained_support
        ),
        "unconstrained_internal_branches_with_ufboot_ge_95": sum(
            float(row["ufboot"]) >= 95 for row in unconstrained_support
        ),
        "all_branch_lengths_nonnegative": all(length >= 0 for length in branch_lengths),
    }
    output_dir = prefix.parent
    hyphy_tree_path = Path(str(prefix) + "_HyPhy.nwk")
    newick = tree_path.read_text().strip()
    newick = re.sub(r"\)([^():,;]+)(?=:)", ")", newick)
    newick = re.sub(r"\)([^():,;]+)(?=;)", ")", newick)
    hyphy_tree_path.write_text(newick + "\n")
    hyphy_tree = Phylo.read(hyphy_tree_path, "newick")
    if {tip.name for tip in hyphy_tree.get_terminals()} != observed:
        raise SystemExit("ERROR: HyPhy tree preparation changed terminal labels")
    hyphy_splits = {
        frozenset(tip.name for tip in clade.get_terminals())
        for clade in hyphy_tree.find_clades()
        if clade is not hyphy_tree.root
    }
    if hyphy_splits != {frozenset(split) for split in edge_splits}:
        raise SystemExit("ERROR: HyPhy tree preparation changed topology")
    if abs(hyphy_tree.total_branch_length() - tree.total_branch_length()) > 1e-10:
        raise SystemExit("ERROR: HyPhy tree preparation changed branch lengths")
    (output_dir / "tree_qc_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_tsv(output_dir / "order_monophyly.tsv", order_rows, ["order", "taxa", "monophyletic"])
    write_tsv(
        output_dir / "longest_terminal_branches.tsv",
        terminal_rows[:25],
        ["safe_id", "terminal_branch_length", "order"],
    )
    write_tsv(
        output_dir / "composition_chi2_failures.tsv",
        composition_rows,
        ["safe_id", "gap_or_ambiguity", "p_value"],
    )
    write_tsv(
        output_dir / "branch_support_qc.tsv",
        support_rows,
        [
            "node_label", "constraint_node", "descendant_tips",
            "branch_length", "sh_alrt", "ufboot",
        ],
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
