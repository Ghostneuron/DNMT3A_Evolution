#!/usr/bin/env python3
"""Remove exact duplicate DNMT3A sequences for HyPhy without changing inference."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]
ALIGNMENT = ROOT / "results/02_alignment/subsets/DNMT3A_mammals.fasta"
TREE = ROOT / "results/03_tree/mammals/DNMT3A_mammals_tree_HyPhy.nwk"
UNIQUE_ALIGNMENT = ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
UNIQUE_TREE = ROOT / "results/03_tree/mammals/DNMT3A_mammals_tree_HyPhy_unique.nwk"
ROOTED_UNIQUE_TREE = (
    ROOT / "results/03_tree/mammals/DNMT3A_mammals_tree_HyPhy_unique_rooted.nwk"
)
ROOTING_AUDIT = ROOT / "results/03_tree/mammals/rooting_audit.json"
AUDIT = ROOT / "results/04_selection/mammals/deduplication_audit.tsv"
COMPOSITION_FAILURES = ROOT / "results/03_tree/mammals/composition_chi2_failures.tsv"
COMPOSITION_ALIGNMENT = (
    ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique_composition_pass.fasta"
)
COMPOSITION_TREE = (
    ROOT / "results/03_tree/mammals/DNMT3A_mammals_tree_HyPhy_unique_composition_pass.nwk"
)
COMPOSITION_AUDIT = ROOT / "results/04_selection/mammals/composition_sensitivity_exclusions.tsv"
MONOTREME_OUTGROUP = ("Ornithorhynchus_anatinus", "Tachyglossus_aculeatus")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                records.append((identifier, "".join(chunks)))
            identifier, chunks = line[1:].split()[0], []
        else:
            if identifier is None:
                raise ValueError("Sequence before FASTA header")
            chunks.append(line)
    if identifier is not None:
        records.append((identifier, "".join(chunks)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start:start + width] + "\n")


def unrooted_splits(tree) -> set[frozenset[str]]:
    all_tips = {tip.name for tip in tree.get_terminals()}
    splits: set[frozenset[str]] = set()
    for clade in tree.find_clades():
        side = {tip.name for tip in clade.get_terminals()}
        other = all_tips - side
        if len(side) < 2 or len(other) < 2:
            continue
        canonical = side if len(side) < len(other) else other
        if len(side) == len(other):
            canonical = min(side, other, key=lambda values: tuple(sorted(values)))
        splits.add(frozenset(canonical))
    return splits


def main() -> None:
    records = read_fasta(ALIGNMENT)
    by_sequence: dict[str, list[str]] = defaultdict(list)
    for identifier, sequence in records:
        by_sequence[sequence].append(identifier)
    duplicate_groups = [
        sorted(members) for members in by_sequence.values() if len(members) > 1
    ]
    remove_to_keep: dict[str, str] = {}
    for members in duplicate_groups:
        kept = members[0]
        for removed in members[1:]:
            remove_to_keep[removed] = kept

    unique_records = [
        (identifier, sequence)
        for identifier, sequence in records
        if identifier not in remove_to_keep
    ]
    if len(unique_records) != len({sequence for _, sequence in unique_records}):
        raise SystemExit("ERROR: deduplicated alignment still contains duplicate sequences")
    write_fasta(UNIQUE_ALIGNMENT, unique_records)

    tree = Phylo.read(TREE, "newick")
    for removed in sorted(remove_to_keep):
        tree.prune(removed)
    UNIQUE_TREE.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, UNIQUE_TREE, "newick", format_branch_length="%1.10f")
    observed = {tip.name for tip in Phylo.read(UNIQUE_TREE, "newick").get_terminals()}
    expected = {identifier for identifier, _ in unique_records}
    if observed != expected:
        raise SystemExit("ERROR: pruned tree and unique alignment tip sets differ")

    unrooted_tree = Phylo.read(UNIQUE_TREE, "newick")
    original_splits = unrooted_splits(unrooted_tree)
    original_distances = {
        (left, right): unrooted_tree.distance(left, right)
        for index, left in enumerate(sorted(expected))
        for right in sorted(expected)[index + 1:]
    }
    outgroup = unrooted_tree.common_ancestor(MONOTREME_OUTGROUP)
    if {tip.name for tip in outgroup.get_terminals()} != set(MONOTREME_OUTGROUP):
        raise SystemExit("ERROR: monotreme outgroup is not a two-tip clade")
    if outgroup.branch_length is None or outgroup.branch_length <= 0:
        raise SystemExit("ERROR: monotreme stem lacks a positive branch length")
    unrooted_tree.root_with_outgroup(
        outgroup, outgroup_branch_length=outgroup.branch_length / 2
    )
    Phylo.write(
        unrooted_tree, ROOTED_UNIQUE_TREE, "newick",
        format_branch_length="%1.10f",
    )
    rooted_tree = Phylo.read(ROOTED_UNIQUE_TREE, "newick")
    rooted_splits = unrooted_splits(rooted_tree)
    if rooted_splits != original_splits:
        raise SystemExit("ERROR: rerooting changed the unrooted topology")
    max_distance_difference = max(
        abs(rooted_tree.distance(left, right) - distance)
        for (left, right), distance in original_distances.items()
    )
    if max_distance_difference > 1e-8:
        raise SystemExit("ERROR: rerooting changed pairwise tip distances")
    root_child_sizes = sorted(len(clade.get_terminals()) for clade in rooted_tree.root.clades)
    if root_child_sizes != [2, len(expected) - 2]:
        raise SystemExit(f"ERROR: unexpected rooted child sizes {root_child_sizes}")
    ROOTING_AUDIT.write_text(json.dumps({
        "outgroup_tips": list(MONOTREME_OUTGROUP),
        "tip_count": len(expected),
        "unrooted_internal_splits": len(original_splits),
        "unrooted_splits_preserved": True,
        "maximum_pairwise_tip_distance_difference": max_distance_difference,
        "root_child_tip_counts": root_child_sizes,
        "interpretation": "root separates sampled monotremes from therians",
    }, indent=2) + "\n")

    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT.open("w", newline="") as handle:
        fields = ["removed_safe_id", "retained_safe_id", "reason"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for removed, retained in sorted(remove_to_keep.items()):
            writer.writerow({
                "removed_safe_id": removed,
                "retained_safe_id": retained,
                "reason": "exactly_identical_clean_codon_sequence;lexicographic_representative",
            })

    with COMPOSITION_FAILURES.open() as handle:
        composition_failures = {
            row["safe_id"] for row in csv.DictReader(handle, delimiter="\t")
        }
    composition_records = [
        (identifier, sequence)
        for identifier, sequence in unique_records
        if identifier not in composition_failures
    ]
    write_fasta(COMPOSITION_ALIGNMENT, composition_records)
    composition_tree = Phylo.read(TREE, "newick")
    for removed in sorted(set(remove_to_keep) | composition_failures):
        composition_tree.prune(removed)
    Phylo.write(
        composition_tree, COMPOSITION_TREE, "newick",
        format_branch_length="%1.10f",
    )
    composition_tip_ids = {
        tip.name for tip in Phylo.read(COMPOSITION_TREE, "newick").get_terminals()
    }
    if composition_tip_ids != {identifier for identifier, _ in composition_records}:
        raise SystemExit("ERROR: composition-pass tree and alignment tip sets differ")
    with COMPOSITION_AUDIT.open("w", newline="") as handle:
        fields = ["excluded_safe_id", "reason"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for identifier in sorted(composition_failures):
            writer.writerow({
                "excluded_safe_id": identifier,
                "reason": "IQ-TREE nucleotide-composition chi-square p<0.05",
            })
    print(
        f"Prepared {len(unique_records)} unique HyPhy sequences; "
        f"removed {len(remove_to_keep)} redundant copies in {len(duplicate_groups)} groups; "
        f"composition sensitivity retains {len(composition_records)} sequences"
    )


if __name__ == "__main__":
    main()
