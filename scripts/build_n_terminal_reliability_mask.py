#!/usr/bin/env python3
"""Build a three-alignment reliability mask for the DNMT3A N terminus."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALIGNMENTS = {
    "primary": (
        ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
    ),
    "mafft": (
        ROOT
        / "results/02_alignment/sensitivity"
        / "DNMT3A_MAFFT_clean_mammals_HyPhy_unique.fasta"
    ),
    "prank": (
        ROOT
        / "results/02_alignment/sensitivity/prank"
        / "DNMT3A_PRANK_clean_mammals_HyPhy_unique.fasta"
    ),
}
CROSSWALKS = {
    "primary": ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv",
    "mafft": (
        ROOT
        / "results/02_alignment/sensitivity/MAFFT_human_coordinate_crosswalk.tsv"
    ),
    "prank": (
        ROOT
        / "results/02_alignment/sensitivity/prank/PRANK_human_coordinate_crosswalk.tsv"
    ),
}
OUT = ROOT / "results/02_alignment/reliability_mask"
N_TERMINAL_END = 277
CANDIDATES = {5: "P5", 6: "S6", 12: "T12", 18: "E18", 34: "G34", 114: "A114"}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            identifier = line[1:].split()[0]
            records[identifier] = ""
        elif identifier is None:
            raise ValueError(f"Sequence before FASTA header in {path}")
        else:
            records[identifier] += line.upper()
    return records


def write_fasta(path: Path, records: dict[str, str], width: int = 80) -> None:
    with path.open("w") as handle:
        for identifier, sequence in records.items():
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start:start + width] + "\n")


def read_crosswalk(path: Path) -> dict[int, dict[str, str]]:
    with path.open() as handle:
        return {
            int(row["human_DNMT3A1_aa"]): row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["human_DNMT3A1_aa"]
        }


def codon(sequence: str, site: int) -> str:
    return sequence[(site - 1) * 3:site * 3]


def is_informative(value: str) -> bool:
    return len(value) == 3 and all(base in "ACGT" for base in value)


def agreement_pattern(values: dict[str, str]) -> str:
    primary, mafft, prank = (
        values["primary"],
        values["mafft"],
        values["prank"],
    )
    if primary == mafft == prank:
        return "all_three_equal"
    if primary == mafft:
        return "primary_mafft_equal"
    if primary == prank:
        return "primary_prank_equal"
    if mafft == prank:
        return "alternatives_equal_primary_differs"
    return "all_three_different"


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def reliability_class(agreement_fraction: float) -> str:
    if agreement_fraction >= 0.99:
        return "high_ge_0.99"
    if agreement_fraction >= 0.95:
        return "moderate_0.95_to_0.99"
    return "low_lt_0.95"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    alignments = {
        name: read_fasta(path) for name, path in ALIGNMENTS.items()
    }
    maps = {
        name: read_crosswalk(path) for name, path in CROSSWALKS.items()
    }
    tip_sets = {frozenset(records) for records in alignments.values()}
    if len(tip_sets) != 1:
        raise RuntimeError("Alignment tip sets differ")
    if any(len(set(map(len, records.values()))) != 1 for records in alignments.values()):
        raise RuntimeError("At least one FASTA is not rectangular")

    shared_positions = sorted(
        set.intersection(*(set(mapping) for mapping in maps.values()))
        & set(range(1, N_TERMINAL_END + 1))
    )
    if not shared_positions:
        raise RuntimeError("No shared N-terminal human coordinates")

    masked = dict(alignments["primary"])
    manifest: list[dict[str, object]] = []
    site_rows: list[dict[str, object]] = []
    taxon_counts = {
        identifier: Counter() for identifier in alignments["primary"]
    }

    for human_position in shared_positions:
        sites = {
            name: int(mapping[human_position]["filtered_codon_site_1based"])
            for name, mapping in maps.items()
        }
        patterns: Counter[str] = Counter()
        masked_count = 0
        already_missing_count = 0
        informative_consensus_count = 0
        for identifier in alignments["primary"]:
            values = {
                name: codon(records[identifier], sites[name])
                for name, records in alignments.items()
            }
            pattern = agreement_pattern(values)
            patterns[pattern] += 1
            taxon_counts[identifier]["positions_audited"] += 1
            taxon_counts[identifier][pattern] += 1
            if pattern == "all_three_equal":
                if is_informative(values["primary"]):
                    informative_consensus_count += 1
                continue

            taxon_counts[identifier]["discordant_positions"] += 1
            primary_informative = is_informative(values["primary"])
            action = "mask_primary_to_NNN" if primary_informative else "primary_already_missing"
            if primary_informative:
                start = (sites["primary"] - 1) * 3
                sequence = masked[identifier]
                masked[identifier] = sequence[:start] + "NNN" + sequence[start + 3:]
                masked_count += 1
                taxon_counts[identifier]["newly_masked_positions"] += 1
            else:
                already_missing_count += 1
                taxon_counts[identifier]["discordant_but_primary_missing"] += 1
            manifest.append({
                "human_DNMT3A1_aa": human_position,
                "human_site": CANDIDATES.get(human_position, ""),
                "primary_filtered_site": sites["primary"],
                "safe_id": identifier,
                "primary_codon": values["primary"],
                "mafft_codon": values["mafft"],
                "prank_codon": values["prank"],
                "agreement_pattern": pattern,
                "action": action,
                "mask": "NNN" if primary_informative else "",
            })

        taxa = len(alignments["primary"])
        agreement_fraction = patterns["all_three_equal"] / taxa
        site_rows.append({
            "human_DNMT3A1_aa": human_position,
            "human_site": CANDIDATES.get(human_position, ""),
            "primary_filtered_site": sites["primary"],
            "taxa": taxa,
            "all_three_equal_count": patterns["all_three_equal"],
            "all_three_equal_fraction": agreement_fraction,
            "informative_three_way_consensus_count": informative_consensus_count,
            "discordant_count": taxa - patterns["all_three_equal"],
            "newly_masked_informative_count": masked_count,
            "discordant_primary_already_missing_count": already_missing_count,
            "primary_mafft_equal_count": patterns["primary_mafft_equal"],
            "primary_prank_equal_count": patterns["primary_prank_equal"],
            "alternatives_equal_primary_differs_count": patterns[
                "alternatives_equal_primary_differs"
            ],
            "all_three_different_count": patterns["all_three_different"],
            "operational_reliability": reliability_class(agreement_fraction),
        })

    taxon_rows: list[dict[str, object]] = []
    for identifier, counts in taxon_counts.items():
        audited = counts["positions_audited"]
        agreement = audited - counts["discordant_positions"]
        taxon_rows.append({
            "safe_id": identifier,
            "positions_audited": audited,
            "all_three_equal_count": agreement,
            "all_three_equal_fraction": agreement / audited,
            "discordant_positions": counts["discordant_positions"],
            "newly_masked_informative_positions": counts["newly_masked_positions"],
            "discordant_but_primary_missing": counts[
                "discordant_but_primary_missing"
            ],
            "operational_reliability": reliability_class(agreement / audited),
        })
    taxon_rows.sort(
        key=lambda row: (
            float(row["all_three_equal_fraction"]),
            str(row["safe_id"]),
        )
    )

    output_alignment = OUT / "DNMT3A_mammals_HyPhy_unique_N_terminal_consensus_masked.fasta"
    write_fasta(output_alignment, masked)
    write_tsv(OUT / "n_terminal_mask_manifest.tsv", manifest)
    write_tsv(OUT / "site_reliability.tsv", site_rows)
    write_tsv(OUT / "taxon_reliability.tsv", taxon_rows)

    # Verify that every sequence change is an intended codon replacement.
    observed_changes = 0
    for identifier, original in alignments["primary"].items():
        revised = masked[identifier]
        if len(original) != len(revised):
            raise RuntimeError(f"Length changed for {identifier}")
        for start in range(0, len(original), 3):
            before, after = original[start:start + 3], revised[start:start + 3]
            if before != after:
                observed_changes += 1
                if after != "NNN" or not is_informative(before):
                    raise RuntimeError(f"Unexpected mask change for {identifier}")
    expected_changes = sum(
        row["action"] == "mask_primary_to_NNN" for row in manifest
    )
    if observed_changes != expected_changes:
        raise RuntimeError(
            f"Observed {observed_changes} changes, expected {expected_changes}"
        )

    candidate_rows = [
        row for row in site_rows if int(row["human_DNMT3A1_aa"]) in CANDIDATES
    ]
    class_counts = Counter(
        str(row["operational_reliability"]) for row in site_rows
    )
    summary = {
        "region": "human_DNMT3A1_aa_1_277_operational_N_terminal",
        "shared_retained_human_positions": len(shared_positions),
        "first_shared_human_position": shared_positions[0],
        "last_shared_human_position": shared_positions[-1],
        "taxa": len(alignments["primary"]),
        "taxon_site_cells_audited": len(shared_positions) * len(alignments["primary"]),
        "discordant_taxon_site_cells": len(manifest),
        "newly_masked_informative_codons": expected_changes,
        "discordant_cells_already_missing_in_primary": sum(
            row["action"] == "primary_already_missing" for row in manifest
        ),
        "site_reliability_class_counts": dict(class_counts),
        "candidate_sites": candidate_rows,
        "lowest_agreement_taxa": taxon_rows[:20],
        "output_alignment": str(output_alignment.relative_to(ROOT)),
        "mask_rule": (
            "At shared human-mapped N-terminal coordinates, replace an informative "
            "primary codon with NNN whenever primary, MAFFT, and PRANK codons are "
            "not all identical. Pre-existing primary gaps or ambiguous codons remain "
            "missing and are recorded without a sequence change."
        ),
        "inference_boundary": (
            "This mask measures alignment-placement reliability. It does not by "
            "itself validate transcript expression or detect every annotation error."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
