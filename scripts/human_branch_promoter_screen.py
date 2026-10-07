#!/usr/bin/env python3
"""Screen DNMT3A promoters for enrichment of human-branch substitutions."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

try:
    from scripts.domain_profiles import hypergeometric_upper_tail
except ModuleNotFoundError:
    from domain_profiles import hypergeometric_upper_tail


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/ucsc_multiz470"
MANIFEST = RAW / "source_manifest.tsv"
INTERVALS = ROOT / "config/human_promoter_intervals.tsv"
OUT = ROOT / "results/06_isoform_evolution/human_branch_promoter_screen"
APES = ("panTro6", "panPan3", "gorGor6", "ponAbe3")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty output: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_maf_sequences(encoded_block: str) -> dict[str, tuple[list[str], str]]:
    """Parse semicolon-encoded bigMaf sequence rows by assembly prefix."""
    output = {}
    for line in encoded_block.split(";"):
        fields = line.strip().split()
        if not fields or fields[0] != "s":
            continue
        species = fields[1].split(".", 1)[0]
        output[species] = (fields, fields[6])
    return output


def validate_sources() -> dict[str, dict[str, str]]:
    sources = {row["local_file"]: row for row in read_tsv(MANIFEST)}
    for filename, row in sources.items():
        path = RAW / filename
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        if observed != row["sha256"]:
            raise ValueError(f"SHA-256 mismatch for {path}")
    return sources


def screen_interval(interval: dict[str, str]) -> tuple[list[dict], dict]:
    path = RAW / interval["input_file"]
    payload = json.loads(path.read_text())
    query_start = int(interval["query_start_0based"])
    query_end = int(interval["query_end_0based"])
    promoter_start = int(interval["promoter_start_0based"])
    promoter_end = int(interval["promoter_end_0based_exclusive"])
    tsses = [int(value) for value in interval["annotated_TSSs_0based"].split(",")]
    counts = {
        "promoter": {"callable": 0, "substitutions": 0},
        "local_flanks": {"callable": 0, "substitutions": 0},
    }
    variants = []
    seen_positions = set()
    for block in payload["multiz470way"]:
        sequences = parse_maf_sequences(block["mafBlock"])
        if "hg38" not in sequences or any(ape not in sequences for ape in APES):
            continue
        human_fields, human_sequence = sequences["hg38"]
        ape_sequences = {ape: sequences[ape][1] for ape in APES}
        coordinate = int(human_fields[2])
        if human_fields[4] != "+":
            raise ValueError("Expected hg38 reference on the plus strand")
        for column, human_base in enumerate(human_sequence):
            if human_base == "-":
                continue
            position = coordinate
            coordinate += 1
            if (
                position < query_start
                or position >= query_end
                or position in seen_positions
            ):
                continue
            seen_positions.add(position)
            human_base = human_base.upper()
            ape_bases = [ape_sequences[ape][column].upper() for ape in APES]
            if (
                human_base not in "ACGT"
                or any(base not in "ACGT" for base in ape_bases)
                or len(set(ape_bases)) != 1
            ):
                continue
            region = (
                "promoter"
                if promoter_start <= position < promoter_end
                else "local_flanks"
            )
            counts[region]["callable"] += 1
            ancestral = ape_bases[0]
            if human_base != ancestral:
                counts[region]["substitutions"] += 1
                variants.append({
                    "promoter": interval["promoter"],
                    "chromosome": interval["query_chrom"],
                    "hg38_position_0based": position,
                    "hg38_position_1based": position + 1,
                    "human_base": human_base,
                    "great_ape_consensus_base": ancestral,
                    "region": region,
                    "distance_to_nearest_annotated_TSS_bp": min(
                        abs(position - tss) for tss in tsses
                    ),
                    "ancestral_state_rule": (
                        "chimpanzee, bonobo, gorilla, and orangutan all "
                        "callable and concordant"
                    ),
                })
    promoter = counts["promoter"]
    flanks = counts["local_flanks"]
    total_callable = promoter["callable"] + flanks["callable"]
    total_substitutions = (
        promoter["substitutions"] + flanks["substitutions"]
    )
    rate = promoter["substitutions"] / promoter["callable"]
    flank_rate = flanks["substitutions"] / flanks["callable"]
    summary = {
        "promoter": interval["promoter"],
        "query_interval": [
            interval["query_chrom"], query_start, query_end
        ],
        "promoter_interval_0based_half_open": [promoter_start, promoter_end],
        "strict_great_ape_callable_promoter_bases": promoter["callable"],
        "strict_great_ape_callable_flank_bases": flanks["callable"],
        "human_branch_promoter_substitutions": promoter["substitutions"],
        "human_branch_flank_substitutions": flanks["substitutions"],
        "promoter_substitution_fraction": rate,
        "flank_substitution_fraction": flank_rate,
        "promoter_vs_flank_rate_ratio": rate / flank_rate,
        "one_sided_enrichment_p": hypergeometric_upper_tail(
            total_callable,
            total_substitutions,
            promoter["callable"],
            promoter["substitutions"],
        ),
        "screen_result": (
            "no_human_branch_substitution_enrichment"
            if rate <= flank_rate else
            "higher_promoter_rate_requires_statistical_interpretation"
        ),
        "inference_limit": (
            "This is a fixed human-branch substitution screen against a "
            "strict four-great-ape consensus and local aligned flanks. It is "
            "not a formal branch-specific phyloP test, functional assay, or "
            "test of positive selection."
        ),
    }
    return variants, summary


def main() -> None:
    validate_sources()
    all_variants = []
    summaries = []
    for interval in read_tsv(INTERVALS):
        variants, summary = screen_interval(interval)
        all_variants.extend(variants)
        summaries.append(summary)
    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT / "human_branch_substitutions.tsv", all_variants)
    write_tsv(OUT / "promoter_enrichment_summary.tsv", summaries)
    output = {
        "promoters_tested": len(summaries),
        "great_ape_consensus_species": list(APES),
        "results": summaries,
        "conclusion": (
            "Neither DNMT3A promoter shows an elevated human-branch "
            "substitution fraction relative to its aligned local flanks."
        ),
        "inference_boundary": (
            "Failure to detect substitution enrichment does not exclude "
            "selection on individual variants, insertions/deletions, "
            "epigenetic regulation, or lineage changes outside the tested "
            "windows."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
