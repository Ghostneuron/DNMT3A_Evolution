#!/usr/bin/env python3
"""Remove only DNMT3A EM-seq intermediates with validated successors."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

try:
    from scripts.run_mch_bismark import DEFAULT_EXTERNAL_ROOT
except ModuleNotFoundError:
    from run_mch_bismark import DEFAULT_EXTERNAL_ROOT


SAMTOOLS = Path("/opt/homebrew/bin/samtools")


def validate_bam(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"missing or empty validated successor: {path}")
    subprocess.run([str(SAMTOOLS), "quickcheck", "-v", str(path)], check=True)


def validate_context_summary(sample_root: Path) -> None:
    output = sample_root / "context_summary"
    report_path = output / "context_summary.json"
    table_path = output / "context_summary.tsv"
    bin_path = output / "autosomal_bin_context_summary.tsv.gz"
    gene_path = output / "gene_body_context_summary.tsv.gz"
    required = (report_path, table_path, bin_path, gene_path)
    if any(not path.is_file() or path.stat().st_size == 0 for path in required):
        raise RuntimeError("context summary is incomplete")
    report = json.loads(report_path.read_text())
    if int(report.get("input_rows", 0)) <= 0:
        raise RuntimeError("context summary reports no input rows")
    table = table_path.read_text()
    for label in ("autosome\tCA\t", "autosome\tCC\t", "autosome\tCT\t"):
        if label not in table:
            raise RuntimeError(f"context summary lacks {label.strip()}")


def remove(paths: list[Path]) -> tuple[list[str], int]:
    removed: list[str] = []
    reclaimed = 0
    for path in paths:
        if path.is_file():
            reclaimed += path.stat().st_size
            path.unlink()
            removed.append(str(path))
    return removed, reclaimed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--level", required=True, choices=("aligned", "summarized")
    )
    parser.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    args = parser.parse_args()

    sample_root = args.external_root / "results/bismark" / args.sample
    dedup_bam = (
        sample_root / "deduplicated" / f"{args.sample}_pe.deduplicated.bam"
    )
    validate_bam(dedup_bam)
    candidates = [sample_root / "alignment" / f"{args.sample}_pe.bam"]

    if args.level == "summarized":
        validate_context_summary(sample_root)
        methylation = sample_root / "methylation"
        candidates.extend(
            [
                sample_root / "trimmed" / f"{args.sample}_R1.trimmed.fastq.gz",
                sample_root / "trimmed" / f"{args.sample}_R2.trimmed.fastq.gz",
                methylation / f"CHG_context_{args.sample}_pe.deduplicated.txt.gz",
                methylation / f"CHH_context_{args.sample}_pe.deduplicated.txt.gz",
                methylation / f"CpG_context_{args.sample}_pe.deduplicated.txt.gz",
                methylation / f"{args.sample}_pe.deduplicated.bedGraph.gz",
            ]
        )

    removed, reclaimed = remove(candidates)
    print(
        json.dumps(
            {
                "sample": args.sample,
                "level": args.level,
                "removed": removed,
                "reclaimed_bytes": reclaimed,
                "reclaimed_gib": round(reclaimed / 2**30, 3),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
