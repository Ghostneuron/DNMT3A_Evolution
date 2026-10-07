#!/usr/bin/env python3
"""Summarize a completed EM-seq sample and start the next queued sample."""

from __future__ import annotations

import argparse
import os
import subprocess
import time
from pathlib import Path

try:
    from scripts.run_mch_bismark import DEFAULT_EXTERNAL_ROOT
except ModuleNotFoundError:
    from run_mch_bismark import DEFAULT_EXTERNAL_ROOT


ROOT = Path(__file__).resolve().parents[1]
PYTHON = Path("/opt/anaconda3/bin/python")
SAMTOOLS = Path("/opt/homebrew/bin/samtools")
PIPELINE = ROOT / "scripts/run_mch_bismark.py"
SUMMARIZE = ROOT / "scripts/summarize_mch_contexts.py"
CLEANUP = ROOT / "scripts/cleanup_mch_intermediates.py"
COMPARE = ROOT / "scripts/compare_mch_genotypes.py"


def process_exists(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def require_nonempty(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"missing or empty required output: {path}")


def validate_sample(external_root: Path, sample: str) -> tuple[Path, Path]:
    sample_root = external_root / "results/bismark" / sample
    dedup = sample_root / "deduplicated" / f"{sample}_pe.deduplicated.bam"
    coverage = (
        sample_root
        / "methylation"
        / f"{sample}_pe.deduplicated.bismark.cov.gz"
    )
    splitting = (
        sample_root
        / "methylation"
        / f"{sample}_pe.deduplicated_splitting_report.txt"
    )
    for path in (dedup, coverage, splitting):
        require_nonempty(path)
    subprocess.run([str(SAMTOOLS), "quickcheck", "-v", str(dedup)], check=True)
    return sample_root, coverage


def run(command: list[str]) -> None:
    print("$ " + " ".join(command), flush=True)
    subprocess.run(command, check=True, cwd=ROOT)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--completed-sample", required=True)
    parser.add_argument("--completed-pid", required=True, type=int)
    parser.add_argument("--next-sample", required=True)
    parser.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--poll-seconds", type=int, default=120)
    args = parser.parse_args()

    print(f"Waiting for sample PID {args.completed_pid}", flush=True)
    while process_exists(args.completed_pid):
        time.sleep(args.poll_seconds)

    sample_root, coverage = validate_sample(
        args.external_root, args.completed_sample
    )
    context = sample_root / "context_summary"
    gene_bed = (
        args.external_root
        / "reference/mm9_annotation/mm9_refgene_longest_protein_coding.bed"
    )
    run(
        [
            str(PYTHON),
            str(SUMMARIZE),
            str(coverage),
            "--external-root",
            str(args.external_root),
            "--output-dir",
            str(context),
            "--gene-bed",
            str(gene_bed),
        ]
    )
    run(
        [
            str(PYTHON),
            str(CLEANUP),
            "--sample",
            args.completed_sample,
            "--level",
            "summarized",
            "--external-root",
            str(args.external_root),
        ]
    )
    run(
        [
            str(PYTHON),
            str(COMPARE),
            "--external-root",
            str(args.external_root),
        ]
    )
    run(
        [
            str(PYTHON),
            str(PIPELINE),
            "all",
            "--external-root",
            str(args.external_root),
            "--sample",
            args.next_sample,
            "--threads",
            str(args.threads),
            "--ignore",
            "10",
            "--ignore-r2",
            "5",
            "--ignore-3prime",
            "10",
            "--ignore-3prime-r2",
            "0",
        ]
    )


if __name__ == "__main__":
    main()
