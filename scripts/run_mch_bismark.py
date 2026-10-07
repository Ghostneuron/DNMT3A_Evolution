#!/usr/bin/env python3
"""Run a resumable Bismark reanalysis of GSE164265 neuronal EM-seq data.

The workflow deliberately preserves CHG and CHH calls separately. It does not
use Bismark's non-conversion filter because methylated CH is the biological
signal being tested in postnatal neurons.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config/mch_samples.tsv"
DEFAULT_EXTERNAL_ROOT = Path(
    os.environ.get("DNMT3A_MCH_ROOT", str(ROOT / "external_data/mCH"))
)
BISMARK = Path("/opt/homebrew/bin/bismark")
BOWTIE2_DIR = Path("/opt/homebrew/bin")
FASTP = Path("/opt/homebrew/bin/fastp")
REFERENCE_RELATIVE = Path("reference/mm9_controls")


def read_samples(path: Path = MANIFEST) -> dict[str, dict[str, str]]:
    """Return one merged manifest record per run accession."""
    samples: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run = row["run_accession"]
            record = samples.setdefault(
                run,
                {
                    "run_accession": run,
                    "sample_accession": row["sample_accession"],
                    "genotype": row["group"],
                    "replicate": row["replicate"],
                },
            )
            record[f"read{row['read_pair']}"] = Path(row["fastq_url"]).name
    for run, record in samples.items():
        if "read1" not in record or "read2" not in record:
            raise ValueError(f"{run} lacks a complete FASTQ pair")
    return samples


def quote_command(command: list[str]) -> str:
    return shlex.join(command)


def run_command(
    command: list[str],
    log_path: Path,
    dry_run: bool = False,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"$ {quote_command(command)}")
    if dry_run:
        return
    with log_path.open("a") as log:
        log.write(f"$ {quote_command(command)}\n")
        log.flush()
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)


def prepare_command(external_root: Path, threads: int) -> list[str]:
    return [
        str(BISMARK),
        "prepare",
        "--bowtie2",
        "--path_to_aligner",
        str(BOWTIE2_DIR),
        "--parallel",
        str(threads),
        "--genomic_composition",
        str(external_root / REFERENCE_RELATIVE),
    ]


def reference_prepared(external_root: Path) -> bool:
    root = external_root / REFERENCE_RELATIVE / "Bisulfite_Genome"
    ct = list((root / "CT_conversion").glob("*.bt2*"))
    ga = list((root / "GA_conversion").glob("*.bt2*"))
    return len(ct) >= 6 and len(ga) >= 6


def align_command(
    external_root: Path,
    sample: dict[str, str],
    threads: int,
    upto: int | None,
) -> list[str]:
    run = sample["run_accession"]
    trimmed_dir = external_root / "results/bismark" / run / "trimmed"
    output_dir = external_root / "results/bismark" / run / "alignment"
    command = [
        str(BISMARK),
        "align",
        "--genome",
        str(external_root / REFERENCE_RELATIVE),
        "--path_to_bowtie2",
        str(BOWTIE2_DIR),
        "--samtools_path",
        str(BOWTIE2_DIR),
        "-p",
        str(threads),
        "--basename",
        run,
        "--output_dir",
        str(output_dir),
        "--temp_dir",
        str(output_dir / "tmp"),
        "-1",
        str(trimmed_dir / f"{run}_R1.trimmed.fastq.gz"),
        "-2",
        str(trimmed_dir / f"{run}_R2.trimmed.fastq.gz"),
    ]
    if upto is not None:
        command[2:2] = ["--upto", str(upto)]
    return command


def trim_command(
    external_root: Path,
    sample: dict[str, str],
    threads: int,
    upto: int | None,
) -> list[str]:
    run = sample["run_accession"]
    fastq_dir = external_root / "raw/GSE164265/fastq"
    output_dir = external_root / "results/bismark" / run / "trimmed"
    command = [
        str(FASTP),
        "--in1",
        str(fastq_dir / sample["read1"]),
        "--in2",
        str(fastq_dir / sample["read2"]),
        "--out1",
        str(output_dir / f"{run}_R1.trimmed.fastq.gz"),
        "--out2",
        str(output_dir / f"{run}_R2.trimmed.fastq.gz"),
        "--detect_adapter_for_pe",
        "--disable_quality_filtering",
        "--disable_length_filtering",
        "--disable_trim_poly_g",
        "--thread",
        str(threads),
        "--json",
        str(output_dir / f"{run}.fastp.json"),
        "--html",
        str(output_dir / f"{run}.fastp.html"),
        "--report_title",
        f"{run} EM-seq adapter trimming",
    ]
    if upto is not None:
        command.extend(["--reads_to_process", str(upto)])
    return command


def dedup_command(
    external_root: Path,
    sample: dict[str, str],
    threads: int,
) -> list[str]:
    run = sample["run_accession"]
    sample_root = external_root / "results/bismark" / run
    return [
        str(BISMARK),
        "dedup",
        "--paired",
        "--parallel",
        str(threads),
        "--output_dir",
        str(sample_root / "deduplicated"),
        str(sample_root / "alignment" / f"{run}_pe.bam"),
    ]


def extract_command(
    external_root: Path,
    sample: dict[str, str],
    threads: int,
    ignore: tuple[int, int, int, int],
) -> list[str]:
    run = sample["run_accession"]
    sample_root = external_root / "results/bismark" / run
    dedup_bam = sample_root / "deduplicated" / f"{run}_pe.deduplicated.bam"
    return [
        str(BISMARK),
        "extract",
        "--paired-end",
        "--comprehensive",
        "--gzip",
        "--no_overlap",
        "--bedGraph",
        "--CX",
        "--counts",
        "--parallel",
        str(threads),
        "--ignore",
        str(ignore[0]),
        "--ignore_r2",
        str(ignore[1]),
        "--ignore_3prime",
        str(ignore[2]),
        "--ignore_3prime_r2",
        str(ignore[3]),
        "--output_dir",
        str(sample_root / "methylation"),
        str(dedup_bam),
    ]


def write_run_metadata(
    external_root: Path,
    sample: dict[str, str],
    args: argparse.Namespace,
) -> None:
    run = sample["run_accession"]
    path = external_root / "results/bismark" / run / "run_parameters.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "sample": sample,
        "stage": args.stage,
        "threads": args.threads,
        "upto_read_pairs": args.upto,
        "extraction_end_trimming": {
            "r1_5prime": args.ignore,
            "r2_5prime": args.ignore_r2,
            "r1_3prime": args.ignore_3prime,
            "r2_3prime": args.ignore_3prime_r2,
        },
        "note": (
            "Extraction trimming must be set from the sample M-bias profile; "
            "zeroes are an initial diagnostic, not a final assumption."
        ),
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "stage", choices=("prepare", "trim", "align", "dedup", "extract", "all")
    )
    parser.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    parser.add_argument("--sample", help="SRR run accession; omit for all samples")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--upto",
        type=int,
        help="Pilot: align only the first N read pairs",
    )
    parser.add_argument("--ignore", type=int, default=0)
    parser.add_argument("--ignore-r2", type=int, default=0)
    parser.add_argument("--ignore-3prime", type=int, default=0)
    parser.add_argument("--ignore-3prime-r2", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.threads < 2:
        parser.error("--threads must be at least 2")
    samples = read_samples()
    if args.sample:
        if args.sample not in samples:
            parser.error(f"unknown sample: {args.sample}")
        selected = [samples[args.sample]]
    else:
        selected = list(samples.values())

    if args.stage in {"prepare", "all"}:
        if args.stage == "all" and reference_prepared(args.external_root):
            print("$ reference index already prepared; skipping preparation")
        else:
            run_command(
                prepare_command(args.external_root, args.threads),
                args.external_root / "logs/bismark_prepare.log",
                args.dry_run,
            )
        if args.stage == "prepare":
            return

    stages = (
        ("trim", "align", "dedup", "extract")
        if args.stage == "all"
        else (args.stage,)
    )
    for sample in selected:
        run = sample["run_accession"]
        if not args.dry_run:
            write_run_metadata(args.external_root, sample, args)
        for stage in stages:
            output_names = {
                "trim": "trimmed",
                "align": "alignment",
                "dedup": "deduplicated",
                "extract": "methylation",
            }
            if not args.dry_run:
                (
                    args.external_root
                    / "results/bismark"
                    / run
                    / output_names[stage]
                ).mkdir(parents=True, exist_ok=True)
            if stage == "trim":
                command = trim_command(
                    args.external_root, sample, args.threads, args.upto
                )
            elif stage == "align":
                if not args.dry_run:
                    (
                        args.external_root
                        / "results/bismark"
                        / run
                        / "alignment"
                        / "tmp"
                    ).mkdir(parents=True, exist_ok=True)
                command = align_command(
                    args.external_root, sample, args.threads, args.upto
                )
            elif stage == "dedup":
                command = dedup_command(args.external_root, sample, args.threads)
            else:
                ignore = (
                    args.ignore,
                    args.ignore_r2,
                    args.ignore_3prime,
                    args.ignore_3prime_r2,
                )
                command = extract_command(
                    args.external_root, sample, args.threads, ignore
                )
            run_command(
                command,
                args.external_root / "logs" / f"{run}.{stage}.log",
                args.dry_run,
            )


if __name__ == "__main__":
    main()
