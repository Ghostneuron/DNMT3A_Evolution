#!/usr/bin/env python3
"""Compare matched DNMT3A1- and DNMT3A2-knockout neuronal methylomes.

The GSE164265 processed BEDGRAPH files contain one methylation fraction per
covered genomic position. This analysis uses only positions present in all six
samples, so group contrasts cannot be driven by unequal site coverage. Sites
are descriptive units, not biological replicates; no site-level P values are
reported.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


ROOT = Path(__file__).resolve().parents[1]
RAW = (
    ROOT
    / "data/raw/isoform_functional_genomics/GSE164265"
    / "neuron_nuclei_methylation"
)
MANIFEST = RAW / "sample_manifest.tsv"
RESULTS = (
    ROOT
    / "results/07_isoform_function"
    / "GSE164265_neuron_nuclei_methylation"
)
GROUP_ORDER = ("WT", "Dnmt3a1_KO", "Dnmt3a2_KO")
WINDOW_SIZE = 100_000
LOSS_THRESHOLD = 0.10
MIN_WINDOW_SITES = 50


def chromosome_key(chromosome: str) -> str:
    """Ordering used by the deposited GSE164265 BEDGRAPH files.

    The files are lexicographically ordered (chr1, chr10, ..., chr19, chr2),
    rather than in natural chromosome order. Preserving the deposited order
    permits a streaming intersection without rewriting ~108 million records.
    """
    return chromosome


def coordinate_key(chromosome: str, start: int) -> tuple[str, int]:
    return chromosome_key(chromosome), start


@dataclass
class BedGraphStream:
    path: Path
    handle: TextIO | None = None
    chromosome: str | None = None
    start: int | None = None
    value: float | None = None
    previous_key: tuple[str, int] | None = None
    records: int = 0

    def open(self) -> None:
        self.handle = gzip.open(self.path, "rt")
        self.advance()

    @property
    def key(self) -> tuple[tuple[int, int | str], int] | None:
        if self.chromosome is None or self.start is None:
            return None
        return coordinate_key(self.chromosome, self.start)

    def advance(self) -> None:
        if self.handle is None:
            raise RuntimeError("stream is not open")
        for line in self.handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 4:
                continue
            chromosome = fields[0]
            start = int(fields[1])
            value = float(fields[3])
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"methylation value outside [0,1] in {self.path}")
            key = coordinate_key(chromosome, start)
            if self.previous_key is not None and key <= self.previous_key:
                raise ValueError(f"BEDGRAPH is not strictly sorted: {self.path}")
            self.previous_key = key
            self.chromosome = chromosome
            self.start = start
            self.value = value
            self.records += 1
            return
        self.chromosome = None
        self.start = None
        self.value = None

    def close(self) -> None:
        if self.handle is not None:
            self.handle.close()


@dataclass
class OnlinePair:
    n: int = 0
    sum_x: float = 0.0
    sum_y: float = 0.0
    sum_x2: float = 0.0
    sum_y2: float = 0.0
    sum_xy: float = 0.0

    def add(self, x: float, y: float) -> None:
        self.n += 1
        self.sum_x += x
        self.sum_y += y
        self.sum_x2 += x * x
        self.sum_y2 += y * y
        self.sum_xy += x * y

    def correlation(self) -> float:
        numerator = self.n * self.sum_xy - self.sum_x * self.sum_y
        denominator = math.sqrt(
            (self.n * self.sum_x2 - self.sum_x**2)
            * (self.n * self.sum_y2 - self.sum_y**2)
        )
        return numerator / denominator if denominator else float("nan")


def read_manifest(path: Path = MANIFEST) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if [row["group"] for row in rows] != [
        "WT", "WT", "Dnmt3a1_KO", "Dnmt3a1_KO", "Dnmt3a2_KO", "Dnmt3a2_KO"
    ]:
        raise ValueError("manifest must contain two replicates for each group")
    return rows


def validate_inputs(rows: list[dict[str, str]], raw_dir: Path = RAW) -> None:
    for row in rows:
        path = raw_dir / row["local_file"]
        if not path.exists():
            raise FileNotFoundError(path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != row["sha256"]:
            raise ValueError(f"SHA-256 mismatch for {path}")


def classify_site(
    wt: tuple[float, float],
    a1: tuple[float, float],
    a2: tuple[float, float],
    threshold: float = LOSS_THRESHOLD,
) -> tuple[str, bool, bool]:
    """Return mean-effect class and replicate-ordered loss flags."""
    wt_mean = sum(wt) / 2
    a1_loss = wt_mean - sum(a1) / 2 >= threshold
    a2_loss = wt_mean - sum(a2) / 2 >= threshold
    if a1_loss and a2_loss:
        effect_class = "shared_loss"
    elif a1_loss:
        effect_class = "Dnmt3a1_KO_only_loss"
    elif a2_loss:
        effect_class = "Dnmt3a2_KO_only_loss"
    else:
        effect_class = "no_mean_loss_ge_threshold"
    a1_ordered = max(a1) <= min(wt) - threshold
    a2_ordered = max(a2) <= min(wt) - threshold
    return effect_class, a1_ordered, a2_ordered


def iter_common_sites(streams: list[BedGraphStream]):
    """Yield coordinates and values present in every sorted input stream."""
    for stream in streams:
        stream.open()
    try:
        while all(stream.key is not None for stream in streams):
            target = max(stream.key for stream in streams if stream.key is not None)
            for stream in streams:
                while stream.key is not None and stream.key < target:
                    stream.advance()
            if any(stream.key is None for stream in streams):
                break
            keys = [stream.key for stream in streams]
            if len(set(keys)) == 1:
                yield (
                    streams[0].chromosome,
                    streams[0].start,
                    tuple(float(stream.value) for stream in streams),
                )
                for stream in streams:
                    stream.advance()
    finally:
        for stream in streams:
            stream.close()


def analyze(
    rows: list[dict[str, str]],
    raw_dir: Path = RAW,
    window_size: int = WINDOW_SIZE,
) -> dict:
    streams = [BedGraphStream(raw_dir / row["local_file"]) for row in rows]
    sample_sums = [0.0] * 6
    pair_stats = [OnlinePair(), OnlinePair(), OnlinePair()]
    effect_counts: defaultdict[str, int] = defaultdict(int)
    ordered_counts: defaultdict[str, int] = defaultdict(int)
    windows: dict[tuple[str, int], list[float]] = {}
    common_sites = 0

    for chromosome, start, values in iter_common_sites(streams):
        common_sites += 1
        for index, value in enumerate(values):
            sample_sums[index] += value
        pair_stats[0].add(values[0], values[1])
        pair_stats[1].add(values[2], values[3])
        pair_stats[2].add(values[4], values[5])

        effect_class, a1_ordered, a2_ordered = classify_site(
            (values[0], values[1]),
            (values[2], values[3]),
            (values[4], values[5]),
        )
        effect_counts[effect_class] += 1
        if a1_ordered and a2_ordered:
            ordered_counts["both_KO_ordered_loss"] += 1
        elif a1_ordered:
            ordered_counts["Dnmt3a1_KO_ordered_loss_only"] += 1
        elif a2_ordered:
            ordered_counts["Dnmt3a2_KO_ordered_loss_only"] += 1
        else:
            ordered_counts["no_ordered_loss"] += 1

        window_start = (start // window_size) * window_size
        key = (str(chromosome), window_start)
        aggregate = windows.setdefault(key, [0.0] * 7)
        aggregate[0] += 1
        for index, value in enumerate(values, start=1):
            aggregate[index] += value

    if common_sites == 0:
        raise ValueError("no genomic positions are shared by all samples")

    sample_means = [
        {
            "sample_accession": row["sample_accession"],
            "group": row["group"],
            "replicate": int(row["replicate"]),
            "common_sites": common_sites,
            "mean_methylation": sample_sums[index] / common_sites,
            "source_records": streams[index].records,
        }
        for index, row in enumerate(rows)
    ]
    group_means = {
        group: sum(
            row["mean_methylation"]
            for row in sample_means
            if row["group"] == group
        )
        / 2
        for group in GROUP_ORDER
    }

    window_rows = []
    for (chromosome, window_start), aggregate in windows.items():
        count = int(aggregate[0])
        if count < MIN_WINDOW_SITES:
            continue
        means = [value / count for value in aggregate[1:]]
        wt = sum(means[0:2]) / 2
        a1 = sum(means[2:4]) / 2
        a2 = sum(means[4:6]) / 2
        window_rows.append({
            "chromosome": chromosome,
            "start": window_start,
            "end": window_start + window_size,
            "common_sites": count,
            "WT_mean": wt,
            "Dnmt3a1_KO_mean": a1,
            "Dnmt3a2_KO_mean": a2,
            "Dnmt3a1_KO_minus_WT": a1 - wt,
            "Dnmt3a2_KO_minus_WT": a2 - wt,
            "KO_difference_A1_minus_A2": a1 - a2,
        })

    return {
        "common_sites": common_sites,
        "sample_means": sample_means,
        "group_means": group_means,
        "group_effects": {
            "Dnmt3a1_KO_minus_WT": group_means["Dnmt3a1_KO"] - group_means["WT"],
            "Dnmt3a2_KO_minus_WT": group_means["Dnmt3a2_KO"] - group_means["WT"],
            "Dnmt3a1_KO_minus_Dnmt3a2_KO": (
                group_means["Dnmt3a1_KO"] - group_means["Dnmt3a2_KO"]
            ),
        },
        "replicate_correlations": {
            group: pair_stats[index].correlation()
            for index, group in enumerate(GROUP_ORDER)
        },
        "mean_effect_class_counts": dict(effect_counts),
        "replicate_ordered_loss_counts": dict(ordered_counts),
        "window_rows": window_rows,
        "source_records": {
            row["sample_accession"]: streams[index].records
            for index, row in enumerate(rows)
        },
    }


def write_results(result: dict, output_dir: Path = RESULTS) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    with (output_dir / "sample_qc.tsv").open("w", newline="") as handle:
        fields = list(result["sample_means"][0])
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(result["sample_means"])

    effect_rows = []
    for group in GROUP_ORDER:
        effect_rows.append({
            "metric": f"{group}_mean_methylation",
            "value": result["group_means"][group],
            "unit": "fraction",
        })
    for metric, value in result["group_effects"].items():
        effect_rows.append({"metric": metric, "value": value, "unit": "fraction"})
    for metric, value in result["replicate_correlations"].items():
        effect_rows.append({
            "metric": f"{metric}_replicate_Pearson_r",
            "value": value,
            "unit": "correlation",
        })
    with (output_dir / "global_effects.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=("metric", "value", "unit")
        )
        writer.writeheader()
        writer.writerows(effect_rows)

    class_rows = []
    for scheme in ("mean_effect_class_counts", "replicate_ordered_loss_counts"):
        for label, count in sorted(result[scheme].items()):
            class_rows.append({
                "classification_scheme": scheme,
                "class": label,
                "sites": count,
                "fraction_of_common_sites": count / result["common_sites"],
            })
    with (output_dir / "site_effect_classes.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=(
                "classification_scheme",
                "class",
                "sites",
                "fraction_of_common_sites",
            ),
        )
        writer.writeheader()
        writer.writerows(class_rows)

    window_fields = (
        "chromosome",
        "start",
        "end",
        "common_sites",
        "WT_mean",
        "Dnmt3a1_KO_mean",
        "Dnmt3a2_KO_mean",
        "Dnmt3a1_KO_minus_WT",
        "Dnmt3a2_KO_minus_WT",
        "KO_difference_A1_minus_A2",
    )
    with gzip.open(output_dir / "window_effects_100kb.tsv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=window_fields)
        writer.writeheader()
        writer.writerows(result["window_rows"])

    ranked = sorted(
        result["window_rows"],
        key=lambda row: abs(row["KO_difference_A1_minus_A2"]),
        reverse=True,
    )[:200]
    with (output_dir / "top_isoform_differential_windows.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=window_fields)
        writer.writeheader()
        writer.writerows(ranked)

    summary = {key: value for key, value in result.items() if key != "window_rows"}
    summary.update({
        "study": "GSE164265",
        "design": (
            "P21 mouse cortical neuron nuclei; two WT, two Dnmt3a1-KO, and "
            "two Dnmt3a2-KO biological replicates"
        ),
        "site_universe": "positions represented in all six processed BEDGRAPH files",
        "loss_threshold": LOSS_THRESHOLD,
        "window_size_bp": WINDOW_SIZE,
        "minimum_common_sites_per_reported_window": MIN_WINDOW_SITES,
        "inference_limit": (
            "Genomic sites are not biological replicates. Effects are descriptive "
            "and require biological-replicate or region-level modeling for formal "
            "inference. BEDGRAPH methylation fractions lack read-depth columns."
        ),
    })
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main() -> None:
    rows = read_manifest()
    validate_inputs(rows)
    result = analyze(rows)
    write_results(result)
    print(json.dumps({
        "common_sites": result["common_sites"],
        "group_means": result["group_means"],
        "group_effects": result["group_effects"],
        "replicate_correlations": result["replicate_correlations"],
    }, indent=2))


if __name__ == "__main__":
    main()
