#!/usr/bin/env python3
"""Score DNMT3A1-tail variants against recurrent natural alternatives."""

from __future__ import annotations

import csv
import json
import math
import re
import statistics
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta  # noqa: E402


HUMAN_PROTEINS = ROOT / "results/01_curated/DNMT3A_representative_proteins.faa"
ALIGNMENT = ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
CROSSWALK = ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv"
VARIANT_PANEL = ROOT / "results/07_functional_prioritization/DNMT3A1_variant_panel.tsv"
OUT = ROOT / "results/07_functional_prioritization/in_silico_idr"
TAIL_START = 1
TAIL_END = 163
MIN_ALTERNATIVE_COUNT = 2
WINDOW_SIZES = (21, 41)

CHARGE = {
    "D": -1.0,
    "E": -1.0,
    "K": 1.0,
    "R": 1.0,
    "H": 0.1,
}
HYDROPATHY = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}
METRICS = (
    "net_charge",
    "fraction_charged",
    "net_charge_per_residue",
    "sequence_charge_decoration",
    "mean_hydropathy",
    "normalized_sequence_entropy",
    "phospho_acceptor_fraction",
    "glycine_proline_fraction",
    "proline_directed_STP_count",
    "RxxST_count",
    "STxxDE_count",
)


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


def human_sequence() -> str:
    records = [
        record for record in read_fasta(HUMAN_PROTEINS)
        if record.identifier == "Homo_sapiens"
    ]
    if len(records) != 1:
        raise RuntimeError(f"Expected one human protein, found {len(records)}")
    sequence = records[0].sequence.rstrip("*")
    if len(sequence) != 912:
        raise RuntimeError(f"Expected 912-aa DNMT3A1, found {len(sequence)} aa")
    return sequence


def normalized_entropy(sequence: str) -> float:
    counts = Counter(sequence)
    total = len(sequence)
    if total == 0 or len(counts) == 1:
        return 0.0
    raw = -sum(
        (count / total) * math.log(count / total)
        for count in counts.values()
    )
    return raw / math.log(20)


def sequence_charge_decoration(sequence: str) -> float:
    """Normalized sequence charge decoration for a local sequence window."""
    charges = [CHARGE.get(residue, 0.0) for residue in sequence]
    if len(charges) < 2:
        return 0.0
    score = sum(
        charges[left] * charges[right] * math.sqrt(right - left)
        for left in range(len(charges))
        for right in range(left + 1, len(charges))
    )
    return score / len(charges)


def sequence_metrics(sequence: str) -> dict[str, float]:
    charges = [CHARGE.get(residue, 0.0) for residue in sequence]
    return {
        "net_charge": sum(charges),
        "fraction_charged": (
            sum(residue in CHARGE for residue in sequence) / len(sequence)
        ),
        "net_charge_per_residue": sum(charges) / len(sequence),
        "sequence_charge_decoration": sequence_charge_decoration(sequence),
        "mean_hydropathy": statistics.fmean(
            HYDROPATHY[residue] for residue in sequence
        ),
        "normalized_sequence_entropy": normalized_entropy(sequence),
        "phospho_acceptor_fraction": (
            sum(residue in "STY" for residue in sequence) / len(sequence)
        ),
        "glycine_proline_fraction": (
            sum(residue in "GP" for residue in sequence) / len(sequence)
        ),
        "proline_directed_STP_count": float(
            len(re.findall(r"(?=[ST]P)", sequence))
        ),
        "RxxST_count": float(len(re.findall(r"(?=R..[ST])", sequence))),
        "STxxDE_count": float(len(re.findall(r"(?=[ST]..[DE])", sequence))),
    }


def window(sequence: str, position: int, size: int) -> tuple[int, int, str]:
    if size % 2 != 1:
        raise ValueError("Window size must be odd")
    radius = size // 2
    start = max(1, position - radius)
    end = min(len(sequence), position + radius)
    return start, end, sequence[start - 1:end]


def mutate(sequence: str, position: int, reference: str, alternative: str) -> str:
    observed = sequence[position - 1]
    if observed != reference:
        raise RuntimeError(
            f"Reference mismatch at {position}: expected {reference}, found {observed}"
        )
    return sequence[:position - 1] + alternative + sequence[position:]


def alternative_states(
    reference_sequence: str,
) -> list[dict[str, object]]:
    records = read_fasta(ALIGNMENT)
    crosswalk = {
        int(row["filtered_codon_site_1based"]): row
        for row in read_tsv(CROSSWALK)
    }
    rows: list[dict[str, object]] = []
    for site, mapped in crosswalk.items():
        position = int(mapped["human_DNMT3A1_aa"])
        if not TAIL_START <= position <= TAIL_END:
            continue
        reference = reference_sequence[position - 1]
        counts: Counter[str] = Counter()
        for record in records:
            codon = record.sequence[(site - 1) * 3:site * 3]
            residue = "-" if codon == "---" else CODON_TABLE.get(codon, "X")
            if residue not in {"-", "X", "*"}:
                counts[residue] += 1
        for alternative, count in sorted(counts.items()):
            if alternative == reference or count < MIN_ALTERNATIVE_COUNT:
                continue
            rows.append({
                "variant": f"{reference}{position}{alternative}",
                "position": position,
                "reference_residue": reference,
                "alternative_residue": alternative,
                "alternative_taxon_count": count,
                "non_gap_taxa": sum(counts.values()),
                "alternative_frequency": count / sum(counts.values()),
                "null_eligible": True,
            })
    return rows


def candidate_variants() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_tsv(VARIANT_PANEL):
        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", row["experimental_substitution"])
        if not match:
            raise RuntimeError(
                f"Unrecognized variant: {row['experimental_substitution']}"
            )
        reference, position, alternative = match.groups()
        rows.append({
            "variant": row["experimental_substitution"],
            "position": int(position),
            "reference_residue": reference,
            "alternative_residue": alternative,
            "priority_tier": int(row["priority_tier"]),
            "change_event_count": int(row["change_event_count"]),
            "supporting_lineages": row["supporting_lineages"],
        })
    return rows


def score_variant(
    reference_sequence: str,
    variant: dict[str, object],
    source: str,
) -> list[dict[str, object]]:
    position = int(variant["position"])
    reference = str(variant["reference_residue"])
    alternative = str(variant["alternative_residue"])
    altered_sequence = mutate(
        reference_sequence, position, reference, alternative
    )
    rows: list[dict[str, object]] = []
    for size in WINDOW_SIZES:
        start, end, reference_window = window(reference_sequence, position, size)
        _, _, altered_window = window(altered_sequence, position, size)
        reference_metrics = sequence_metrics(reference_window)
        altered_metrics = sequence_metrics(altered_window)
        for metric in METRICS:
            delta = altered_metrics[metric] - reference_metrics[metric]
            rows.append({
                "variant": variant["variant"],
                "source": source,
                "position": position,
                "reference_residue": reference,
                "alternative_residue": alternative,
                "window_size": size,
                "window_start": start,
                "window_end": end,
                "metric": metric,
                "reference_value": reference_metrics[metric],
                "alternative_value": altered_metrics[metric],
                "delta": delta,
                "absolute_delta": abs(delta),
            })
    return rows


def bh_adjust(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg adjustment preserving input order."""
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [1.0] * len(p_values)
    running = 1.0
    total = len(p_values)
    for rank_index in range(total - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * total / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def empirical_effect_summary(
    candidate_scores: list[dict[str, object]],
    null_scores: list[dict[str, object]],
    candidate_metadata: dict[str, dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    null_by_metric: dict[tuple[int, str], list[dict[str, object]]] = {}
    for row in null_scores:
        key = (int(row["window_size"]), str(row["metric"]))
        null_by_metric.setdefault(key, []).append(row)

    detail_rows: list[dict[str, object]] = []
    percentiles_by_variant: dict[str, list[tuple[float, float, int]]] = {}
    for row in candidate_scores:
        key = (int(row["window_size"]), str(row["metric"]))
        null = null_by_metric[key]
        reference = str(row["reference_residue"])
        matched = [
            item for item in null
            if item["reference_residue"] == reference
            and item["variant"] != row["variant"]
        ]
        observed = float(row["absolute_delta"])

        def compare(pool: list[dict[str, object]]) -> tuple[float, float]:
            values = [float(item["absolute_delta"]) for item in pool]
            # Use a strict percentile so a zero effect is not ranked highly
            # merely because many null variants also leave the metric unchanged.
            percentile = sum(value < observed for value in values) / len(values)
            upper_p = (1 + sum(value >= observed for value in values)) / (
                1 + len(values)
            )
            return percentile, upper_p

        global_percentile, global_p = compare(null)
        if matched:
            matched_percentile, matched_p = compare(matched)
        else:
            matched_percentile, matched_p = float("nan"), float("nan")
        detail_rows.append({
            **row,
            "global_null_variants": len(null),
            "global_effect_percentile": global_percentile,
            "global_empirical_upper_p": global_p,
            "reference_residue_matched_null_variants": len(matched),
            "matched_effect_percentile": matched_percentile,
            "matched_empirical_upper_p": matched_p,
        })
        percentiles_by_variant.setdefault(str(row["variant"]), []).append(
            (global_percentile, matched_percentile, len(matched))
        )

    detail_by_variant: dict[str, list[dict[str, object]]] = {}
    for row in detail_rows:
        detail_by_variant.setdefault(str(row["variant"]), []).append(row)
    for rows in detail_by_variant.values():
        global_q = bh_adjust(
            [float(row["global_empirical_upper_p"]) for row in rows]
        )
        matched_indices = [
            index for index, row in enumerate(rows)
            if not math.isnan(float(row["matched_empirical_upper_p"]))
        ]
        matched_q_values = bh_adjust([
            float(rows[index]["matched_empirical_upper_p"])
            for index in matched_indices
        ])
        matched_q_by_index = dict(zip(matched_indices, matched_q_values))
        for index, row in enumerate(rows):
            row["global_BH_q_within_variant"] = global_q[index]
            row["matched_BH_q_within_variant"] = matched_q_by_index.get(
                index, ""
            )

    summary_rows: list[dict[str, object]] = []
    for variant, values in percentiles_by_variant.items():
        metadata = candidate_metadata[variant]
        variant_detail = detail_by_variant[variant]
        global_values = [item[0] for item in values]
        matched_values = [
            item[1] for item in values if not math.isnan(item[1])
        ]
        reference = str(metadata["reference_residue"])
        alternative = str(metadata["alternative_residue"])
        if reference in "STY" and alternative not in "STY":
            phospho_change = "removes_STY_acceptor"
        elif reference not in "STY" and alternative in "STY":
            phospho_change = "adds_STY_acceptor"
        else:
            phospho_change = "retains_STY_acceptor_class"
        summary_rows.append({
            "variant": variant,
            "priority_tier": metadata["priority_tier"],
            "position": metadata["position"],
            "reference_residue": reference,
            "alternative_residue": alternative,
            "change_event_count": metadata["change_event_count"],
            "phospho_acceptor_class_change": phospho_change,
            "mean_global_effect_percentile": statistics.fmean(global_values),
            "maximum_global_effect_percentile": max(global_values),
            "metrics_at_or_above_95th_global_percentile": sum(
                value >= 0.95 for value in global_values
            ),
            "mean_reference_matched_effect_percentile": (
                statistics.fmean(matched_values) if matched_values else ""
            ),
            "maximum_reference_matched_effect_percentile": (
                max(matched_values) if matched_values else ""
            ),
            "metrics_with_reference_matched_null": len(matched_values),
            "minimum_global_empirical_upper_p": min(
                float(row["global_empirical_upper_p"])
                for row in variant_detail
            ),
            "minimum_global_BH_q_within_variant": min(
                float(row["global_BH_q_within_variant"])
                for row in variant_detail
            ),
            "global_metrics_BH_q_below_0.05": sum(
                float(row["global_BH_q_within_variant"]) < 0.05
                for row in variant_detail
            ),
            "minimum_reference_matched_empirical_upper_p": min(
                float(row["matched_empirical_upper_p"])
                for row in variant_detail
                if not math.isnan(float(row["matched_empirical_upper_p"]))
            ),
            "minimum_reference_matched_BH_q_within_variant": min(
                float(row["matched_BH_q_within_variant"])
                for row in variant_detail
                if row["matched_BH_q_within_variant"] != ""
            ),
            "reference_matched_metrics_BH_q_below_0.05": sum(
                row["matched_BH_q_within_variant"] != ""
                and float(row["matched_BH_q_within_variant"]) < 0.05
                for row in variant_detail
            ),
            "interpretation": (
                "sequence_property_screen_only; not evidence of biochemical "
                "function or adaptation"
            ),
        })
    summary_rows.sort(
        key=lambda row: (
            int(row["priority_tier"]),
            -float(row["mean_global_effect_percentile"]),
            str(row["variant"]),
        )
    )
    return detail_rows, summary_rows


def main() -> None:
    reference_sequence = human_sequence()
    natural = alternative_states(reference_sequence)
    candidates = candidate_variants()
    candidate_metadata = {
        str(row["variant"]): row for row in candidates
    }

    null_scores = [
        scored
        for variant in natural
        for scored in score_variant(reference_sequence, variant, "natural_null")
    ]
    candidate_scores = [
        scored
        for variant in candidates
        for scored in score_variant(reference_sequence, variant, "candidate")
    ]
    detail, candidate_summary = empirical_effect_summary(
        candidate_scores, null_scores, candidate_metadata
    )

    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT / "recurrent_natural_alternatives.tsv", natural)
    write_tsv(OUT / "natural_alternative_property_effects.tsv", null_scores)
    write_tsv(OUT / "candidate_property_effects.tsv", detail)
    write_tsv(OUT / "candidate_property_summary.tsv", candidate_summary)

    tier_one = [
        row for row in candidate_summary if int(row["priority_tier"]) == 1
    ]
    summary = {
        "human_reference": "DNMT3A1_912aa",
        "tail_interval": [TAIL_START, TAIL_END],
        "mammal_unique_sequences": len(read_fasta(ALIGNMENT)),
        "minimum_alternative_taxon_count": MIN_ALTERNATIVE_COUNT,
        "recurrent_natural_alternatives": len(natural),
        "candidate_variants": len(candidates),
        "window_sizes": list(WINDOW_SIZES),
        "metrics": list(METRICS),
        "tier_one_ranked_by_mean_global_effect_percentile": [
            {
                "variant": row["variant"],
                "mean_global_effect_percentile": row[
                    "mean_global_effect_percentile"
                ],
                "maximum_global_effect_percentile": row[
                    "maximum_global_effect_percentile"
                ],
            }
            for row in tier_one
        ],
        "inference_boundary": (
            "Property changes and empirical percentiles identify unusual "
            "sequence perturbations within this alignment-derived null. They "
            "do not validate phosphorylation, binding, methylation activity, "
            "brain function, or adaptive evolution."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
