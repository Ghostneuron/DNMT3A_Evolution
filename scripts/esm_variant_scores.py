#!/usr/bin/env python3
"""Score DNMT3A1-tail substitutions with ESM-2 masked marginals."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import esm
import torch

from idr_variant_effects import human_sequence


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/07_functional_prioritization/in_silico_idr"
OUT = ROOT / "results/07_functional_prioritization/in_silico_esm"
MODEL_LOADERS = {
    "esm2_t12_35M_UR50D": esm.pretrained.esm2_t12_35M_UR50D,
    "esm2_t30_150M_UR50D": esm.pretrained.esm2_t30_150M_UR50D,
}


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


def choose_device(requested: str) -> torch.device:
    if requested != "auto":
        return torch.device(requested)
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def masked_position_log_probabilities(
    sequence: str,
    positions: list[int],
    model: torch.nn.Module,
    alphabet,
    device: torch.device,
    batch_size: int,
) -> dict[int, torch.Tensor]:
    converter = alphabet.get_batch_converter()
    mask = alphabet.mask_idx
    results: dict[int, torch.Tensor] = {}
    model.eval()
    for offset in range(0, len(positions), batch_size):
        batch_positions = positions[offset:offset + batch_size]
        entries = []
        for position in batch_positions:
            masked = sequence[:position - 1] + "<mask>" + sequence[position:]
            entries.append((str(position), masked))
        _, _, tokens = converter(entries)
        tokens = tokens.to(device)
        with torch.no_grad():
            logits = model(tokens, repr_layers=[], return_contacts=False)["logits"]
            log_probs = torch.log_softmax(logits, dim=-1)
        for index, position in enumerate(batch_positions):
            # ESM prepends a beginning-of-sequence token, so one-based protein
            # coordinates equal token coordinates.
            results[position] = log_probs[index, position].detach().cpu()
    return results


def score_rows(
    rows: list[dict[str, str]],
    log_probs: dict[int, torch.Tensor],
    alphabet,
    model_name: str,
    source: str,
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for row in rows:
        position = int(row["position"])
        reference = row["reference_residue"]
        alternative = row["alternative_residue"]
        reference_index = alphabet.get_idx(reference)
        alternative_index = alphabet.get_idx(alternative)
        reference_logp = float(log_probs[position][reference_index])
        alternative_logp = float(log_probs[position][alternative_index])
        llr = alternative_logp - reference_logp
        output.append({
            "model": model_name,
            "source": source,
            "variant": row["variant"],
            "position": position,
            "reference_residue": reference,
            "alternative_residue": alternative,
            "reference_log_probability": reference_logp,
            "alternative_log_probability": alternative_logp,
            "alternative_vs_reference_log_likelihood_ratio": llr,
            "model_prefers_alternative": llr > 0,
            "priority_tier": row.get("priority_tier", ""),
            "alternative_taxon_count": row.get("alternative_taxon_count", ""),
        })
    return output


def empirical_candidate_comparison(
    candidates: list[dict[str, object]],
    null: list[dict[str, object]],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for candidate in candidates:
        observed = float(
            candidate["alternative_vs_reference_log_likelihood_ratio"]
        )
        matched = [
            row for row in null
            if row["reference_residue"] == candidate["reference_residue"]
            and row["variant"] != candidate["variant"]
        ]

        def compare(pool: list[dict[str, object]]) -> tuple[float, float, float]:
            values = [
                float(row["alternative_vs_reference_log_likelihood_ratio"])
                for row in pool
            ]
            lower_p = (1 + sum(value <= observed for value in values)) / (
                len(values) + 1
            )
            upper_p = (1 + sum(value >= observed for value in values)) / (
                len(values) + 1
            )
            percentile = sum(value < observed for value in values) / len(values)
            return percentile, lower_p, upper_p

        global_percentile, global_lower, global_upper = compare(null)
        if matched:
            matched_percentile, matched_lower, matched_upper = compare(matched)
        else:
            matched_percentile = matched_lower = matched_upper = float("nan")
        output.append({
            **candidate,
            "global_null_variants": len(null),
            "global_LLR_percentile": global_percentile,
            "global_empirical_lower_p": global_lower,
            "global_empirical_upper_p": global_upper,
            "reference_matched_null_variants": len(matched),
            "reference_matched_LLR_percentile": matched_percentile,
            "reference_matched_empirical_lower_p": matched_lower,
            "reference_matched_empirical_upper_p": matched_upper,
            "interpretation": (
                "ESM masked-marginal sequence compatibility; "
                "not a biochemical or adaptive-effect measurement"
            ),
        })
    output.sort(
        key=lambda row: (
            int(row["priority_tier"]),
            float(row["alternative_vs_reference_log_likelihood_ratio"]),
        )
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=sorted(MODEL_LOADERS), required=True)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--batch-size", type=int, default=4)
    args = parser.parse_args()

    natural = read_tsv(INPUT / "recurrent_natural_alternatives.tsv")
    candidate_property = read_tsv(INPUT / "candidate_property_summary.tsv")
    positions = sorted({
        int(row["position"]) for row in natural + candidate_property
    })
    sequence = human_sequence()
    device = choose_device(args.device)
    model, alphabet = MODEL_LOADERS[args.model]()
    model = model.to(device)
    log_probs = masked_position_log_probabilities(
        sequence,
        positions,
        model,
        alphabet,
        device,
        args.batch_size,
    )

    null_scores = score_rows(
        natural, log_probs, alphabet, args.model, "natural_null"
    )
    candidate_scores = score_rows(
        candidate_property, log_probs, alphabet, args.model, "candidate"
    )
    comparison = empirical_candidate_comparison(candidate_scores, null_scores)

    model_out = OUT / args.model
    model_out.mkdir(parents=True, exist_ok=True)
    write_tsv(model_out / "natural_alternative_scores.tsv", null_scores)
    write_tsv(model_out / "candidate_scores.tsv", comparison)
    summary = {
        "model": args.model,
        "device": str(device),
        "scoring_method": (
            "wild-type human DNMT3A1 full-sequence masked-marginal "
            "log-likelihood ratio"
        ),
        "protein_length": len(sequence),
        "masked_positions": len(positions),
        "natural_null_variants": len(null_scores),
        "candidate_variants": len(candidate_scores),
        "tier_one_candidates": [
            {
                "variant": row["variant"],
                "log_likelihood_ratio": row[
                    "alternative_vs_reference_log_likelihood_ratio"
                ],
                "reference_matched_percentile": row[
                    "reference_matched_LLR_percentile"
                ],
                "reference_matched_lower_p": row[
                    "reference_matched_empirical_lower_p"
                ],
            }
            for row in comparison
            if int(row["priority_tier"]) == 1
        ],
        "inference_boundary": (
            "A masked-marginal score measures compatibility with patterns "
            "learned from protein sequences. It does not measure DNMT3A "
            "activity, phosphorylation, chromatin binding, brain function, "
            "or natural selection."
        ),
    }
    (model_out / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
