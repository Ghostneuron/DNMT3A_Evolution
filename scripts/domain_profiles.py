#!/usr/bin/env python3
"""Summarize DNMT3A mammalian evolution by human DNMT3A1 domain."""

from __future__ import annotations

import csv
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dnmt3a_pipeline import CODON_TABLE, read_fasta  # noqa: E402


RESULTS = ROOT / "results/04_selection/mammals"
CROSSWALK = ROOT / "results/02_alignment/human_coordinate_crosswalk.tsv"
MAFFT_CROSSWALK = ROOT / "results/02_alignment/sensitivity/MAFFT_human_coordinate_crosswalk.tsv"
ALIGNMENT = ROOT / "results/02_alignment/subsets/DNMT3A_mammals_HyPhy_unique.fasta"
PRIMARY_FUBAR = RESULTS / "DNMT3A.FUBAR.json"
COMPOSITION_FUBAR = RESULTS / "DNMT3A.composition_pass.FUBAR.json"
MAFFT_FUBAR = RESULTS / "DNMT3A.mafft_sensitivity.FUBAR.json"
DOMAINS = ROOT / "config/human_domains.tsv"
POSTERIOR_THRESHOLD = 0.90


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def fubar_rows(path: Path) -> list[dict[str, float]]:
    data = json.loads(path.read_text())
    headers = [
        entry[0] if isinstance(entry, list) else str(entry)
        for entry in data["MLE"]["headers"]
    ]
    return [
        {"site": site, **{header: values[index] for index, header in enumerate(headers)}}
        for site, values in enumerate(data["MLE"]["content"]["0"], start=1)
    ]


def amino_acid(codon: str) -> str:
    return "-" if codon == "---" else CODON_TABLE.get(codon, "X")


def normalized_entropy(counts: Counter[str]) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = -sum(
        (count / total) * math.log(count / total)
        for count in counts.values() if count
    )
    return entropy / math.log(20)


def hypergeometric_upper_tail(
    population: int, successes: int, draws: int, observed: int
) -> float:
    denominator = math.comb(population, draws)
    upper = min(successes, draws)
    lower = max(observed, draws - (population - successes))
    return sum(
        math.comb(successes, value)
        * math.comb(population - successes, draws - value)
        for value in range(lower, upper + 1)
    ) / denominator


def main() -> None:
    crosswalk_rows = read_tsv(CROSSWALK)
    crosswalk = {
        int(row["filtered_codon_site_1based"]): row for row in crosswalk_rows
    }
    domain_order = [row["domain"] for row in read_tsv(DOMAINS)]
    records = read_fasta(ALIGNMENT)
    primary = {int(row["site"]): row for row in fubar_rows(PRIMARY_FUBAR)}
    composition = {int(row["site"]): row for row in fubar_rows(COMPOSITION_FUBAR)}
    if set(primary) != set(crosswalk) or set(composition) != set(crosswalk):
        raise RuntimeError("Primary FUBAR site coordinates do not match the crosswalk")

    mafft_map = {
        int(row["filtered_codon_site_1based"]): row for row in read_tsv(MAFFT_CROSSWALK)
    }
    mafft = {int(row["site"]): row for row in fubar_rows(MAFFT_FUBAR)}
    mafft_by_human = {
        int(mafft_map[site]["human_DNMT3A1_aa"]): row
        for site, row in mafft.items()
        if mafft_map[site]["human_DNMT3A1_aa"]
    }

    site_rows: list[dict[str, object]] = []
    for site, mapped in crosswalk.items():
        state_counts: Counter[str] = Counter()
        for record in records:
            codon = record.sequence[(site - 1) * 3:site * 3]
            state_counts[amino_acid(codon)] += 1
        non_gap = Counter({state: count for state, count in state_counts.items() if state != "-"})
        total_non_gap = sum(non_gap.values())
        ordered_counts = sorted(non_gap.values(), reverse=True)
        nonmajor_count = total_non_gap - (ordered_counts[0] if ordered_counts else 0)
        human_position = int(mapped["human_DNMT3A1_aa"])
        mafft_row = mafft_by_human.get(human_position)
        site_rows.append({
            "filtered_codon_site_1based": site,
            "human_DNMT3A1_aa": human_position,
            "human_residue": mapped["human_residue"],
            "domain": mapped["domain"],
            "full_panel_occupancy": float(mapped["occupancy"]),
            "mammal_unique_occupancy": total_non_gap / len(records),
            "amino_acid_states": len(non_gap),
            "nonmajor_residue_count": nonmajor_count,
            "polymorphic": len(non_gap) >= 2,
            "recurrently_variable": nonmajor_count >= 2,
            "normalized_amino_acid_entropy": normalized_entropy(non_gap),
            "primary_alpha": primary[site]["alpha"],
            "primary_beta": primary[site]["beta"],
            "primary_beta_minus_alpha": primary[site]["beta-alpha"],
            "primary_negative_posterior": primary[site]["Prob[alpha>beta]"],
            "primary_positive_posterior": primary[site]["Prob[alpha<beta]"],
            "composition_positive_posterior": composition[site]["Prob[alpha<beta]"],
            "mafft_positive_posterior": (
                mafft_row["Prob[alpha<beta]"] if mafft_row else ""
            ),
            "primary_candidate": primary[site]["Prob[alpha<beta]"] >= POSTERIOR_THRESHOLD,
            "composition_candidate": (
                composition[site]["Prob[alpha<beta]"] >= POSTERIOR_THRESHOLD
            ),
            "mafft_candidate": (
                bool(mafft_row)
                and mafft_row["Prob[alpha<beta]"] >= POSTERIOR_THRESHOLD
            ),
        })

    by_domain: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in site_rows:
        by_domain[str(row["domain"])].append(row)
    domain_rows: list[dict[str, object]] = []
    for domain in domain_order:
        rows = by_domain[domain]
        domain_rows.append({
            "domain": domain,
            "retained_sites": len(rows),
            "mean_full_panel_occupancy": statistics.fmean(
                float(row["full_panel_occupancy"]) for row in rows
            ),
            "mean_mammal_unique_occupancy": statistics.fmean(
                float(row["mammal_unique_occupancy"]) for row in rows
            ),
            "polymorphic_sites": sum(bool(row["polymorphic"]) for row in rows),
            "recurrently_variable_sites": sum(
                bool(row["recurrently_variable"]) for row in rows
            ),
            "mean_normalized_amino_acid_entropy": statistics.fmean(
                float(row["normalized_amino_acid_entropy"]) for row in rows
            ),
            "median_primary_beta_minus_alpha": statistics.median(
                float(row["primary_beta_minus_alpha"]) for row in rows
            ),
            "mean_primary_positive_posterior": statistics.fmean(
                float(row["primary_positive_posterior"]) for row in rows
            ),
            "primary_candidates_ge_0_90": sum(
                bool(row["primary_candidate"]) for row in rows
            ),
            "composition_candidates_ge_0_90": sum(
                bool(row["composition_candidate"]) for row in rows
            ),
            "mafft_candidates_ge_0_90": sum(
                bool(row["mafft_candidate"]) for row in rows
            ),
            "union_candidates_ge_0_90": sum(
                bool(row["primary_candidate"])
                or bool(row["composition_candidate"])
                or bool(row["mafft_candidate"])
                for row in rows
            ),
            "negative_selection_sites_posterior_ge_0_90": sum(
                float(row["primary_negative_posterior"]) >= POSTERIOR_THRESHOLD
                for row in rows
            ),
        })

    candidate_sets = {
        "primary_fubar_ge_0.90": {
            int(row["filtered_codon_site_1based"])
            for row in site_rows if bool(row["primary_candidate"])
        },
        "composition_pass_fubar_ge_0.90": {
            int(row["filtered_codon_site_1based"])
            for row in site_rows if bool(row["composition_candidate"])
        },
        "mafft_projection_fubar_ge_0.90": {
            int(row["filtered_codon_site_1based"])
            for row in site_rows if bool(row["mafft_candidate"])
        },
    }
    candidate_sets["union_three_fubar_scans"] = set().union(*candidate_sets.values())
    universe_filters = {
        "all_retained_sites": lambda row: True,
        "polymorphic_sites": lambda row: bool(row["polymorphic"]),
        "recurrently_variable_sites": lambda row: bool(row["recurrently_variable"]),
    }
    enrichment_rows: list[dict[str, object]] = []
    for candidate_name, candidate_sites in candidate_sets.items():
        for universe_name, include in universe_filters.items():
            universe = [row for row in site_rows if include(row)]
            universe_sites = {
                int(row["filtered_codon_site_1based"]) for row in universe
            }
            candidates = candidate_sites & universe_sites
            n_terminal = {
                int(row["filtered_codon_site_1based"])
                for row in universe if row["domain"] == "N_terminal_regulatory"
            }
            observed = len(candidates & n_terminal)
            p_value = hypergeometric_upper_tail(
                len(universe), len(n_terminal), len(candidates), observed
            )
            enrichment_rows.append({
                "candidate_definition": candidate_name,
                "site_universe": universe_name,
                "universe_sites": len(universe),
                "n_terminal_sites_in_universe": len(n_terminal),
                "candidate_sites": len(candidates),
                "n_terminal_candidates": observed,
                "one_sided_hypergeometric_p": p_value,
                "candidate_human_positions": ";".join(
                    str(crosswalk[site]["human_DNMT3A1_aa"])
                    for site in sorted(candidates)
                ),
            })

    n_terminal_sites = [
        row for row in site_rows if row["domain"] == "N_terminal_regulatory"
    ]
    other_sites = [
        row for row in site_rows if row["domain"] != "N_terminal_regulatory"
    ]
    variability_rows: list[dict[str, object]] = []
    for label, field in (
        ("polymorphic", "polymorphic"),
        ("recurrently_variable", "recurrently_variable"),
    ):
        total_variable = sum(bool(row[field]) for row in site_rows)
        observed_n_terminal = sum(bool(row[field]) for row in n_terminal_sites)
        other_variable = sum(bool(row[field]) for row in other_sites)
        n_terminal_nonvariable = len(n_terminal_sites) - observed_n_terminal
        other_nonvariable = len(other_sites) - other_variable
        odds_ratio = (
            observed_n_terminal * other_nonvariable
            / (n_terminal_nonvariable * other_variable)
            if n_terminal_nonvariable and other_variable else float("inf")
        )
        variability_rows.append({
            "variability_definition": label,
            "all_retained_sites": len(site_rows),
            "all_variable_sites": total_variable,
            "n_terminal_sites": len(n_terminal_sites),
            "n_terminal_variable_sites": observed_n_terminal,
            "n_terminal_variable_fraction": observed_n_terminal / len(n_terminal_sites),
            "other_sites": len(other_sites),
            "other_variable_sites": other_variable,
            "other_variable_fraction": other_variable / len(other_sites),
            "n_terminal_vs_other_odds_ratio": odds_ratio,
            "one_sided_hypergeometric_p": hypergeometric_upper_tail(
                len(site_rows), total_variable, len(n_terminal_sites),
                observed_n_terminal,
            ),
        })

    top_rows = sorted(
        site_rows,
        key=lambda row: float(row["primary_positive_posterior"]),
        reverse=True,
    )[:25]
    site_fields = list(site_rows[0])
    domain_fields = list(domain_rows[0])
    enrichment_fields = list(enrichment_rows[0])
    write_tsv(RESULTS / "site_domain_metrics.tsv", site_rows, site_fields)
    write_tsv(RESULTS / "domain_evolution_summary.tsv", domain_rows, domain_fields)
    write_tsv(
        RESULTS / "candidate_enrichment_tests.tsv",
        enrichment_rows,
        enrichment_fields,
    )
    write_tsv(
        RESULTS / "n_terminal_variability_tests.tsv",
        variability_rows,
        list(variability_rows[0]),
    )
    write_tsv(RESULTS / "top_primary_fubar_sites.tsv", top_rows, site_fields)
    summary = {
        "mammal_unique_sequences": len(records),
        "retained_sites": len(site_rows),
        "posterior_candidate_threshold": POSTERIOR_THRESHOLD,
        "domain_count": len(domain_rows),
        "candidate_sets": {
            name: sorted(values) for name, values in candidate_sets.items()
        },
    }
    (RESULTS / "domain_profile_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
