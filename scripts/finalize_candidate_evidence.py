#!/usr/bin/env python3
"""Assemble the final evidence matrix for current DNMT3A coding candidates."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/04_selection/mammals"
MASK = ROOT / "results/02_alignment/reliability_mask"
SITES = ["P5", "S6", "T12", "A16", "E18", "G34", "S97", "A114"]
DECISIONS = {
    "G34": (
        1,
        "qualified_primary_coding_hypothesis",
        "retain_for_functional_and_comparative_followup",
    ),
    "T12": (
        2,
        "qualified_secondary_coding_hypothesis",
        "retain_with_composition_and_branch_posterior_qualification",
    ),
    "S6": (
        4,
        "alignment_and_annotation_sensitive_lead",
        "exploratory_only",
    ),
    "S97": (
        3,
        "mammal_scope_alignment_robust_tertiary_hypothesis",
        "retain_for_lineage_and_comparative_followup",
    ),
    "A16": (
        5,
        "mask_sensitivity_emergent_without_MEME_support",
        "do_not_promote",
    ),
    "P5": (
        6,
        "artifact_dependent_primary_discovery",
        "deprioritize",
    ),
    "E18": (
        7,
        "artifact_dependent_primary_discovery",
        "deprioritize",
    ),
    "A114": (
        8,
        "unsupported",
        "stop_followup_absent_new_evidence",
    ),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def indexed(path: Path, key: str = "human_site") -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_tsv(path)}


def main() -> None:
    full_rows = read_tsv(RESULTS / "full_meme/all_sites.tsv")
    full = {
        f"{row['human_residue']}{row['human_DNMT3A1_aa']}": row
        for row in full_rows
    }
    mammal_scope_full = indexed(
        RESULTS / "mammal_scope/full_meme/all_sites.tsv"
    )
    fubar = indexed(MASK / "candidate_fubar_comparison.tsv")
    discovery_mask = indexed(
        RESULTS
        / "full_meme/consensus_masked/consensus_masked_summary.tsv"
    )
    region_mask_meme = indexed(MASK / "candidate_meme_comparison.tsv")
    lineage_summary = {
        row["human_site"]: row
        for row in json.loads(
            (
                RESULTS
                / "reliability_mask/lineages/summary.json"
            ).read_text()
        )["sites"]
    }
    cobalt = indexed(
        RESULTS / "cobalt/cobalt_candidate_sensitivity.tsv",
        key="human_site",
    )
    s97_selection = {
        row["alignment"]: row
        for row in read_tsv(
            RESULTS / "s97_homology_audit/s97_selection_sensitivity.tsv"
        )
    }
    s97_lineages = json.loads(
        (RESULTS / "s97_homology_audit/lineages/summary.json").read_text()
    )
    scope_audit = json.loads(
        (RESULTS / "mammal_scope_filter_audit/summary.json").read_text()
    )
    biological_context = indexed(
        RESULTS / "candidate_biological_context/retained_candidate_context.tsv"
    )
    phenotype_coverage = json.loads(
        (
            ROOT
            / "results/05_brain_integration"
            / "candidate_brain_phenotype_coverage_summary.json"
        ).read_text()
    )

    rows: list[dict[str, object]] = []
    for site in SITES:
        full_row = full.get(site, {
            "human_DNMT3A1_aa": "97",
            "domain": "N_terminal_regulatory",
            "p_value": "",
            "bh_q_value_901_sites": "",
        })
        fubar_row = fubar.get(site, {
            "primary_fubar_posterior": "",
            "masked_fubar_posterior": "",
            "mafft_fubar_posterior": "",
            "prank_fubar_posterior": s97_selection["PRANK"][
                "fubar_positive_posterior"
            ],
        })
        if site in discovery_mask:
            mask_test = discovery_mask[site]
            masked_meme_p: object = mask_test["consensus_masked_p"]
            masked_meme_adjusted: object = mask_test[
                "consensus_masked_holm_p_3_candidates"
            ]
            mask_test_type = "site_specific_Puma_Vombatus_consensus_mask"
        elif site in region_mask_meme:
            mask_test = region_mask_meme[site]
            masked_meme_p = mask_test["meme_p"]
            masked_meme_adjusted = mask_test["holm_p_3_mask_followups"]
            mask_test_type = "taxon_wide_N_terminal_reliability_mask"
        else:
            masked_meme_p = ""
            masked_meme_adjusted = ""
            mask_test_type = "not_run"

        lineage = lineage_summary.get(site)
        if site == "G34":
            annotation_effect = "Carlito_artifact_masked;signal_persists"
            branch_test = (
                "Node145_Rhinolophus_and_Node288_lorisiform_retained;"
                "both_targeted_aBSREL_nonsignificant"
            )
        elif site == "T12":
            annotation_effect = "Puma_and_Vombatus_artifacts_masked;signal_persists"
            branch_test = (
                "nine_change_coupled_branches_retained;"
                "raw_MEME_branch_posterior_saturates_475_branches;"
                "no_formal_branch_confirmation"
            )
        elif site in {"P5", "S6", "E18"}:
            annotation_effect = (
                "Puma_first_exon_nonhomology_and_unresolved_Vombatus_drive_signal"
            )
            branch_test = "not_warranted"
        elif site == "A16":
            annotation_effect = "Puma_and_Vombatus_masking_creates_FUBAR_threshold_crossing"
            branch_test = "not_warranted"
        elif site == "S97":
            annotation_effect = (
                "original_omission_caused_by_475_vertebrate_filter;"
                "100_percent_mammal_occupancy_and_four_alignment_agreement"
            )
            branch_test = (
                "eight_MEME_posterior_change_coupled_branches;"
                "descriptive_not_branch_wide_tests"
            )
        else:
            annotation_effect = "no_alignment_discordance_at_site"
            branch_test = "not_warranted"

        priority, evidence_class, action = DECISIONS[site]
        context = biological_context.get(site, {})
        mammal_scope_row = mammal_scope_full[site]
        rows.append({
            "priority_rank": priority,
            "human_site": site,
            "human_DNMT3A1_aa": full_row["human_DNMT3A1_aa"],
            "domain": full_row["domain"],
            "full_scan_meme_p": full_row["p_value"],
            "full_scan_meme_bh_q_901": full_row["bh_q_value_901_sites"],
            "mammal_scope_full_meme_p": mammal_scope_row["p_value"],
            "mammal_scope_full_meme_bh_q_909": mammal_scope_row[
                "bh_q_value_909_sites"
            ],
            "mammal_scope_full_meme_bonferroni_p_909": mammal_scope_row[
                "bonferroni_p_909_sites"
            ],
            "primary_fubar_posterior": fubar_row["primary_fubar_posterior"],
            "masked_fubar_posterior": fubar_row["masked_fubar_posterior"],
            "mafft_fubar_posterior": fubar_row["mafft_fubar_posterior"],
            "prank_fubar_posterior": fubar_row["prank_fubar_posterior"],
            "mammal_scope_macse_fubar_posterior": (
                s97_selection["MACSE"]["fubar_positive_posterior"]
                if site == "S97" else ""
            ),
            "mammal_scope_macse_meme_p": (
                s97_selection["MACSE"]["meme_p"] if site == "S97" else ""
            ),
            "cobalt_fubar_posterior": cobalt[site][
                "cobalt_fubar_positive_posterior"
            ],
            "cobalt_meme_p": cobalt[site]["cobalt_meme_p"],
            "cobalt_meme_holm_p_8_candidates": cobalt[site][
                "cobalt_meme_holm_p_8_candidates"
            ],
            "masked_meme_test_type": mask_test_type,
            "masked_meme_p": masked_meme_p,
            "masked_meme_adjusted_p": masked_meme_adjusted,
            "annotation_alignment_effect": annotation_effect,
            "masked_prioritized_branch_count": (
                len(lineage["prioritized_nonsynonymous_change_branches"])
                if lineage else 0
            ),
            "s97_prioritized_branch_count": (
                len(s97_lineages["prioritized_change_coupled_branches"])
                if site == "S97" else 0
            ),
            "branch_level_status": branch_test,
            "evidence_class": evidence_class,
            "recommended_action": action,
            "uniprot_disordered_region": context.get(
                "uniprot_disordered_region", ""
            ),
            "dnmt3a1_specific_absent_from_isoform2": context.get(
                "dnmt3a1_specific_absent_from_isoform2", ""
            ),
            "alphafold_plddt": context.get("alphafold_plddt", ""),
            "exact_uniprot_modified_residue": context.get(
                "exact_uniprot_modified_residue", ""
            ),
            "nearest_uniprot_modified_residue": context.get(
                "nearest_uniprot_modified_residue", ""
            ),
            "brain_phenotype_inference": "none",
        })
    rows.sort(key=lambda row: int(row["priority_rank"]))

    with (RESULTS / "final_candidate_evidence.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "candidate_order": [row["human_site"] for row in rows],
        "candidates_for_followup": ["G34", "T12", "S97"],
        "exploratory_only": ["S6"],
        "do_not_promote": ["A16"],
        "deprioritized_or_stopped": ["P5", "E18", "A114"],
        "four_alignment_mammal_scope_fubar_candidates": ["T12", "G34", "S97"],
        "mammal_scope_full_meme_bh_discoveries": ["P5", "E18", "S6"],
        "scope_filter_audit": scope_audit["conclusion"],
        "brain_phenotype_coverage": phenotype_coverage["conclusion"],
        "rows": rows,
        "statistical_boundary": (
            "The corrected 909-codon mammal-scope MACSE scan finds P5, E18, "
            "and S6 at BH q < 0.05, but all three fail conservative "
            "alignment/annotation sensitivity tests and are not promoted. "
            "G34 and T12 are nominal in the 909-codon scan (both BH q=0.0871) "
            "and retain FUBAR plus masked nominal MEME support. S97 is "
            "alignment-robust but its full-scan p=0.0224 becomes BH q=0.7206. "
            "No candidate has confirmed branch-wide selection or a "
            "demonstrated brain phenotype effect."
        ),
        "s97_status": (
            "S97 has 100% occupancy and identical residue assignments across "
            "mammal-scope MACSE, MAFFT, PRANK, and COBALT. Its omission from "
            "the original MACSE/MAFFT subsets was caused by filtering across "
            "475 vertebrates before mammal subsetting. It is promoted to a "
            "qualified tertiary hypothesis, but its full mammal-scope MEME "
            "result (p=0.0224; BH q=0.7206 across 909 codons) is not a "
            "corrected discovery."
        ),
    }
    (RESULTS / "final_candidate_evidence.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
