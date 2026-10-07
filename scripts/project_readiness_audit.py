#!/usr/bin/env python3
"""Assemble a claim-level readiness audit for the DNMT3A project."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/08_synthesis"


def load_json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text())


def load_tsv(relative: str) -> list[dict[str, str]]:
    with (ROOT / relative).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    curated = load_json("results/01_curated/curation_summary.json")
    alignment = load_json("results/02_alignment/alignment_qc_summary.json")
    tree = load_json("results/03_tree/mammals/tree_qc_summary.json")
    meme = load_json(
        "results/04_selection/mammals/mammal_scope/full_meme/summary.json"
    )
    candidates = load_json(
        "results/04_selection/mammals/final_candidate_evidence.json"
    )
    phenotype = load_json(
        "results/05_brain_integration/"
        "candidate_brain_phenotype_coverage_summary.json"
    )
    atlas = load_json(
        "results/05_brain_integration/human_cell_atlas_import_summary.json"
    )
    dunnart = load_json(
        "results/06_isoform_evolution/"
        "dunnart_neocortex_stage_junctions_summary.json"
    )
    robustness = load_json(
        "results/06_isoform_evolution/"
        "dunnart_neocortex_stage_robustness_summary.json"
    )
    tss = load_json(
        "results/06_isoform_evolution/direct_tss_evidence_summary.json"
    )
    sra = load_json(
        "results/06_isoform_evolution/marsupial_sra_inventory_summary.json"
    )
    functional = load_json("results/07_functional_prioritization/summary.json")
    tail = load_json(
        "results/04_selection/mammals/regulatory_tail/summary.json"
    )
    idr = load_json(
        "results/07_functional_prioritization/in_silico_idr/summary.json"
    )
    esm = load_json(
        "results/07_functional_prioritization/in_silico_esm/"
        "model_sensitivity/summary.json"
    )
    dnmt3a1_promoter = load_json(
        "results/06_isoform_evolution/canonical_dnmt3a1_promoter/summary.json"
    )
    dnmt3a1_tissue = load_json(
        "results/06_isoform_evolution/canonical_dnmt3a1_promoter/"
        "canonical_fantom5_tissue_context_summary.json"
    )
    broad_isoform = load_json(
        "results/06_isoform_evolution/broad_isoform_phylogeny/summary.json"
    )
    human_promoter = load_json(
        "results/06_isoform_evolution/human_branch_promoter_screen/summary.json"
    )
    developmental = load_json(
        "results/07_isoform_function/GSE295720_developmental_methylation/summary.json"
    )
    mca_effects = {
        row["genotype"]: row
        for row in load_tsv("results/mch_global_summary/genotype_effect_summary.tsv")
    }
    structural_sensitivity = load_tsv(
        "structural_modeling/8QZM/md/matched_box/analysis/"
        "cross_preparation_comparison.tsv"
    )

    rows = [
        {
            "inference_layer": "ortholog_curation",
            "status": "complete",
            "supported_claim": (
                f"One exact protein/CDS-matched representative was selected "
                f"for {curated['selected_species']} species."
            ),
            "not_supported": "Every public model is biologically complete.",
            "next_requirement": "None for the current mammal coding analysis.",
        },
        {
            "inference_layer": "alignment_and_tree",
            "status": "complete_with_documented_caveats",
            "supported_claim": (
                f"Four alignment strategies and a {tree['tips']}-tip mammal "
                f"tree support site-level sensitivity analysis."
            ),
            "not_supported": (
                "Every N-terminal predicted residue is homologous or every "
                "branch is strongly resolved."
            ),
            "next_requirement": (
                "Independent transcript evidence for unresolved predicted "
                "N termini only if those sites are reconsidered."
            ),
        },
        {
            "inference_layer": "regulatory_tail_architecture",
            "status": "complete",
            "supported_claim": (
                f"The experimentally supported 164-219 engagement region has "
                f"{tail['engagement_region']['polymorphic_sites']} polymorphic "
                f"sites among {tail['engagement_region']['retained_sites']} "
                f"retained codons and is strongly conserved relative to the "
                f"upstream tail."
            ),
            "not_supported": (
                "Positive adaptation of the upstream tail or a functional "
                "effect of T12, G34, or S97."
            ),
            "next_requirement": (
                "A prespecified molecular assay or independent comparative "
                "phenotype."
            ),
        },
        {
            "inference_layer": "coding_selection",
            "status": "hypotheses_ready_not_discoveries",
            "supported_claim": (
                f"{meme['codon_sites_tested']} codons were tested across "
                f"{meme['sequences']} mammals; G34, T12, and S97 are the "
                f"qualified follow-up candidates."
            ),
            "not_supported": (
                "A corrected, alignment-independent positive-selection "
                "discovery at G34, T12, or S97."
            ),
            "next_requirement": (
                "Functional testing or a prespecified independent "
                "comparative hypothesis."
            ),
        },
        {
            "inference_layer": "functional_prioritization",
            "status": "in_silico_screen_complete_no_exceptional_candidates",
            "supported_claim": (
                f"{len(functional['tier_1_substitutions'])} recurrent "
                f"lineage-supported substitutions define a prespecified "
                f"panel; none has a corrected sequence-property result and "
                f"the ESM-2 natural-null rankings are strongly concordant "
                f"across model sizes (Spearman rho "
                f"{esm['natural_alternative_concordance']['spearman_rho']:.3f})."
            ),
            "not_supported": (
                "Any substitution has a demonstrated phenotype or an "
                "exceptional sequence-only perturbation."
            ),
            "next_requirement": (
                "Independent comparative phenotype data or a future "
                "molecular assay; additional sequence-only scoring is "
                "unlikely to resolve function."
            ),
        },
        {
            "inference_layer": "transposable_element_association",
            "status": "prespecified_blocked_on_source_workbook",
            "supported_claim": (
                "A young-TE association is biologically motivated and has a "
                "prespecified phylogenetic analysis plan."
            ),
            "not_supported": (
                "A TE association or TE-driven origin of T12, G34, or S97."
            ),
            "next_requirement": (
                "Obtain the original Osmanski et al. table S4 workbook and "
                "validate species-level matches before analysis."
            ),
        },
        {
            "inference_layer": "comparative_brain_association",
            "status": "blocked_on_new_species_data",
            "supported_claim": phenotype["conclusion"],
            "not_supported": (
                "Association between candidate residue state and brain mCH, "
                "developmental expression, or brain adaptation."
            ),
            "next_requirement": (
                "Comparable neuronal mCH or developmental-expression data "
                "from replicated non-reference residue lineages."
            ),
        },
        {
            "inference_layer": "human_cell_type_mch_context",
            "status": "context_complete",
            "supported_claim": (
                f"In {atlas['motor_cortex_cells']} adult human motor-cortex "
                f"nuclei, neuronal median mCH is "
                f"{atlas['motor_cortex_neuronal_mCH_median'] / atlas['motor_cortex_non_neuronal_mCH_median']:.2f}-fold "
                f"above non-neuronal median mCH, consistently in two donors."
            ),
            "not_supported": (
                "Cross-species evolution, developmental accumulation, or "
                "DNMT3A1-versus-DNMT3A2 expression."
            ),
            "next_requirement": (
                "Matched biological replicates in variant experiments; cells "
                "must not be treated as replicate organisms."
            ),
        },
        {
            "inference_layer": "marsupial_developmental_isoform_use",
            "status": "replicated_directional_evidence",
            "supported_claim": (
                f"Dunnart internal-junction signal is "
                f"{dunnart['descriptive_contrasts']['P12_over_P20']['internal_pairs_per_million_fold']:.2f}-fold "
                f"higher at P12 than P20 and the direction survives every "
                f"single-animal deletion."
            ),
            "not_supported": (
                "An internal-over-full-length isoform switch, exact capped "
                "initiation, or exclusion of sex/batch confounding."
            ),
            "next_requirement": (
                "Larger staged cohorts with sex/batch metadata and "
                "cap-enriched 5-prime assays."
            ),
        },
        {
            "inference_layer": "lineage_specific_promoter_isoform_evolution",
            "status": "therian_ancestry_supported_human_acceleration_not_detected",
            "supported_claim": (
                "Downstream-start-compatible products occur in "
                f"{broad_isoform['eutherian_orders_with_any_downstream_core_product']} "
                f"of {broad_isoform['eutherian_orders_total']} placental "
                "orders and two marsupial orders; neither tested promoter has "
                "an elevated fixed human-branch substitution fraction "
                "relative to local flanks."
            ),
            "not_supported": (
                "Human-specific promoter acceleration, adaptive brain "
                "evolution, or lineage-specific isoform losses."
            ),
            "next_requirement": (
                "Cap-enriched monotreme and marsupial TSS data, or an "
                "independent comparative promoter-usage phenotype."
            ),
        },
        {
            "inference_layer": "dnmt3a1_promoter_architecture",
            "status": "human_mouse_cluster_orthology_supported",
            "supported_claim": (
                "Human and mouse full-length DNMT3A1 TSS clusters have direct "
                "CAGE support; mapped annotated starts approach within "
                f"{dnmt3a1_promoter['human_mouse_synteny']['minimum_mapped_human_annotated_TSS_distance_to_mouse_TSS_bp']} "
                "bp and mapped CAGE CTSSs include an exact cross-species match."
            ),
            "not_supported": (
                "Basewise promoter conservation across Theria, brain-specific "
                "promoter activity, or conserved regulatory effect."
            ),
            "next_requirement": (
                "Cap-enriched TSS data and orthology-aware regulatory "
                "annotation in marsupials and monotremes."
            ),
        },
        {
            "inference_layer": "dnmt3a2_promoter_ancestry",
            "status": "architecture_supported_function_unresolved",
            "supported_claim": (
                f"Direct CAGE evidence supports {tss['models_with_direct_CAGE_peak_within_500bp']} "
                f"eutherian internal starts; marsupial transcript architecture "
                f"is independently supported in developing neocortex."
            ),
            "not_supported": (
                "Exact marsupial TSSs, promoter orthology across Theria, or "
                "conserved promoter function."
            ),
            "next_requirement": (
                f"Cap-enriched marsupial 5-prime data; the current inventory "
                f"contains {sra['exact_5prime_capable_runs']} exact-TSS-capable runs."
            ),
        },
        {
            "inference_layer": "developmental_isoform_function",
            "status": "stage_specific_effect_switch_supported",
            "supported_claim": (
                "In GSE295720 brain, Dnmt3a2 loss has the larger mean "
                f"methylation effect at E15.5 "
                f"({developmental['brain_stage_isoform_effects']['E15.5_Dnmt3a2-/-']:.5f}), "
                "whereas Dnmt3a1 loss has the larger effect at P21 "
                f"({developmental['brain_stage_isoform_effects']['PD21_Dnmt3a1-/-']:.5f})."
            ),
            "not_supported": (
                "Isoform exclusivity, a conserved cross-species handoff, or "
                "evolutionary causality."
            ),
            "next_requirement": (
                "Matched isoform-resolved developmental methylomes from an "
                "independent mammalian lineage."
            ),
        },
        {
            "inference_layer": "p21_neuronal_mca",
            "status": "large_replicate_consistent_effects_with_small_animal_n",
            "supported_claim": (
                "P21 neuronal corrected mCA relative to WT is "
                f"{float(mca_effects['Dnmt3a1_KO']['corrected_mCA_relative_vs_WT']):.4f} "
                "after Dnmt3a1 knockout and "
                f"{float(mca_effects['Dnmt3a1_delta_N']['corrected_mCA_relative_vs_WT']):.3f} "
                "after N-terminal deletion; Dnmt3a2 knockout retains the bulk signal."
            ),
            "not_supported": (
                "Universal postnatal dispensability of DNMT3A2 or assignment "
                "of the deletion phenotype to residues 164-219 alone."
            ),
            "next_requirement": (
                "Independent animal-level replication and finer N-terminal "
                "perturbations for causal localization."
            ),
        },
        {
            "inference_layer": "l188h_structural_modeling",
            "status": "preparation_sensitive_hypothesis",
            "supported_claim": (
                f"Two 500 ps WT-L188H preparations were compared across "
                f"{len(structural_sensitivity)} structural metrics; alternative "
                "local interface states were sampled."
            ),
            "not_supported": (
                "A reproducible distance or mobility direction, stronger "
                "binding, altered mCA, or brain adaptation."
            ),
            "next_requirement": (
                "Multiple longer independent replicas or a validated "
                "alchemical free-energy design."
            ),
        },
    ]

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "project_readiness.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "computational_discovery_phase": "complete",
        "candidate_ranking": candidates["candidate_order"],
        "coding_candidates_for_followup": candidates["candidates_for_followup"],
        "brain_association_ready": bool(
            phenotype["association_ready_combinations"]
        ),
        "tier_1_functional_substitutions": functional["tier_1_substitutions"],
        "in_silico_functional_screen": {
            "property_metrics_per_variant": (
                len(idr["metrics"]) * len(idr["window_sizes"])
            ),
            "tier_one_candidates_with_global_BH_q_below_0.05": 0,
            "esm_natural_null_spearman_rho": esm[
                "natural_alternative_concordance"
            ]["spearman_rho"],
            "interpretation": (
                "No tier-1 candidate is exceptional in the current "
                "sequence-property or ESM-2 screens."
            ),
        },
        "dnmt3a1_promoter_screen": {
            "representative_promoters": dnmt3a1_promoter[
                "representative_promoters"
            ],
            "canonical_tss_cluster_models": dnmt3a1_promoter[
                "canonical_tss_cluster_models"
            ],
            "human_mouse_upstream_shared_15mers": dnmt3a1_promoter[
                "human_mouse_upstream_shared_15mers"
            ],
            "human_brain_CNS_cluster_detection_fraction": dnmt3a1_tissue[
                "species"
            ][0]["cluster_detected_brain_samples"]["fraction"],
            "mouse_brain_CNS_cluster_detection_fraction": dnmt3a1_tissue[
                "species"
            ][1]["cluster_detected_brain_samples"]["fraction"],
            "interpretation": (
                "The human and mouse DNMT3A1 promoter regions are best "
                "represented as syntenically corresponding active TSS "
                "clusters, not one-to-one transcript-variant-1 starts."
            ),
        },
        "dense_promoter_isoform_evolution": {
            "mammal_species": broad_isoform["mammal_species"],
            "orders": broad_isoform["orders"],
            "eutherian_orders_with_any_downstream_core_product": (
                broad_isoform[
                    "eutherian_orders_with_any_downstream_core_product"
                ]
            ),
            "metatherian_orders_with_any_downstream_core_product": (
                broad_isoform[
                    "metatherian_orders_with_any_downstream_core_product"
                ]
            ),
            "human_branch_promoter_screen_conclusion": human_promoter[
                "conclusion"
            ],
            "interpretation": (
                "The leading model is therian ancestry with regulatory "
                "remodeling; the current data do not support a human-specific "
                "promoter innovation."
            ),
        },
        "readiness_rows": rows,
        "highest_information_next_actions": [
            (
                "Acquire comparable neuronal mCH or developmental-expression "
                "data in replicated non-reference candidate lineages."
            ),
            (
                "Run the prespecified phylogenetic young-TE association after "
                "obtaining the original Zoonomia table S4 workbook."
            ),
            (
                "Generate cap-enriched marsupial 5-prime data and an "
                "orthology-aware therian promoter analysis."
            ),
        ],
        "stop_rule": (
            "Additional post hoc selection tests on the same alignments should "
            "not change candidate status without new homology, functional, or "
            "phenotype evidence."
        ),
        "revised_manuscript_focus": (
            "Regional constraint and developmental isoform specialization "
            "in mammalian DNMT3A."
        ),
        "regulatory_tail_interpretation": tail["interpretation"],
        "robustness_note": robustness["interpretation"],
        "vertebrate_scope_retained_alignment_codons": alignment[
            "retained_codon_columns"
        ],
        "mammal_scope_selection_codons": meme["codon_sites_tested"],
    }
    (OUT / "project_readiness.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
