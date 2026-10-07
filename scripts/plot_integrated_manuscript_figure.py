#!/usr/bin/env python3
"""Create the main synthesis figure and a chromatin-context supplement."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = ROOT / "figures" / "manuscript"

COLORS = {
    "upstream": "#4C78A8", "engagement": "#E45756", "a1": "#D95F02",
    "a2": "#1B9E77", "wt": "#4D4D4D", "delta": "#7570B3",
    "single": "#2F5597", "joint": "#C44E52", "null": "#666666",
    "shared": "#8DA0CB",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_figure(fig, stem: str) -> None:
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}",
                    dpi=600 if extension == "png" else None, bbox_inches="tight")
    plt.close(fig)


def panel_constraint(ax, region_rows) -> None:
    names = ["Residues 1-163\nvariable tail", "Residues 164-219\nengagement region"]
    rows = [region_rows["upstream_disordered_tail"],
            region_rows["nucleosome_H2AK119ub_engagement"]]
    polymorphic = [100 * float(row["polymorphic_fraction"]) for row in rows]
    recurrent = [100 * float(row["recurrently_variable_fraction"]) for row in rows]
    x = np.arange(2)
    colors = [COLORS["upstream"], COLORS["engagement"]]
    ax.bar(x - 0.18, polymorphic, 0.36, color=colors, label="Polymorphic")
    ax.bar(x + 0.18, recurrent, 0.36, color=colors, alpha=0.45, hatch="//",
           label="Recurrently variable")
    ax.set_xticks(x, names)
    ax.set_ylabel("Mammalian sites (%)")
    ax.set_ylim(0, 100)
    ax.text(0.50, 0.54, "Polymorphic depletion\n$p$ = 3.26 × 10$^{-24}$",
            transform=ax.transAxes, ha="center", va="center", fontsize=8.5)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.set_title("A  Constraint changes within the DNMT3A1 tail", loc="left",
                 fontweight="bold", fontsize=10.5)


def panel_isoform_architecture(ax, phylogeny_summary) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    for y, label in ((0.78, "DNMT3A1"), (0.58, "DNMT3A2-like")):
        ax.text(0.02, y + 0.045, label, ha="left", va="center", fontweight="bold")
    ax.add_patch(Rectangle((0.20, 0.78), 0.22, 0.09,
                           color=COLORS["upstream"], alpha=0.9))
    ax.text(0.31, 0.825, "N-terminal tail", color="white", ha="center",
            va="center", fontsize=8)
    for y in (0.78, 0.58):
        for x, width, label in ((0.44, 0.13, "PWWP"), (0.59, 0.12, "ADD"),
                                (0.73, 0.23, "Catalytic")):
            ax.add_patch(Rectangle((x, y), width, 0.09,
                                   color=COLORS["shared"], alpha=0.95))
            ax.text(x + width / 2, y + 0.045, label, ha="center", va="center",
                    fontsize=8)
    ax.annotate("internal start", xy=(0.44, 0.68), xytext=(0.30, 0.68),
                arrowprops={"arrowstyle": "->", "lw": 1.0, "color": COLORS["a2"]},
                color=COLORS["a2"], ha="center", va="center", fontsize=8.5)
    placental_present = phylogeny_summary["eutherian_orders_with_any_downstream_core_product"]
    placental_total = phylogeny_summary["eutherian_orders_total"]
    marsupial_present = len(phylogeny_summary["metatherian_orders_with_any_downstream_core_product"])
    ax.text(0.04, 0.36, "Downstream-start-compatible products",
            fontweight="bold", fontsize=9)
    ax.text(0.08, 0.25,
            f"Placental mammals: {placental_present}/{placental_total} sampled orders",
            fontsize=8.5)
    ax.text(0.08, 0.16, f"Marsupials: {marsupial_present} separated orders", fontsize=8.5)
    ax.text(0.08, 0.07, "Monotreme non-detection is inconclusive", fontsize=8.2,
            color="#555555")
    ax.set_title("B  Internal-transcript architecture is broadly therian",
                 loc="left", fontweight="bold", fontsize=10.5)


def panel_developmental_handoff(ax, developmental_rows) -> None:
    rows = {(row["stage"], row["knockout"]): row for row in developmental_rows
            if row["tissue"] == "brain" and row["sex"] == "pooled"}
    stages = ["E15.5", "PD21"]
    x = np.arange(2)
    a1 = [float(rows[(stage, "Dnmt3a1-/-")]["mean_KO_minus_WT"]) for stage in stages]
    a2 = [float(rows[(stage, "Dnmt3a2-/-")]["mean_KO_minus_WT"]) for stage in stages]
    ax.bar(x - 0.19, a1, 0.38, color=COLORS["a1"], label="Dnmt3a1 KO")
    ax.bar(x + 0.19, a2, 0.38, color=COLORS["a2"], label="Dnmt3a2 KO")
    ax.axhline(0, color="#666666", lw=0.8)
    ax.set_xticks(x, ["E15.5", "P21"])
    ax.set_ylabel("Mean brain methylation\n(KO minus WT beta)")
    ax.set_ylim(-0.082, 0.012)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.text(0.19, a2[0] - 0.004, f"{a2[0]:.3f}", ha="center", va="top", fontsize=8)
    ax.text(0.81, a1[1] - 0.004, f"{a1[1]:.3f}", ha="center", va="top", fontsize=8)
    ax.set_title("C  Isoform effects switch across development", loc="left",
                 fontweight="bold", fontsize=10.5)


def panel_mca(ax, global_rows, replicate_rows) -> None:
    order = ["WT", "Dnmt3a1_KO", "Dnmt3a2_KO", "Dnmt3a1_delta_N"]
    labels = ["WT", "Dnmt3a1 KO", "Dnmt3a2 KO", "Dnmt3a1 ΔN"]
    colors = [COLORS["wt"], COLORS["a1"], COLORS["a2"], COLORS["delta"]]
    values = {row["genotype"]: float(row["corrected_mCA_relative_vs_WT"])
              for row in global_rows}
    wt_mean = float(next(row["corrected_mCA_mean"] for row in global_rows
                         if row["genotype"] == "WT"))
    bars = ax.bar(np.arange(4), [values[name] for name in order], color=colors, width=0.68)
    for position, genotype in enumerate(order):
        points = [float(row["autosomal_corrected_mCA"]) / wt_mean
                  for row in replicate_rows if row["genotype"] == genotype]
        jitter = np.linspace(-0.07, 0.07, len(points))
        ax.scatter(position + jitter, points, color="white", edgecolor="#222222",
                   s=28, zorder=3, lw=0.8)
    for bar, genotype in zip(bars, order):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.035,
                f"{values[genotype]:.3f}", ha="center", va="bottom", fontsize=8)
    ax.axhline(1, color="#777777", lw=0.8, ls="--")
    ax.set_xticks(np.arange(4), labels, rotation=18, ha="right")
    ax.set_ylabel("Corrected autosomal mCA relative to WT")
    ax.set_ylim(0, 1.3)
    ax.text(0.02, 0.97, "$n$ = 2 animals per genotype", transform=ax.transAxes,
            ha="left", va="top", fontsize=8)
    ax.set_title("D  DNMT3A1 dominates P21 neuronal mCA", loc="left",
                 fontweight="bold", fontsize=10.5)


def forest(ax, labels, single, joint, xlabel, title, null=None) -> None:
    y = np.arange(len(labels))[::-1]
    tick_y, tick_labels = list(y), list(labels)
    for position, values in enumerate((single, joint)):
        offset = 0.12 if position == 0 else -0.12
        color = COLORS["single"] if position == 0 else COLORS["joint"]
        marker = "o" if position == 0 else "s"
        estimates = np.array([value[0] for value in values])
        low = np.array([value[1] for value in values])
        high = np.array([value[2] for value in values])
        ax.errorbar(estimates, y + offset, xerr=[estimates - low, high - estimates],
                    fmt=marker, color=color, ecolor=color, capsize=2.5,
                    markersize=5, lw=1.3,
                    label="Single-mark" if position == 0 else "Three-mark joint")
    if null is not None:
        estimate, low, high, label = null
        extra_y = -1.1
        ax.errorbar([estimate], [extra_y], xerr=[[estimate - low], [high - estimate]],
                    fmt="D", color=COLORS["null"], capsize=2.5, markersize=5, lw=1.3)
        tick_y.append(extra_y)
        tick_labels.append(label)
        ax.set_ylim(extra_y - 0.45, len(labels) - 0.45)
    ax.axvline(0, color="#777777", lw=0.8, ls="--")
    ax.set_yticks(tick_y, tick_labels)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=10.5)
    ax.grid(axis="x", color="#DDDDDD", lw=0.6)
    ax.legend(frameon=False, fontsize=8, loc="best")


def chromatin_supplement(mca_single, mca_joint, occupancy) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    feature_map = [("h2ak119ub", "H2AK119ub"), ("h3k27me3", "H3K27me3"),
                   ("h3k4me3", "H3K4me3")]
    selected_single = {row["feature"]: row for row in mca_single
                       if row["contrast"] == "Dnmt3a1_delta_N-WT"
                       and row["sensitivity"] == "WT_mCA_ge_0.005"}
    selected_joint = next(row for row in mca_joint
                          if row["contrast"] == "Dnmt3a1_delta_N-WT"
                          and row["sensitivity"] == "WT_mCA_ge_0.005"
                          and row["region"] == "gene_body"
                          and row["model"] == "three_marks")
    single_values, joint_values = [], []
    for feature, _ in feature_map:
        field = f"{feature}_gene_body_enrichment"
        row = selected_single[field]
        single_values.append((float(row[field + "_standardized_beta"]),
                              float(row[field + "_beta_block_interval_low"]),
                              float(row[field + "_beta_block_interval_high"])))
        joint_values.append((float(selected_joint[field + "_standardized_beta"]),
                             float(selected_joint[field + "_beta_block_interval_low"]),
                             float(selected_joint[field + "_beta_block_interval_high"])))
    null_field = "dnmt3a_occupancy_change_gene_body"
    null_row = selected_single[null_field]
    null = (float(null_row[null_field + "_standardized_beta"]),
            float(null_row[null_field + "_beta_block_interval_low"]),
            float(null_row[null_field + "_beta_block_interval_high"]),
            "DNMT3A occupancy change")
    forest(axes[0], [label for _, label in feature_map], single_values, joint_values,
           "Standardized β for ΔN minus WT mCA",
           "A  Chromatin context and ΔN mCA loss", null)
    occupancy_index = {(row["region"], row["model"]): row for row in occupancy}
    occ_single, occ_joint = [], []
    three_mark = occupancy_index[("gene_body", "three_marks")]
    for feature, _ in feature_map:
        field = f"{feature}_gene_body_enrichment"
        row = occupancy_index[("gene_body", feature)]
        occ_single.append((float(row[field + "_standardized_beta"]),
                           float(row[field + "_beta_block_interval_low"]),
                           float(row[field + "_beta_block_interval_high"])))
        occ_joint.append((float(three_mark[field + "_standardized_beta"]),
                          float(three_mark[field + "_beta_block_interval_low"]),
                          float(three_mark[field + "_beta_block_interval_high"])))
    forest(axes[1], [label for _, label in feature_map], occ_single, occ_joint,
           "Standardized β for ΔN minus WT occupancy",
           "B  Chromatin context and DNMT3A redistribution")
    fig.suptitle("Chromatin-context sensitivity analyses", fontsize=13,
                 fontweight="bold")
    save_figure(fig, "Figure_S1_chromatin_context")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    region_rows = {row["region"]: row for row in read_tsv(
        RESULTS / "04_selection/mammals/regulatory_tail/regulatory_subregion_summary.tsv")}
    with (RESULTS / "06_isoform_evolution/broad_isoform_phylogeny/summary.json").open() as handle:
        phylogeny_summary = json.load(handle)
    developmental_rows = read_tsv(
        RESULTS / "07_isoform_function/GSE295720_developmental_methylation/developmental_contrasts.tsv")
    global_rows = read_tsv(RESULTS / "mch_global_summary/genotype_effect_summary.tsv")
    replicate_rows = read_tsv(RESULTS / "mch_global_summary/replicate_qc_and_methylation.tsv")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlepad": 7, "figure.dpi": 150})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)
    panel_constraint(axes[0, 0], region_rows)
    panel_isoform_architecture(axes[0, 1], phylogeny_summary)
    panel_developmental_handoff(axes[1, 0], developmental_rows)
    panel_mca(axes[1, 1], global_rows, replicate_rows)
    fig.suptitle("Regional constraint and developmental isoform specialization in mammalian DNMT3A",
                 fontsize=14, fontweight="bold")
    save_figure(fig, "Figure_integrated_DNMT3A1_tail")
    chromatin_supplement(
        read_tsv(RESULTS / "mch_gene_body_comparison/p21_polycomb_mca_single_models.tsv"),
        read_tsv(RESULTS / "mch_gene_body_comparison/p21_polycomb_mca_joint_models.tsv"),
        read_tsv(RESULTS / "mch_gene_body_comparison/p21_dnmt3a_occupancy_redistribution_models.tsv"),
    )


if __name__ == "__main__":
    main()
