#!/usr/bin/env python3
"""Generate the submission-oriented DNMT3A main-figure package."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures/manuscript"

BLUE = "#3274A1"
ORANGE = "#E1812C"
TEAL = "#3A923A"
RED = "#C44E52"
PURPLE = "#8172B3"
GRAY = "#777777"
LIGHT = "#E7E7E7"
DARK = "#262626"


def read_tsv(path: str) -> list[dict[str, str]]:
    with (ROOT / path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{suffix}", dpi=350, bbox_inches="tight")
    plt.close(fig)


def setup() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.dpi": 150,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel(ax: plt.Axes, label: str, title: str) -> None:
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top")
    ax.set_title(title, loc="left", fontweight="bold", pad=9)


def clean(ax: plt.Axes, grid_axis: str | None = None) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color="#DDDDDD", linewidth=0.6, zorder=0)


def figure1() -> None:
    regions = read_tsv(
        "results/04_selection/mammals/regulatory_tail/"
        "regulatory_subregion_summary.tsv"
    )
    candidates = read_tsv(
        "results/04_selection/mammals/final_candidate_evidence.tsv"
    )
    fig = plt.figure(figsize=(12.2, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.82, 1.18])
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    panel(ax_a, "A", "DNMT3A1 domain architecture and regional variation")
    ax_a.set_xlim(0, 912)
    ax_a.set_ylim(-0.62, 1.35)
    ax_a.axis("off")
    domain_defs = [
        (1, 163, "Variable tail", ORANGE),
        (164, 219, "Engage.", RED),
        (220, 277, "pre-\nPWWP", "#D9A441"),
        (278, 427, "PWWP", BLUE),
        (428, 475, "linker", LIGHT),
        (476, 614, "ADD", PURPLE),
        (615, 633, "", LIGHT),
        (634, 912, "Catalytic", TEAL),
    ]
    for start, end, label, color in domain_defs:
        ax_a.add_patch(Rectangle((start, 0.45), end - start + 1, 0.33,
                                 facecolor=color, edgecolor="white", lw=0.8))
        if label:
            ax_a.text((start + end) / 2, 0.615, label, ha="center", va="center",
                      fontsize=8, color="white" if color != LIGHT else DARK,
                      fontweight="bold")
    ax_a.text(1, 0.91, "N terminus", ha="left", va="bottom", fontsize=8)
    ax_a.text(912, 0.91, "C terminus", ha="right", va="bottom", fontsize=8)
    by_name = {r["region"]: r for r in regions}
    for name, y, color in [
        ("upstream_disordered_tail", 0.05, ORANGE),
        ("nucleosome_H2AK119ub_engagement", -0.22, RED),
    ]:
        row = by_name[name]
        start, end = float(row["start_aa"]), float(row["end_aa"])
        frac = 100 * float(row["polymorphic_fraction"])
        ax_a.plot([start, end], [y, y], color=color, lw=5, solid_capstyle="butt")
        ax_a.text((start + end) / 2, y - 0.08, f"{frac:.1f}% variable",
                  ha="center", va="top", fontsize=8, color=color)
    for site in (12, 34, 97):
        ax_a.plot([site, site], [0.80, 1.05], color=DARK, lw=0.8)
        ax_a.text(site, 1.10, {12: "T12", 34: "G34", 97: "S97"}[site],
                  ha="center", fontsize=8)

    names = [r["region"] for r in regions]
    labels = ["1–163", "164–219", "220–277", "PWWP", "PWWP–ADD",
              "ADD", "ADD–MTase", "Catalytic"]
    y = np.arange(len(regions))
    poly = np.array([float(r["polymorphic_fraction"]) for r in regions])
    recur = np.array([float(r["recurrently_variable_fraction"]) for r in regions])
    height = 0.35
    ax_b.barh(y + height / 2, poly, height, color=BLUE, label="Variable across species")
    ax_b.barh(y - height / 2, recur, height, color=ORANGE, label="Recurrent")
    ax_b.set_yticks(y, labels)
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 1)
    ax_b.set_xlabel("Fraction of retained sites")
    ax_b.legend(frameon=False, loc="lower right")
    panel(ax_b, "B", "Variation changes sharply within the N terminus")
    clean(ax_b, "x")
    ax_b.text(0.98, 0.97, "p = 3.26 × 10⁻²⁴",
              transform=ax_b.transAxes, ha="right", va="top", fontsize=8,
              color=RED, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))

    tier = [r for r in candidates if r["human_site"] in {"T12", "G34", "S97"}]
    tier.sort(key=lambda r: ["T12", "G34", "S97"].index(r["human_site"]))
    cols = [
        ("MACSE", "primary_fubar_posterior", "mammal_scope_macse_fubar_posterior"),
        ("Masked", "masked_fubar_posterior", None),
        ("MAFFT", "mafft_fubar_posterior", None),
        ("PRANK", "prank_fubar_posterior", None),
        ("COBALT", "cobalt_fubar_posterior", None),
    ]
    matrix = np.full((3, len(cols)), np.nan)
    for i, row in enumerate(tier):
        for j, (_, primary, fallback) in enumerate(cols):
            value = row.get(primary, "")
            if not value and fallback:
                value = row.get(fallback, "")
            if value:
                matrix[i, j] = float(value)
    cmap = LinearSegmentedColormap.from_list("fubar", ["#F2F2F2", "#A9CCE3", BLUE])
    cmap.set_bad("white")
    im = ax_c.imshow(matrix, aspect="auto", vmin=0.8, vmax=1.0, cmap=cmap)
    ax_c.set_xticks(range(len(cols)), [c[0] for c in cols])
    ax_c.set_yticks(range(3), [r["human_site"] for r in tier])
    for i in range(3):
        for j in range(len(cols)):
            txt = "NA" if np.isnan(matrix[i, j]) else f"{matrix[i, j]:.3f}"
            ax_c.text(j, i, txt, ha="center", va="center", fontsize=8,
                      color="white" if not np.isnan(matrix[i, j]) and matrix[i, j] > 0.96 else DARK)
    ax_c.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
    ax_c.set_yticks(np.arange(-0.5, 3, 1), minor=True)
    ax_c.grid(which="minor", color="white", linewidth=1.5)
    ax_c.tick_params(which="minor", bottom=False, left=False)
    panel(ax_c, "C", "Candidate FUBAR support is alignment-aware")
    cbar = fig.colorbar(im, ax=ax_c, fraction=0.046, pad=0.03)
    cbar.set_label("Posterior probability")
    for i, row in enumerate(tier):
        q = float(row["mammal_scope_full_meme_bh_q_909"])
        ax_c.text(len(cols) - 0.5, i + 0.40, f"full-scan MEME q = {q:.3f}",
                  ha="right", va="bottom", fontsize=7, color=GRAY)

    fig.suptitle("Figure 1 | Regional evolutionary constraint in mammalian DNMT3A",
                 x=0.01, ha="left", fontsize=14, fontweight="bold")
    save(fig, "Figure_1_regional_constraint")


def figure2() -> None:
    clades = read_tsv(
        "results/06_isoform_evolution/broad_isoform_phylogeny/"
        "isoform_clade_distribution.tsv"
    )
    depth = read_tsv(
        "results/06_isoform_evolution/broad_isoform_phylogeny/"
        "annotation_depth_sensitivity.tsv"
    )
    promoter = read_tsv(
        "results/06_isoform_evolution/human_branch_promoter_screen/"
        "promoter_enrichment_summary.tsv"
    )
    dunnart = read_tsv(
        "results/06_isoform_evolution/dunnart_neocortex_stage_junctions.tsv"
    )
    fig = plt.figure(figsize=(12.5, 10.5), constrained_layout=True)
    gs = fig.add_gridspec(3, 3, height_ratios=[0.7, 1.05, 1.05],
                          width_ratios=[1.25, 1, 1])
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1:, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[1, 2])
    ax_e = fig.add_subplot(gs[2, 1:])

    panel(ax_a, "A", "Full-length and internal DNMT3A transcript architectures")
    ax_a.set_xlim(0, 100)
    ax_a.set_ylim(-0.2, 3.15)
    ax_a.axis("off")
    for y, name, start, evidence in [
        (2.25, "DNMT3A1", 3, "human/mouse CAGE-supported promoter clusters"),
        (1.25, "DNMT3A2-like", 27, "internal transcript enters shared coding region"),
    ]:
        ax_a.text(0, y + 0.18, name, ha="left", va="center", fontweight="bold")
        if start == 3:
            blocks = [(3, 38, "N-terminal region", "#D9A441"),
                      (41, 16, "PWWP", BLUE), (59, 13, "ADD", PURPLE),
                      (75, 20, "Catalytic", TEAL)]
        else:
            blocks = [(20, 6, "leader", ORANGE), (27, 14, "PWWP", BLUE),
                      (44, 13, "ADD", PURPLE), (60, 35, "Catalytic", TEAL)]
        for block_start, width, block_label, color in blocks:
            ax_a.add_patch(Rectangle((block_start, y), width, 0.35,
                                     facecolor=color, edgecolor="white"))
            ax_a.text(block_start + width / 2, y + 0.18, block_label,
                      ha="center", va="center", fontsize=7.5, color="white",
                      fontweight="bold")
        ax_a.text(98, y - 0.13, evidence, ha="right", va="top", fontsize=8,
                  color=GRAY)
    ax_a.text(50, 0.42,
              "Downstream-start-compatible products: 13/18 placental orders + two separated marsupial orders",
              ha="center", va="center", color=DARK, fontsize=9,
              bbox=dict(boxstyle="round,pad=0.35", fc="#F6F6F6", ec="#BBBBBB"))

    order_rows = [r for r in clades if r["group_type"] == "order"]
    marsupial = {"Dasyuromorphia", "Didelphimorphia", "Diprotodontia", "Peramelemorphia"}
    monotreme = {"Monotremata"}
    placental = [r for r in order_rows if r["group"] not in marsupial | monotreme]
    meta = [r for r in order_rows if r["group"] in marsupial]
    proto = [r for r in order_rows if r["group"] in monotreme]
    ordered = placental + meta + proto
    mat = np.array([[int(r["strict_A2_like_species"]) > 0,
                     int(r["any_downstream_core_species"]) > 0] for r in ordered], dtype=float)
    ax_b.imshow(mat, cmap=LinearSegmentedColormap.from_list("binary", ["#F0F0F0", BLUE]),
                vmin=0, vmax=1, aspect="auto")
    ax_b.set_xticks([0, 1], ["Strict A2-like", "Any downstream\ncore"])
    ax_b.set_yticks(range(len(ordered)), [r["group"] for r in ordered])
    for i, r in enumerate(ordered):
        for j in range(2):
            ax_b.text(j, i, "present" if mat[i, j] else "—", ha="center", va="center",
                      fontsize=6.8, color="white" if mat[i, j] else GRAY)
    ax_b.axhline(len(placental) - 0.5, color=RED, lw=1.2)
    ax_b.axhline(len(placental) + len(meta) - 0.5, color=RED, lw=1.2)
    ax_b.text(1.62, (len(placental) - 1) / 2, "Placental", rotation=90,
              va="center", ha="center", fontsize=8)
    ax_b.text(1.62, len(placental) + (len(meta) - 1) / 2, "Marsupial", rotation=90,
              va="center", ha="center", fontsize=8)
    ax_b.text(1.62, len(ordered) - 1, "Monotreme", rotation=90,
              va="center", ha="center", fontsize=8)
    panel(ax_b, "B", "Order-level annotation evidence")
    ax_b.tick_params(length=0)

    panel(ax_c, "C", "Developing dunnart neocortex")
    for i, stage in enumerate(["P12", "P20"]):
        vals = [float(r["internal_067_first_junction_pairs_per_million"])
                for r in dunnart if r["stage"] == stage]
        offsets = np.linspace(-0.08, 0.08, len(vals))
        ax_c.scatter(np.full(len(vals), i) + offsets, vals, s=42,
                     color=[ORANGE, BLUE][i], edgecolor="white", linewidth=0.7, zorder=3)
        ax_c.plot([i - 0.18, i + 0.18], [np.mean(vals)] * 2, color=DARK, lw=2)
    ax_c.set_xticks([0, 1], ["P12", "P20"])
    ax_c.set_ylabel("Internal-junction pairs per million")
    ax_c.text(0.98, 0.96, "3.40-fold pooled\nexact permutation p = 0.05",
              transform=ax_c.transAxes, ha="right", va="top", fontsize=8)
    clean(ax_c, "y")

    panel(ax_d, "D", "Detection depends on annotation depth")
    x = np.arange(2)
    strict = [float(r["strict_A2_like_detection_fraction"]) for r in depth]
    broad = [float(r["any_downstream_core_detection_fraction"]) for r in depth]
    ax_d.bar(x - 0.18, strict, 0.36, color=ORANGE, label="Strict A2-like")
    ax_d.bar(x + 0.18, broad, 0.36, color=BLUE, label="Any downstream core")
    ax_d.set_xticks(x, ["<8 products", "≥8 products"])
    ax_d.set_ylim(0, 1)
    ax_d.set_ylabel("Species detection fraction")
    ax_d.legend(frameon=False, loc="upper left")
    clean(ax_d, "y")

    panel(ax_e, "E", "No promoter enrichment for candidate human-lineage differences")
    x = np.arange(2)
    pvals = [100 * float(r["promoter_substitution_fraction"]) for r in promoter]
    fvals = [100 * float(r["flank_substitution_fraction"]) for r in promoter]
    ax_e.bar(x - 0.18, pvals, 0.36, color=PURPLE, label="Promoter")
    ax_e.bar(x + 0.18, fvals, 0.36, color=GRAY, label="Local flank")
    ax_e.set_xticks(x, ["DNMT3A1 promoter cluster", "DNMT3A2 internal promoter"])
    ax_e.set_ylabel("Human-reference differences (%)")
    ax_e.legend(frameon=False, loc="upper right")
    for i, row in enumerate(promoter):
        ax_e.text(i, max(pvals[i], fvals[i]) + 0.025,
                  f"one-sided p = {float(row['one_sided_enrichment_p']):.3f}",
                  ha="center", fontsize=8)
    ax_e.set_ylim(0, max(fvals) * 1.35)
    clean(ax_e, "y")

    fig.suptitle("Figure 2 | Mammalian evolution of DNMT3A transcript architecture",
                 x=0.01, ha="left", fontsize=14, fontweight="bold")
    save(fig, "Figure_2_transcript_architecture")


def figure3() -> None:
    pooled = read_tsv(
        "results/07_isoform_function/GSE295720_developmental_methylation/"
        "developmental_contrasts.tsv"
    )
    sex = read_tsv(
        "results/07_isoform_function/GSE295720_developmental_methylation/"
        "developmental_contrasts_by_sex.tsv"
    )
    cpg = read_tsv(
        "results/07_isoform_function/GSE164265_neuron_nuclei_methylation/"
        "global_effects.tsv"
    )
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.4), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.flat

    brain = [r for r in pooled if r["tissue"] == "brain"]
    x = np.arange(2)
    for offset, knockout, color, label in [(-0.18, "Dnmt3a1-/-", BLUE, r"$\mathit{Dnmt3a1}$ KO"),
                                           (0.18, "Dnmt3a2-/-", ORANGE, r"$\mathit{Dnmt3a2}$ KO")]:
        vals = [float(next(r["mean_KO_minus_WT"] for r in brain
                           if r["stage"] == stage and r["knockout"] == knockout))
                for stage in ["E15.5", "PD21"]]
        ax_a.bar(x + offset, vals, 0.36, color=color, label=label)
    ax_a.axhline(0, color=DARK, lw=0.8)
    ax_a.set_xticks(x, ["E15.5", "P21"])
    ax_a.set_ylabel("Mean methylation difference (KO − WT)")
    ax_a.legend(frameon=False)
    panel(ax_a, "A", "Isoform effects reverse across brain development")
    clean(ax_a, "y")

    positions, values, colors, labels = [], [], [], []
    idx = 0
    for stage in ["E15.5", "PD21"]:
        for tissue in ["brain", "liver"]:
            for knockout, color in [("Dnmt3a1-/-", BLUE), ("Dnmt3a2-/-", ORANGE)]:
                row = next(r for r in pooled if r["stage"] == stage and
                           r["tissue"] == tissue and r["knockout"] == knockout)
                positions.append(idx)
                values.append(float(row["mean_KO_minus_WT"]))
                colors.append(color)
                labels.append(f"{stage.replace('PD','P')}\n{tissue}")
                idx += 0.34
            idx += 0.36
        idx += 0.4
    ax_b.bar(positions, values, 0.28, color=colors)
    centers = [(positions[i] + positions[i + 1]) / 2 for i in range(0, len(positions), 2)]
    ax_b.set_xticks(centers, [labels[i] for i in range(0, len(labels), 2)])
    ax_b.axhline(0, color=DARK, lw=0.8)
    ax_b.set_ylabel("Mean methylation difference")
    panel(ax_b, "B", "Stage-specific effects are stronger in brain")
    clean(ax_b, "y")

    key = [
        ("E15.5", "Dnmt3a1-/-", "E15.5\n$\\mathit{Dnmt3a1}$ KO"),
        ("E15.5", "Dnmt3a2-/-", "E15.5\n$\\mathit{Dnmt3a2}$ KO"),
        ("PD21", "Dnmt3a1-/-", "P21\n$\\mathit{Dnmt3a1}$ KO"),
        ("PD21", "Dnmt3a2-/-", "P21\n$\\mathit{Dnmt3a2}$ KO"),
    ]
    x = np.arange(4)
    for offset, sx, marker, color in [(-0.15, "female", "o", PURPLE),
                                      (0.15, "male", "s", TEAL)]:
        vals = [float(next(r["mean_KO_minus_WT"] for r in sex
                           if r["stage"] == stage and r["tissue"] == "brain"
                           and r["knockout"] == ko and r["sex"] == sx))
                for stage, ko, _ in key]
        ax_c.scatter(x + offset, vals, s=70, marker=marker, color=color,
                     label=sx.capitalize(), edgecolor="white", linewidth=0.7, zorder=3)
    ax_c.axhline(0, color=DARK, lw=0.8)
    ax_c.set_xticks(x, [k[2] for k in key])
    ax_c.set_ylabel("Mean methylation difference")
    ax_c.legend(frameon=False)
    panel(ax_c, "C", "Sex-stratified estimates show the within-stage ordering")
    clean(ax_c, "y")

    metrics = {r["metric"]: float(r["value"]) for r in cpg}
    groups = ["WT", r"$\mathit{Dnmt3a1}$ KO", r"$\mathit{Dnmt3a2}$ KO"]
    vals = [metrics["WT_mean_methylation"], metrics["Dnmt3a1_KO_mean_methylation"],
            metrics["Dnmt3a2_KO_mean_methylation"]]
    ax_d.bar(np.arange(3), vals, color=[GRAY, BLUE, ORANGE], width=0.65)
    ax_d.set_xticks(np.arange(3), groups)
    ax_d.set_ylim(0.55, 0.70)
    ax_d.set_ylabel("Mean CpG methylation")
    for i, value in enumerate(vals):
        ax_d.text(i, value + 0.003, f"{value:.3f}", ha="center", fontsize=8)
    panel(ax_d, "D", "Independent P21 neuronal CpG comparison")
    clean(ax_d, "y")

    fig.suptitle("Figure 3 | Developmental allocation of DNMT3A isoform function",
                 x=0.01, ha="left", fontsize=14, fontweight="bold")
    save(fig, "Figure_3_developmental_isoform_function")


def figure4() -> None:
    reps = read_tsv("results/mch_global_summary/replicate_qc_and_methylation.tsv")
    gene = read_tsv("results/mch_gene_body_comparison/gene_body_mCA_summary.tsv")
    mecp2 = read_tsv("results/mch_gene_body_comparison/mecp2_gene_set_enrichment.tsv")
    expr = read_tsv("results/mch_gene_body_comparison/mCA_expression_effect_associations.tsv")
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.6), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.flat

    group_order = ["WT", "Dnmt3a1_KO", "Dnmt3a2_KO", "Dnmt3a1_delta_N"]
    labels = ["WT", r"$\mathit{Dnmt3a1}$ KO", r"$\mathit{Dnmt3a2}$ KO", "DNMT3A1 ΔN"]
    colors = [GRAY, BLUE, ORANGE, TEAL]
    group_means = {}
    for i, group in enumerate(group_order):
        vals = [float(r["autosomal_corrected_mCA"]) for r in reps if r["genotype"] == group]
        group_means[group] = float(np.mean(vals))
        offsets = np.linspace(-0.08, 0.08, len(vals))
        ax_a.scatter(np.full(len(vals), i) + offsets, vals, s=58, color=colors[i],
                     edgecolor="white", linewidth=0.7, zorder=3)
        ax_a.plot([i - 0.19, i + 0.19], [np.mean(vals)] * 2, color=DARK, lw=2)
    wt_mean = group_means["WT"]
    for i, group in enumerate(group_order):
        ax_a.text(i, group_means[group] + 0.00065,
                  f"{group_means[group] / wt_mean:.3f}× WT",
                  ha="center", va="bottom", fontsize=7.5)
    ax_a.set_xticks(range(4), labels)
    ax_a.set_ylabel("Lambda-corrected autosomal mCA")
    panel(ax_a, "A", "DNMT3A1 dominates bulk P21 neuronal mCA")
    clean(ax_a, "y")

    lookup = {r["contrast"]: r for r in gene}
    contrasts = ["Dnmt3a1_KO-WT", "Dnmt3a2_KO-WT", "Dnmt3a1_delta_N-WT"]
    vals = [float(lookup[c]["median_relative_to_WT"]) for c in contrasts]
    ax_b.bar(np.arange(3), vals, color=[BLUE, ORANGE, TEAL], width=0.65)
    ax_b.axhline(1, color=DARK, ls="--", lw=0.9)
    ax_b.set_xticks(np.arange(3), [r"$\mathit{Dnmt3a1}$ KO", r"$\mathit{Dnmt3a2}$ KO", "DNMT3A1 ΔN"])
    ax_b.set_ylabel("Median gene-body mCA / WT")
    ax_b.set_ylim(0, 1.25)
    for i, value in enumerate(vals):
        ax_b.text(i, value + 0.035, f"{value:.3f}", ha="center", fontsize=8)
    panel(ax_b, "B", "Gene bodies reproduce the global ordering")
    clean(ax_b, "y")

    x = np.arange(3)
    long_vals = [float(lookup[c]["median_absolute_difference_genes_ge_100kb"])
                 for c in contrasts]
    short_vals = [float(lookup[c]["median_absolute_difference_genes_lt_100kb"])
                  for c in contrasts]
    ax_c.bar(x - 0.18, long_vals, 0.36, color=PURPLE, label="≥100 kb")
    ax_c.bar(x + 0.18, short_vals, 0.36, color="#A9A9A9", label="<100 kb")
    ax_c.axhline(0, color=DARK, lw=0.8)
    ax_c.set_xticks(x, [r"$\mathit{Dnmt3a1}$ KO", r"$\mathit{Dnmt3a2}$ KO", "DNMT3A1 ΔN"])
    ax_c.set_ylabel("Median mCA difference vs WT")
    ax_c.legend(frameon=False)
    panel(ax_c, "C", r"Absolute loss is greatest in long genes after $\mathit{Dnmt3a1}$ loss")
    clean(ax_c, "y")

    selected_sets = ["core_MR", "recurrent_subclass_MR"]
    selected_contrasts = ["Dnmt3a1_KO-WT", "Dnmt3a1_delta_N-WT"]
    x = np.arange(2)
    for offset, contrast, color, label in [(-0.18, selected_contrasts[0], BLUE, r"$\mathit{Dnmt3a1}$ KO"),
                                           (0.18, selected_contrasts[1], TEAL, "DNMT3A1 ΔN")]:
        vals = [float(next(r["most_negative_decile_odds_ratio"] for r in mecp2
                           if r["contrast"] == contrast and r["gene_set"] == gs))
                for gs in selected_sets]
        ax_d.bar(x + offset, vals, 0.36, color=color, label=label)
        for i, value in enumerate(vals):
            ax_d.text(i + offset, value + 0.25, f"{value:.2f}", ha="center", fontsize=8)
    ax_d.axhline(1, color=DARK, ls="--", lw=0.9)
    ax_d.set_xticks(x, ["Core MeCP2-\nrepressed", "Recurrent MeCP2-\nrepressed"])
    ax_d.set_ylabel("Odds ratio for strongest mCA-loss decile")
    ax_d.legend(frameon=False)
    panel(ax_d, "D", "MeCP2-sensitive genes carry large absolute mCA losses")
    clean(ax_d, "y")

    fig.suptitle("Figure 4 | DNMT3A1-dependent neuronal mCA and gene context",
                 x=0.01, ha="left", fontsize=14, fontweight="bold")
    save(fig, "Figure_4_neuronal_mCA")


def figure_s2() -> None:
    rows = read_tsv(
        "structural_modeling/8QZM/md/matched_box/analysis/"
        "cross_preparation_comparison.tsv"
    )
    lookup = {r["metric"]: r for r in rows}
    previous = "previous_variant_box_500ps_delta_L188H_minus_WT_A"
    common = "common_solvent_box_500ps_delta_L188H_minus_WT_A"
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), constrained_layout=True)
    ax_a, ax_b = axes
    distance_metrics = ["R188_H2A_min_distance", "R188_H3_min_distance",
                        "R188_ubiquitin_min_distance"]
    distance_labels = ["H2A", "H3", "Ubiquitin"]
    x = np.arange(3)
    ax_a.bar(x - 0.18, [float(lookup[m][previous]) for m in distance_metrics],
             0.36, color=PURPLE, label="Variant-specific boxes")
    ax_a.bar(x + 0.18, [float(lookup[m][common]) for m in distance_metrics],
             0.36, color=TEAL, label="Common-solvent boxes")
    ax_a.axhline(0, color=DARK, lw=0.8)
    ax_a.axhspan(-0.1, 0.1, color=LIGHT, alpha=0.65, zorder=0)
    ax_a.set_xticks(x, distance_labels)
    ax_a.set_ylabel("L188H − WT minimum distance (Å)")
    ax_a.legend(frameon=False, loc="lower left")
    panel(ax_a, "A", "Interface-distance changes depend on preparation")
    clean(ax_a, "y")

    mobility_metrics = ["UDR_backbone_RMSD", "R188_CA_RMSF"]
    mobility_labels = ["UDR backbone RMSD", "Residue-188 RMSF"]
    x = np.arange(2)
    ax_b.bar(x - 0.18, [float(lookup[m][previous]) for m in mobility_metrics],
             0.36, color=PURPLE, label="Variant-specific boxes")
    ax_b.bar(x + 0.18, [float(lookup[m][common]) for m in mobility_metrics],
             0.36, color=TEAL, label="Common-solvent boxes")
    ax_b.axhline(0, color=DARK, lw=0.8)
    ax_b.set_xticks(x, mobility_labels)
    ax_b.set_ylabel("L188H − WT (Å)")
    panel(ax_b, "B", "Mobility directions reverse between preparations")
    clean(ax_b, "y")
    fig.suptitle("Figure S2 | Preparation sensitivity of the L188H simulations",
                 x=0.01, ha="left", fontsize=13, fontweight="bold")
    save(fig, "Figure_S2_L188H_preparation_sensitivity")


def main() -> None:
    setup()
    figure1()
    figure2()
    figure3()
    figure4()
    figure_s2()
    print(OUT)


if __name__ == "__main__":
    main()
