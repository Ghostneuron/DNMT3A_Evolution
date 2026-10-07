# Final neuronal mCH reanalysis: DNMT3A1 versus DNMT3A2

## Status and scope

All eight P21 NeuN-positive cortical-neuron EM-seq libraries from GSE164265
completed the same checksum-controlled, M-bias-trimmed Bismark workflow. The
design contains two biological replicates each for WT, `Dnmt3a1_KO`,
`Dnmt3a2_KO`, and `Dnmt3a1_delta_N`. The pre-specified primary endpoint is the
lambda-corrected autosomal mCA fraction.

The complete machine-readable outputs are on the external drive under:

`results/bismark/final_summary/`

- `replicate_qc_and_methylation.tsv`
- `genotype_effect_summary.tsv`
- `final_report.json`

The underlying sample-context and WT-contrast tables are under:

`results/bismark/genotype_comparison/`

## Primary result

| Genotype | n | Corrected mCA, mean (SD) | Relative to WT | Interpretation |
|---|---:|---:|---:|---|
| WT | 2 | 0.01483 (0.00086) | 1.000 | Reference |
| `Dnmt3a1_KO` | 2 | 0.000085 (0.000120) | 0.0057 | Approximately 99.4% below WT |
| `Dnmt3a2_KO` | 2 | 0.01657 (0.00178) | 1.117 | Bulk mCA retained |
| `Dnmt3a1_delta_N` | 2 | 0.00593 (0.00096) | 0.400 | Approximately 60.0% below WT |

Both `Dnmt3a1_KO` replicates had no, or almost no, autosomal mCA remaining
above their matched lambda background. In contrast, both `Dnmt3a2_KO`
replicates retained the bulk neuronal mCA signal. Deleting the DNMT3A1
N-terminal region produced a reproducible intermediate phenotype, retaining
about 40% of WT corrected mCA.

Corrected aggregate mCH gave the same ordering: `Dnmt3a1_KO` had no signal
above matched lambda background after correction, `Dnmt3a1_delta_N` retained
about 41% of WT, and `Dnmt3a2_KO` retained the signal. Values reported as zero
are background-corrected estimates clamped at zero; they are not proof that
every biological mCH molecule is absent.

## Quality audit

Mapping efficiency ranged from 55.4% to 61.0%, and duplicate fractions ranged
from 19.9% to 22.5%. There is no obvious genotype-specific alignment or
duplication failure that explains the biological ordering. The two replicates
within each genotype agree in direction. CpG methylation also follows a graded
pattern: WT 66.0--66.1%, `Dnmt3a1_KO` 58.6--58.9%, `Dnmt3a2_KO` 66.2--66.5%,
and `Dnmt3a1_delta_N` 63.8--64.0%.

The apparent elevation of corrected mCA/mCH in `Dnmt3a2_KO` should not be
called a gain of function. Its matched lambda background is lower than in WT,
so the size of the positive difference depends on background correction. The
robust conclusion is retention, not enhancement, of bulk neuronal mCH.

## Biological conclusion

In this P21 cortical-neuron experiment, DNMT3A1 is the dominant isoform for
bulk postnatal mCA deposition, DNMT3A2 is dispensable for the bulk endpoint,
and the DNMT3A1-specific N-terminal region contributes substantially to full
activity or targeting. This provides functional support for focusing the
evolutionary study on the origin and regulation of the long DNMT3A1 form and
its N-terminal tail, while treating DNMT3A2 as a developmentally regulated
isoform whose importance may be locus-, stage-, activity-, or cell-state-
specific rather than required for the global P21 mCH level.

The experiment does not establish that DNMT3A isoform evolution changed brain
size, intelligence, or species differences in mCH. It supplies a mechanistic
bridge: evolutionary changes in DNMT3A1 expression, promoter use, or the
N-terminal tail are plausible candidates for changing neuronal methylation,
but cross-species regulatory and methylome evidence is still needed.

## Statistical limits and next analysis

With only two biological replicates per genotype, these are large, replicated
descriptive effect sizes rather than definitive population-level estimates.
Genomic bins or cytosines must not be treated as independent biological
replicates. The next defensible analyses are gene-body and chromosome-block
effects, emphasizing effect-size consistency and block resampling. Priority
comparisons are long neuronal genes, MeCP2-sensitive genes, and chromatin
contexts proposed to recruit DNMT3A1. These analyses can localize the retained
signal but cannot replace additional animals for genotype-level inference.

## Gene-body localization (secondary)

A primary-autosome complete-case analysis retained 19,702 genes having at
least 100 covered CA sites in every sample; 2,230 were at least 100 kb long.
Lambda correction was performed separately for every library before the two
replicates were averaged. Median gene-body mCA relative to WT was 0.021 for `Dnmt3a1_KO`,
1.109 for `Dnmt3a2_KO`, and 0.409 for `Dnmt3a1_delta_N`, independently
reproducing the global ordering.

The absolute DNMT3A1-dependent loss was larger in long genes. The median
`Dnmt3a1_KO - WT` difference was -0.01599 among genes at least 100 kb, versus
-0.00626 among shorter genes. For `Dnmt3a1_delta_N - WT`, the corresponding
medians were -0.00889 and -0.00372. DNMT3A2 loss retained a small positive
median difference in both length classes, again interpreted as retention
rather than enhancement because it is background-correction-sensitive.

This long-gene pattern is consistent with the concentration of neuronal mCA
across long transcription units, but it is not yet evidence that a particular
gene pathway is selectively targeted. Gene length, WT mCA abundance, coverage,
and neuronal expression are correlated. A pathway claim requires covariate-
aware analysis and an independently defined gene set.

The reproducible analysis is `scripts/compare_mch_gene_bodies.py`; outputs are
in `results/mch_gene_body_comparison/`. No gene-level p-values are reported,
because genes and cytosites are not biological replicates.

## Covariate-aware result

All 16 processed whole-cortex RNA count tables from the source study were
restored and checksum-recorded. Eight matched WT libraries provide
experiment-specific expression covariates. After NCBI Ensembl-to-symbol
mapping, 17,179 autosomal genes were shared by the methylation and expression
tables. Expression was normalized as log2(CPM + 1) and averaged only within
the WT group matched to each perturbation experiment. Whole cortex is not
cell-matched to NeuN-positive nuclei, so expression remains a supporting
covariate.

For each contrast, a standardized model predicted genotype-minus-WT corrected
gene-body mCA using WT mCA, log gene length, covered CA-site density, and
matched WT expression. Chromosome-block resampling used chr1--chr19 and tests
genomic robustness, not animal-level uncertainty.

- `Dnmt3a1_KO`: the unadjusted length coefficient was -0.317. It remained
  negative after adjustment (beta -0.094; chromosome-block interval -0.102 to
  -0.086; length partial R-squared 0.080). Thus longer genes show somewhat
  greater DNMT3A1-dependent loss even after accounting for their higher WT
  mCA, callability, and cortex expression.
- `Dnmt3a1_delta_N`: adjustment reversed the weak unadjusted trend. The length
  coefficient was +0.044 (block interval +0.030 to +0.058; partial R-squared
  0.007). Conditional on WT mCA and covariates, longer genes retain slightly
  more—not less—mCA after N-terminal deletion. This does not support selective
  failure of the truncated protein at long genes.
- `Dnmt3a2_KO`: length and expression coefficients are exploratory because the
  full model explains only 6.5% of gene-level variation and the positive
  global contrast is sensitive to lambda correction. The defensible result
  remains preservation of bulk mCA.

Matched WT expression contributed essentially no independent effect in the
DNMT3A1 knockout or N-terminal-deletion models. The broad long-gene result is
therefore not explained by baseline whole-cortex abundance, although a
cell-type-matched expression dataset would be a stronger test. These models
remain observational and are based on only two methylome animals per genotype.

The covariate workflow is `scripts/prepare_mch_cortex_expression.py` followed
by `scripts/model_mch_gene_body_covariates.py`. Results are in
`results/mch_gene_body_comparison/gene_body_covariate_models.tsv`.

## Does mCA loss track transcriptional change?

The matched perturbation RNA tables permit a cross-assay, gene-level test.
For 17,179 shared autosomal genes, treatment-minus-WT whole-cortex expression
was compared with treatment-minus-WT corrected gene-body mCA. Models adjusted
for WT expression, WT mCA, log gene length, and covered CA-site density.

The direction is compatible with mCA-associated repression, but the magnitude
is very small:

- `Dnmt3a1_KO`: Spearman rho -0.081; adjusted standardized mCA-effect beta
  -0.067; chromosome-block interval -0.100 to -0.024; partial R-squared
  0.00037.
- `Dnmt3a1_delta_N`: Spearman rho -0.084; adjusted beta -0.087; block interval
  -0.116 to -0.059; partial R-squared 0.00171.
- `Dnmt3a2_KO`: Spearman rho +0.014 and adjusted beta +0.019, providing no
  comparable inverse pattern.

Thus genes with greater DNMT3A1-dependent mCA loss tend, very weakly, toward
increased expression. Gene-body mCA change explains less than 0.2% of adjusted
expression-effect variation, so this is not a strong genome-wide functional
coupling and should not be presented as causal. Dilution in whole cortex,
cell-composition differences, CpG changes, and other direct or indirect
DNMT3A effects remain plausible explanations. The value of this analysis is
its directional consistency and its clear upper bound on the broad effect.

The reproducible workflow is
`scripts/summarize_mch_cortex_expression_effects.py` followed by
`scripts/associate_mch_expression_effects.py`. Results are in
`results/mch_gene_body_comparison/mCA_expression_effect_associations.tsv`.

## Independent MeCP2-target enrichment

To test a biologically focused hypothesis, neuronal MeCP2-regulated gene sets
were taken from Moore et al. (2025), DOI 10.1038/s41593-025-01947-w. These
classes were defined independently in L4, L5, PV, and SST neurons. The shared
analysis universe contains 386 core MeCP2-repressed genes and 147 genes called
MeCP2-repressed in at least two neuronal subclasses.

These targets are strongly enriched for high WT mCA and therefore for large
absolute DNMT3A1-dependent losses:

- Core MeCP2-repressed genes have median WT mCA 0.00930 above unchanged genes.
  Their `Dnmt3a1_KO - WT` absolute loss is 0.00980 larger, and they are
  3.63-fold enriched in the most-negative mCA-effect decile (BH q
  1.15e-22).
- Recurrently MeCP2-repressed genes have median WT mCA 0.01623 above unchanged
  genes. Their DNMT3A1-knockout loss is 0.01694 larger, with 10.11-fold
  enrichment in the most-negative decile (BH q 6.65e-22).
- The DNMT3A1 N-terminal deletion shows the same absolute ordering: core and
  recurrent MeCP2-repressed genes have losses 0.00559 and 0.01106 larger than
  controls, respectively, with 2.21- and 6.16-fold decile enrichment.

This does not demonstrate preferential fractional targeting. Gene-level
retention ratios are unstable when WT corrected mCA is near zero. In the
sensitivity subset requiring WT corrected mCA at least 0.005, block intervals
for the core and recurrent DNMT3A1-knockout retention differences include
zero. After N-terminal deletion, core MeCP2-repressed genes actually retain a
slightly larger fraction than controls (+0.0545; block interval +0.0301 to
+0.0806), while recurrent targets show no clear difference (+0.0092; interval
-0.0313 to +0.0457).

Thus DNMT3A1 supplies exceptionally abundant mCA at independently defined
MeCP2-sensitive neuronal genes, but current data do not show that the
DNMT3A1-specific tail preferentially determines the fractional methylation of
those targets. Nor do the P21 whole-cortex RNA data show their expected broad
de-repression: median expression differences for core MeCP2-repressed genes
are small and negative in both `Dnmt3a1_KO` and `Dnmt3a1_delta_N`. This is
consistent with cell-type dilution, developmental timing, or weak/contextual
MeCP2 action, and argues against an immediate genome-wide transcriptional
switch.

The source workbook, extraction audit, and enrichment outputs are recorded in
`data/raw/mecp2_gene_sets/`,
`results/mch_gene_body_comparison/moore_2025_mecp2_gene_classes.tsv`, and
`results/mch_gene_body_comparison/mecp2_gene_set_enrichment.tsv`.
The reproducible entry points are `scripts/extract_mecp2_gene_sets.py` and
`scripts/analyze_mecp2_gene_set_enrichment.py`.

## Cortex chromatin context

Four source-study P18 cerebral-cortex bigWigs (input, WT DNMT3A-FLAG,
H3K4me3, and H3K27me3) were checksum-validated and summarized over the same
mm9 gene bodies and TSS +/- 2 kb promoters. ChIP features were defined as
`log1p(ChIP mean) - log1p(input mean)`. Models predicted genotype-minus-WT
corrected gene-body mCA after adjustment for WT mCA, log gene length, covered
CA-site density, and matched WT expression. The primary sensitivity analysis
required WT corrected mCA at least 0.005.

The strongest result is a division between DNMT3A2 loss and deletion of the
DNMT3A1 N terminus:

- After `Dnmt3a2` loss, WT DNMT3A-FLAG gene-body enrichment strongly predicts
  mCA preservation (standardized beta +0.398; chromosome-block interval
  +0.349 to +0.447; partial R-squared 0.115). This is consistent with DNMT3A1
  maintaining mCA in DNMT3A-occupied domains when DNMT3A2 is absent.
- After DNMT3A1 N-terminal deletion, H3K27me3 gene-body enrichment predicts
  modestly greater mCA loss (beta -0.048; interval -0.069 to -0.032; partial
  R-squared 0.0061), whereas DNMT3A-FLAG enrichment predicts modest
  preservation (beta +0.024; interval +0.004 to +0.041).
- In a joint gene-body model, both signals remain: for the N-terminal deletion,
  DNMT3A-FLAG beta is +0.046 (+0.023 to +0.069) and H3K27me3 beta is -0.065
  (-0.089 to -0.041), with joint partial R-squared 0.0110. Thus the Polycomb-
  context association is not explained solely by its moderate correlation
  with WT DNMT3A occupancy.
- For complete `Dnmt3a1` knockout, the corresponding joint coefficients are
  small and their block intervals include zero after WT mCA adjustment. This
  is expected because the knockout removes nearly all measurable mCA, leaving
  little locus-selective residual variation beyond baseline WT mCA.

These results provide an in-silico mechanistic bridge between the DNMT3A1
N-terminal region and chromatin context, but not proof of direct H3K27me3
recognition. The P18 ChIP profiles are single source-study tracks, the P21
methylomes have two animals per genotype, and chromosome-block intervals test
genomic robustness rather than animal-level uncertainty. H3K27me3 may mark a
correlated Polycomb domain feature. The P18 subset does not contain a matched
H2AK119ub track; a directly age- and cell-matched P21 analysis is reported
below. Neither subset has sufficient biological replication for a causal
claim.

The reproducible entry points are
`scripts/download_cortex_chromatin_tracks.py`,
`scripts/extract_cortex_chromatin_signals.py` and
`scripts/model_mch_chromatin_context.py`. Outputs are
`results/mch_gene_body_comparison/cortex_chromatin_gene_signals.tsv.gz`,
`results/mch_gene_body_comparison/mch_chromatin_context_models.tsv`, and
`results/mch_gene_body_comparison/mch_dnmt3a_h3k27me3_joint_models.tsv`.

## Matched P21 neuronal H2AK119ub test

The source study also deposited P21 NeuN-positive cortical-nuclei input,
H2AK119ub, H3K27me3, and H3K4me3 tracks. These are precisely age- and
cell-fraction-matched to the neuronal methylomes. Separate P21 whole-cortex
anti-DNMT3A tracks compare WT with the DNMT3A1 N-terminal deletion. All six
bigWigs were size- and SHA-256-validated before gene-body and promoter
summarization.

Among genes with WT corrected mCA at least 0.005, single-feature models show
that both Polycomb marks identify greater loss after N-terminal deletion:

- H2AK119ub gene-body enrichment: standardized beta -0.0394, chromosome-block
  interval -0.0562 to -0.0241, partial R-squared 0.00445.
- H3K27me3 gene-body enrichment: beta -0.0572, interval -0.0774 to -0.0390,
  partial R-squared 0.00923.

This does not isolate H2AK119ub recognition. H2AK119ub and H3K27me3 gene-body
features are strongly correlated (Pearson r 0.736). In their joint model,
H2AK119ub contributes essentially no independent coefficient (beta -0.0020;
interval -0.0258 to +0.0195), while H3K27me3 remains negative (beta -0.0558;
interval -0.0821 to -0.0298). Adding H3K4me3 leaves independent negative
coefficients for H3K27me3 (-0.0697) and H3K4me3 (-0.0575), but not H2AK119ub.
Thus the result supports sensitivity to chromatin context, with the strongest
signal in H3K27me3-rich domains, rather than H2AK119ub-specific targeting.

The whole-cortex DNMT3A tracks show a related redistribution. Separately,
H2AK119ub and H3K27me3 predict lower delta-N-versus-WT DNMT3A occupancy
(betas -0.0959 and -0.1635). In a joint model, H3K27me3 remains strongly
negative (beta -0.2430), whereas the conditional H2AK119ub coefficient reverses
to positive (+0.1175), reflecting their overlap and distinct conditional
components. Most importantly, gene-wise DNMT3A occupancy change does not
predict the N-terminal-deletion mCA effect in the mCA-rich subset (beta
+0.0039; interval -0.0103 to +0.0194). The data therefore do not establish an
H2AK119ub-to-DNMT3A-occupancy-to-mCA mediation chain.

These are single-profile genomic associations. The neuronal histone tracks
are well matched in age and cell fraction, but the DNMT3A occupancy comparison
uses whole cortex and has one track per genotype. Block resampling measures
robustness across chromosomes, not biological replication.

Reproducible outputs are
`results/mch_gene_body_comparison/p21_polycomb_mca_single_models.tsv`,
`results/mch_gene_body_comparison/p21_polycomb_mca_joint_models.tsv`,
`results/mch_gene_body_comparison/p21_dnmt3a_occupancy_redistribution_models.tsv`,
and `results/mch_gene_body_comparison/p21_polycomb_feature_correlations.tsv`.
The entry points are `scripts/extract_p21_polycomb_signals.py` and
`scripts/model_p21_polycomb_context.py`.
The deleted raw tracks can be restored and checksum-validated with
`scripts/download_cortex_chromatin_tracks.py --manifest
config/p21_neuronal_polycomb_tracks.tsv --output-dir
data/raw/isoform_functional_genomics/GSE164265/p21_polycomb_chromatin`.
