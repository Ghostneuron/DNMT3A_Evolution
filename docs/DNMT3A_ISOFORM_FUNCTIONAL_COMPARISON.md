# In-silico comparison of DNMT3A1 and DNMT3A2 in brain development

## Question

The analysis tests whether DNMT3A1 and DNMT3A2 make distinguishable
contributions to the developing mouse brain methylome. Because the isoforms
share the PWWP, ADD, and catalytic domains, sequence or structure prediction
alone cannot resolve most functional differences. The informative evidence is
isoform-specific perturbation in matched tissues and developmental stages.

Three public functional-genomics studies were audited:

- GSE164265 contains matched P21 cortex and sorted cortical-neuron data from
  wild-type, `Dnmt3a1`-knockout, and `Dnmt3a2`-knockout mice.
- GSE295720 contains MM285 methylation-array measurements from E15.5 and P21
  brain and liver for both isoform-specific knockouts.
- GSE96529 contains controlled DNMT3A1 and DNMT3A2 rescue and localization
  experiments in mouse embryonic stem cells.

The audit includes 231 GEO samples. Machine-readable sample and design tables
are in `results/07_isoform_function`.

## P21 cortical-neuron CpG methylation

The first analysis reprocessed the six matched EM-seq BEDGRAPH files in
GSE164265: two wild-type, two `Dnmt3a1`-knockout, and two
`Dnmt3a2`-knockout P21 cortical-neuron-nuclei samples. The comparison was
restricted to 13,159,458 CpGs represented in all six samples.

Mean CpG methylation was:

| Genotype | Mean methylation |
| --- | ---: |
| Wild type | 0.67558 |
| `Dnmt3a1` knockout | 0.59722 |
| `Dnmt3a2` knockout | 0.67913 |

Relative to wild type, `Dnmt3a1` loss reduced mean methylation by 0.07836.
`Dnmt3a2` loss changed the mean by +0.00355. At a descriptive 0.10
methylation-loss threshold, 26.77% of common CpGs showed a DNMT3A1-only loss,
8.16% showed a DNMT3A2-only loss, and 13.55% showed a loss in both knockouts.
A stricter criterion, requiring both knockout replicates to fall at least 0.10
below both wild-type replicates, classified 14.94% as DNMT3A1-only losses and
3.01% as DNMT3A2-only losses.

The strongest 100-kb DNMT3A1-dependent windows overlap genes including
`Prdm16`, `Msi2`, `Auts2`, `Nfia`, `Gli3`, `Zfp423`, `Notch2`, and `Meis2`.
This list is exploratory. Gene length, CpG density, and window selection have
not yet been included in an enrichment model.

## Developmental comparison

The processed GSE295720 MM285 matrix contains 296,070 probes from 138 animals.
Contrasts were estimated within stage and tissue.

### Brain

| Stage | Knockout | Mean KO minus WT | Probes with loss of at least 0.10 |
| --- | --- | ---: | ---: |
| E15.5 | `Dnmt3a1` | +0.00060 | 0.48% |
| E15.5 | `Dnmt3a2` | -0.02265 | 8.80% |
| P21 | `Dnmt3a1` | -0.07031 | 25.30% |
| P21 | `Dnmt3a2` | -0.00256 | 1.20% |

The direction reverses across development. DNMT3A2 has the larger methylation
effect at E15.5, whereas DNMT3A1 has the larger effect at P21. The change from
E15.5 to P21 is -0.07091 for the DNMT3A1 knockout effect and +0.02009 for the
DNMT3A2 knockout effect.

The embryonic DNMT3A2 effect is stronger in brain than liver. At E15.5, mean
KO-minus-WT methylation is -0.02265 in brain and -0.00821 in liver. The
postnatal DNMT3A1 effect is also stronger in brain: -0.07031 in brain and
-0.03356 in liver at P21.

Sex stratification preserves the stage-specific pattern. The E15.5 brain
DNMT3A2-knockout effect is -0.02621 in females and -0.01875 in males. The P21
brain DNMT3A1-knockout effect is -0.07209 in females and -0.06828 in males.
The result is therefore not explained by the unequal female-to-male ratio in
the pooled embryonic wild-type group.

## Evidence-supported model

The two independent comparisons support a developmental division of labor:

1. DNMT3A2 contributes disproportionately to embryonic methylation
   establishment, including regulatory regions identified by the source
   study.
2. Much of the DNMT3A2-dependent embryonic deficit is repaired by P21.
3. DNMT3A1 becomes the dominant isoform for the postnatal cortical-neuron CpG
   methylome and, in the all-context raw-read reanalysis, for bulk neuronal
   mCA.
4. The DNMT3A1-specific N terminus provides a plausible targeting mechanism
   through H2AK119ub and Polycomb-associated neurodevelopmental regions.

The data are consistent with a temporal handoff from early DNMT3A2-dependent
methylation establishment to DNMT3A1-dependent postnatal targeting or
maintenance. They do not show that one isoform is universally more active.
Isoform function depends on developmental stage, cell type, and genomic
context.

This model provides a functional target for the evolutionary analysis. The
comparative question is whether mammalian lineages altered the timing or
magnitude of DNMT3A2 promoter deployment relative to the conserved DNMT3A1
program.

## Boundaries

- The deposited GSE164265 BEDGRAPH files report CpG methylation only. The
  separate raw-read Bismark reanalysis tests neuronal mCH across all eight
  libraries and is documented in `docs/MCH_FINAL_RESULTS.md`.
- The GSE295720 array samples a selected subset of the methylome and is not
  equivalent to whole-genome bisulfite sequencing.
- Probe and CpG counts describe genomic effects. They are not biological
  replicate counts and are not used as independent observations for P values.
- The analyses use processed public data. Read-depth information is absent
  from the GSE164265 BEDGRAPH files.
- The results do not connect either isoform to brain size, IQ, or
  species-specific cognitive ability.
- This functional reanalysis does not itself establish evolutionary change.
  Cross-species promoter usage and matched developmental data remain required.

The all-context analysis additionally shows that `Dnmt3a1` knockout removes
approximately 99.4% of corrected autosomal mCA, `Dnmt3a2` knockout retains the
bulk signal, and DNMT3A1 N-terminal deletion retains approximately 40% of WT.
Matched cortex RNA changes have only a very weak inverse association with
gene-body mCA loss, so a broad causal transcriptional effect is not supported.

Source-study P18 cortex chromatin tracks add a locus-level refinement. Among
mCA-rich genes, WT DNMT3A-FLAG gene-body enrichment predicts preservation of
mCA after `Dnmt3a2` loss (adjusted partial R-squared 0.115). After deletion of
the DNMT3A1 N terminus, H3K27me3-rich gene bodies show modestly greater mCA
loss; this association remains in a joint model containing WT DNMT3A occupancy
(H3K27me3 standardized beta -0.065, chromosome-block interval -0.089 to
-0.041). This is consistent with a special requirement for the N-terminal
region in Polycomb-associated chromatin, but is not evidence that DNMT3A1
directly reads H3K27me3. The single P18 ChIP profiles and two P21 methylome
animals per genotype support descriptive genomic association only.

The directly matched P21 neuronal H2AK119ub track refines, rather than proves,
this model. H2AK119ub enrichment alone predicts greater mCA loss after
N-terminal deletion, but H2AK119ub adds no independent coefficient when
H3K27me3 is included. H3K27me3 remains associated with both greater mCA loss
and relative DNMT3A occupancy loss, while measured occupancy redistribution
does not predict mCA loss gene by gene. The result is therefore best described
as N-terminal chromatin-context dependence, not a demonstrated direct
H2AK119ub recruitment mechanism.

## Reproducible outputs

- `scripts/audit_isoform_functional_datasets.py`
- `scripts/compare_isoform_neuron_methylomes.py`
- `scripts/compare_isoform_developmental_methylation.py`
- `scripts/annotate_isoform_methylation_windows.py`
- `scripts/compare_mch_gene_bodies.py`
- `scripts/model_mch_gene_body_covariates.py`
- `scripts/associate_mch_expression_effects.py`
- `scripts/download_cortex_chromatin_tracks.py`
- `scripts/extract_cortex_chromatin_signals.py`
- `scripts/model_mch_chromatin_context.py`
- `scripts/extract_p21_polycomb_signals.py`
- `scripts/model_p21_polycomb_context.py`
- `results/07_isoform_function/functional_dataset_registry.tsv`
- `results/07_isoform_function/GSE164265_neuron_nuclei_methylation/summary.json`
- `results/07_isoform_function/GSE295720_developmental_methylation/summary.json`

Primary sources:

- Gu et al. 2022, *Nature Genetics*, doi:10.1038/s41588-022-01063-6.
- Liu et al. 2026, *Communications Biology*,
  doi:10.1038/s42003-025-09311-1.
- Manzo et al. 2017, *EMBO Journal*, doi:10.15252/embj.201797038.
