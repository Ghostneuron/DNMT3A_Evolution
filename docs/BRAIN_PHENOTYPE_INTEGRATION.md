# Brain-development phenotype integration

## Current status

The phenotype layer has started with one acquired quantitative dataset and one
metadata-verified developmental transcriptome resource. It is deliberately
kept downstream of the coding-evolution results.

The acquired BrainZoo table provides global brain WGBS summaries for 13 samples
from 12 distant animal lineages. Six species occur in the DNMT3A sequence panel:
human, mouse, opossum, platypus, chicken, and great tit. Four are mammals.
These values are suitable for a descriptive pilot, not a powered mammalian
phylogenetic association.

Cardoso-Moreira et al. provide developmental RNA-seq for brain/cerebrum and
cerebellum in human, rhesus macaque, mouse, rat, rabbit, opossum, and chicken.
DNMT3A CPM has now been extracted for all six mammals in the mammalian DNMT3A
panel: 555 brain/cerebellum records and 176 tissue-age summaries. The processed
tables and sample-description files are recorded under the ArrayExpress
accessions in `data/provenance/brain_phenotype_sources.tsv`.

The 2026 Human Cell Epigenome Atlas now supplies an orthogonal adult-human
cell-type reference: 86,689 nuclei from 16 tissues, including 8,225 primary
motor-cortex nuclei. The project imports only the compact processed metadata
and subtype annotations. It reports tissue-, motor-cortex class-, and
motor-cortex subtype-level mCH/mCG summaries under
`results/05_brain_integration/`. This is useful mechanistic context for the
neuronal mCH phenotype, but it is intentionally excluded from the
cross-species association table.

All 206 published subtype labels map to the companion annotation table. In
primary motor cortex, 4,987 neuronal nuclei have median mCH fraction 0.05571,
whereas the 3,238 non-neuronal nuclei have median 0.01358 (4.10-fold
separation). Excitatory and inhibitory neurons show similar class medians
(0.05506 and 0.05721), both well above astrocytes (0.01966),
oligodendrocytes (0.01340), and microglia (0.00974). These are descriptive
cell-level distributions from two adult donors. The neuronal/non-neuronal
median contrast is present separately in both donors: 3.82-fold in H1930001
and 4.61-fold in H1930002. This donor-stratified consistency supports the
cell-type interpretation, but two donors are still not independent species
and the comparison is not a test of any DNMT3A amino-acid state.

An independent calibration against the BrainZoo adult human middle-frontal
gyrus WGBS summary is directionally coherent. Its coverage-weighted bulk mCH
fraction (0.02069) lies between the atlas non-neuronal median (0.01358) and
neuronal median (0.05571). This is compatible with a mixed-cell bulk sample,
but it cannot be converted into a neuronal-fraction estimate because the
studies differ in cortical region, donors, assay, coverage weighting, and
summary statistic.

## Important measurement boundaries

- BrainZoo WGBS mCA/mCH values may be compared only with explicit adjustment
  for brain region, age, and source study.
- The 580-species RRBS atlas is primarily a CpG resource and must not be treated
  as a broad neuronal mCH dataset.
- The Cardoso-Moreira gene-level expression matrices measure total `DNMT3A`.
  They do not by themselves distinguish DNMT3A1 from DNMT3A2. Isoform-specific
  inference requires transcript- or promoter-resolved quantification.
- Adult global mCH abundance and developmental mCH accumulation rate are
  different phenotypes and must remain separate.
- The Human Cell Epigenome Atlas is adult and human-only. Its single-cell
  resolution strengthens neuronal-versus-glial interpretation but adds zero
  independent species to a residue-state association. It also contains no
  isoform-resolved DNMT3A expression measurement.
- Species-specific developmental clocks are retained. In particular, the
  opossum processed matrix expresses postnatal stages on a conception-relative
  scale; its 14-day gestation offset is reversed before SDRF age matching.
- Candidate substitutions T12, G34, and S97 are DNMT3A1-specific N-terminal
  positions. Any connection to expression or methylation remains a hypothesis,
  not a consequence of co-occurrence in a lineage.

## Reproducible commands

```bash
python scripts/import_brainzoo_mch.py
python scripts/import_human_cell_epigenome_atlas.py
python scripts/calibrate_human_mch_context.py
python scripts/build_brain_phenotypes.py
python scripts/validate_brain_phenotypes.py
python scripts/summarize_brain_expression.py
python scripts/expression_phase_contrasts.py
python scripts/expression_trajectory_tests.py
python scripts/candidate_expression_overlay.py
python scripts/candidate_brain_phenotype_coverage.py
python -m unittest discover -s tests
```

The importer preserves the public table under `data/raw/brain`, maps only taxa
present in the sequence panel, and records excluded taxa in
`results/05_brain_integration/brainzoo_import_audit.tsv`. The validator checks
schema completeness, finite values, fraction ranges, exact taxon matching, and
duplicate measurements.

## Analysis gate

The validator requires at least 10 matched mammalian species before labeling
even a pilot phylogenetic association for a specific phenotype as ready. This
is a workflow gate, not a claim that 10 species guarantees adequate power. The
current brain mCH phenotypes each have four mammals, and total-DNMT3A
developmental expression has six; all fail that gate.

The next analysis should model within-species standardized trajectories while
keeping cerebrum and cerebellum separate and representing age on comparable
developmental scales. DNMT3A1/DNMT3A2 analysis should follow as a separate
transcript-aware module; the present gene-level CPM matrices cannot answer that
isoform question.

`expression_phase_contrasts.tsv` provides only a descriptive check of broad
prenatal-to-postnatal shifts on log2(CPM + 1). It is not a phylogenetic model:
age distributions and the anatomical definition of the brain series differ
between phases and species. Ten of twelve species-by-region contrasts decrease;
the two modest increases are the mouse brain and cerebellum series. This pattern
is a trajectory-QC observation, not evidence for adaptive DNMT3A evolution.

The monotonicity analysis uses within-species ordinal stage order, stage-median
log2(CPM + 1), 20,000 deterministic label permutations per tissue trajectory,
and Benjamini-Hochberg correction across the 12 species-by-region tests. It
does not compare raw CPM magnitudes between species. Ten trajectories have
negative Spearman correlations and remain significant at BH q < 0.05. Mouse
brain and cerebellum show small, nonsignificant positive correlations. These
are total-DNMT3A developmental patterns, not isoform-specific trajectories.

Candidate overlay is possible only descriptively. G34 and S97 are invariant
among the six expression species, while the non-reference T12 state is
represented only by opossum. Mammal-scope MACSE and PRANK agree on all three
states, but there is no replicated residue-state contrast suitable for
association testing.
Opossum has decreasing trajectories in both regions, a pattern also seen in
four of the five T12-reference species, so the available data do not support a
T12-specific expression interpretation.

Across all current phenotype types, no candidate/phenotype combination has
both the ten-mammal pilot minimum and at least two species in each of two
residue-state groups. Of the 11 sampled mammals with non-serine S97 states,
only platypus has any current brain phenotype measurement. The explicit
acquisition gaps are listed in
`results/05_brain_integration/s97_nonreference_phenotype_gaps.tsv`; paired
lineages such as *Acomys*, stenodermatine bats, and monotremes would be
especially informative if comparable neuronal mCH or developmental-expression
data become available.
