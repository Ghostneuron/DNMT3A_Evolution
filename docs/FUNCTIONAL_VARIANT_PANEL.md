# DNMT3A1 functional variant panel

This panel is a follow-up module, not the manuscript's primary result. The
primary coding result is the contrast between the variable upstream tail and
the conserved 164–219 nucleosome/H2AK119ub-engagement region. All variants
below lie upstream of that known functional interval.

Because wet-lab assays are not currently available, the panel was screened
with local disordered-tail properties and two ESM-2 models. None of the tier-1
substitutions was exceptional relative to 207 recurrent natural alternatives.
The panel is therefore retained for hypothesis definition and lineage
description, not as evidence of molecular effects.

The first-pass panel prioritizes substitutions reconstructed more than once
and supported on at least one internal branch. This avoids selecting variants
solely because an isolated terminal sequence is unusual.

## Tier 1

| Construct | Evolutionary support |
|---|---|
| G34A | independent sampled Rhinolophus and Pecora internal stems |
| G34S | lorisiform internal stem and Tasmanian devil |
| T12N | sampled didelphine stem plus two rodent terminals |
| T12A | cetacean stem and armadillo |
| T12S | sampled Acomys stem and *Mus caroli* |
| S97G | independent sampled Stenodermatinae and Acomys stems |

Tier 2 contains S97N, S97T, and S97A. Tier 3 contains singleton terminal
changes G34T, G34V, T12I, and S97K.

## Interpretation of the in-silico screen

T12N and G34S receive negative ESM-2 log-likelihood ratios in both model sizes,
but neither lies in an unusual lower tail of the reference-residue-matched
natural distribution. T12A, T12S, and S97G receive positive ratios in both
models. G34A is close to neutral and changes sign between model sizes. No
property perturbation survives within-variant false-discovery correction.

These results do not establish neutrality. They show that broad sequence
compatibility and simple disordered-region chemistry do not distinguish the
candidates. The lineage trajectories remain valid; the functional
interpretation remains unresolved.

## Staged experimental design if future assays become available

### Stage 0: assay qualification

First compare wild-type DNMT3A1, DNMT3A2, and empty vector. The assay should
detect a reproducible DNMT3A1-tail-dependent difference before screening
individual substitutions. A published DNMT3A1-tail deletion or
H2AK119ub-binding-defective construct can provide a mechanistic benchmark if
its exact design is reproduced. It should not be inferred from the current
sequence analysis.

### Stage 1: evolutionary substitutions

Each substitution should be introduced individually into the same human
DNMT3A1 construct. Required controls are wild-type human DNMT3A1, human
DNMT3A2 as the control lacking the DNMT3A1-specific tail, and a matched empty
vector.

The first assays should measure:

1. steady-state protein abundance and turnover;
2. nuclear and chromatin localization;
3. tail-dependent interaction or proximity profiles;
4. global and locus-resolved methylation activity.

The first screen can use one representative substitution per site, such as
G34A, T12A, and S97G. The remaining recurrent substitutions should be added
only after a reproducible site-level effect or a clear assay phenotype. This
staging separates the question of whether a site matters from the question of
whether different lineage states have different effects.

### Stage 2: neuronal context

Because all retained sites are in a disordered tail, the experiment should not
be framed as testing a predicted rigid structural contact. For S97 variants,
nearby S105 phosphorylation can be measured as a secondary mechanistic readout,
without assuming S97 itself is phosphorylated.

Neuronal differentiation, developmental timing, or mCH assays should follow
only after Stage 1 yields a reproducible molecular phenotype. A generic
molecular phenotype would establish biochemical relevance, not brain
adaptation.

The adult human atlas provides a calibration target for that later stage:
primary motor-cortex neuronal nuclei have approximately fourfold higher median
mCH than non-neuronal nuclei, and the direction is reproduced in both
available donors. A candidate-variant experiment should therefore compare
matched neuronal and non-neuronal differentiation contexts with biological
replicates; it should not treat thousands of single cells from one donor as
independent replicate organisms.

## Outputs

- `results/07_functional_prioritization/DNMT3A1_variant_panel.tsv`
- `results/07_functional_prioritization/control_panel.tsv`
- `results/07_functional_prioritization/human_motor_cortex_assay_context.tsv`
- `results/07_functional_prioritization/summary.json`
