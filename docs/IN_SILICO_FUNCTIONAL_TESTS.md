# In-silico assessment of DNMT3A1 tail candidates

## Question

Can sequence-based computation provide independent support for functional
effects of the recurrent T12, G34, and S97 substitutions when wet-lab assays
are not available?

Two complementary screens were used. The first measures local sequence
properties of the disordered tail. The second uses masked-marginal
probabilities from two ESM-2 protein language models. Both screens compare the
prespecified candidates with 207 recurrent natural alternatives in human
DNMT3A1 residues 1–163. This natural null controls for the unusually variable
sequence context of the upstream tail.

## Disordered-tail property screen

For 21- and 41-residue windows, each substitution was evaluated for changes in
net charge, charged fraction, NCPR, charge decoration, hydropathy, sequence
entropy, serine/threonine/tyrosine fraction, glycine/proline fraction, and
simple motif proxies. Empirical comparisons were made against all recurrent
natural alternatives and alternatives matched by the human reference residue.
False-discovery correction was applied across the 22 correlated metrics within
each variant.

None of the six tier-1 substitutions had a corrected property result
(`q < 0.05`). The minimum global within-variant values were:

| Variant | Minimum global BH q |
|---|---:|
| T12S | 0.2115 |
| T12A | 0.5288 |
| T12N | 0.6346 |
| G34A | 1.0000 |
| G34S | 1.0000 |
| S97G | 1.0000 |

The substitutions change local chemistry by definition. For example, T12A,
T12N, and S97G remove a serine/threonine/tyrosine acceptor, whereas G34S adds
one. Those changes were not exceptional relative to recurrent natural
variation in the same tail.

## ESM-2 sequence-compatibility screen

Each alternative was scored in the full 912-residue human DNMT3A1 sequence as
the masked-marginal log-likelihood ratio (LLR) of the alternative relative to
the human reference residue. The analysis was repeated with
`esm2_t12_35M_UR50D` and `esm2_t30_150M_UR50D`.

The models were strongly concordant across the 207 natural alternatives
(Pearson `r = 0.9755`; Spearman `rho = 0.9729`; preference direction agreed for
196 of 207 alternatives). Candidate-panel scores were also concordant
(Spearman `rho = 0.9231`).

| Variant | 35M LLR | 150M LLR | Matched lower-tail p, 35M | Matched lower-tail p, 150M |
|---|---:|---:|---:|---:|
| T12N | -1.1356 | -1.2003 | 0.350 | 0.350 |
| G34S | -0.6876 | -0.7222 | 0.294 | 0.353 |
| G34A | -0.1703 | 0.0823 | 0.471 | 0.588 |
| S97G | 0.6052 | 0.4379 | 0.968 | 0.871 |
| T12A | 0.8143 | 1.2321 | 0.750 | 0.850 |
| T12S | 0.9827 | 1.3480 | 0.900 | 0.900 |

Neither model classified a tier-1 candidate as unusually incompatible relative
to reference-residue-matched natural alternatives. Positive LLRs for T12A,
T12S, and S97G mean only that the model assigns the alternative higher
sequence compatibility than the human residue in the human context. They do
not imply increased activity or adaptive benefit.

## Interpretation

These analyses do not provide independent evidence that T12, G34, or S97
substitutions alter DNMT3A1 function. They also do not rule out effects on
phosphorylation, turnover, chromatin association, partner binding, or
methylation. Protein language models are limited for lineage-specific
regulation in disordered regions, and the property screen cannot represent
cellular context.

The defensible conclusion is negative but informative: the candidate
substitutions are recurrent evolutionary states, not computationally
exceptional perturbations. They should remain trajectory examples and
hypotheses. They should not organize the manuscript in the absence of an
independent molecular or comparative phenotype.

## Reproducible outputs

- `scripts/idr_variant_effects.py`
- `scripts/esm_variant_scores.py`
- `scripts/summarize_esm_sensitivity.py`
- `results/07_functional_prioritization/in_silico_idr/`
- `results/07_functional_prioritization/in_silico_esm/`

