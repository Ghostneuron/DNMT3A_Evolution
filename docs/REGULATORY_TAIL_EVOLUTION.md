# Evolutionary architecture of the DNMT3A1 regulatory tail

## Question

The DNMT3A1 N terminus is intrinsically disordered, but it should not be
treated as one evolutionary unit. Residues 164–219 contain experimentally
supported nucleosome- and H2AK119ub-engagement functions, whereas T12, G34,
and S97 lie farther upstream. This analysis asks whether those parts of the
tail have different evolutionary profiles across mammals.

The regional coordinates are intentionally broad. They capture the interval
supported by deletion, chromatin-occupancy, and nucleosome-interaction studies;
they do not imply that every residue from 164 to 219 directly contacts
H2AK119ub or the nucleosome.

## Mammalian profile

The analysis uses 250 unique mammalian DNMT3A coding sequences and the 901
human-mapped codons retained in the mammal-scope alignment.

| Region | Retained sites | Polymorphic | Recurrently variable | Mean normalized amino-acid entropy | Primary FUBAR candidates |
|---|---:|---:|---:|---:|---:|
| Upstream disordered tail, 1–163 | 154 | 132 (85.7%) | 107 (69.5%) | 0.0766 | 2 |
| Nucleosome/H2AK119ub engagement, 164–219 | 56 | 6 (10.7%) | 5 (8.9%) | 0.0050 | 0 |
| Proximal pre-PWWP segment, 220–277 | 56 | 23 (41.1%) | 20 (35.7%) | 0.0350 | 0 |
| PWWP extended, 278–427 | 150 | 42 (28.0%) | 16 (10.7%) | 0.0106 | 0 |
| ADD, 476–614 | 139 | 16 (11.5%) | 4 (2.9%) | 0.0041 | 0 |
| DNA methyltransferase, 634–912 | 279 | 45 (16.1%) | 6 (2.2%) | 0.0021 | 0 |

Relative to residues 1–163, the 164–219 engagement region is depleted for
polymorphic sites (odds ratio 0.020, one-sided hypergeometric
`p = 3.26e-24`) and recurrently variable sites (odds ratio 0.043,
`p = 6.83e-16`). Fifty-three of its 56 retained sites also have primary FUBAR
posterior probability of at least 0.90 for negative selection.

These results support a regional model. Most mammalian amino-acid variation is
concentrated in the upstream disordered tail, while the experimentally
supported engagement region is strongly conserved. The N terminus is therefore
not uniformly labile.

## Consequences for candidate interpretation

T12, G34, and S97 all fall in residues 1–163. None lies within the supported
164–219 engagement interval. Their alignment sensitivity and evolutionary
trajectories justify targeted experiments, but their positions do not connect
them to the known H2AK119ub/nucleosome interface.

The candidates should be reported as secondary functional hypotheses:

- T12 and G34 have reproducible site-level evidence but do not pass the
  mammal-scope MEME false-discovery threshold.
- S97 is alignment-robust in mammals but has `BH q = 0.7206`.
- No candidate has a demonstrated molecular, methylomic, developmental, or
  brain phenotype.

The central coding result is the contrast between a variable upstream tail and
a conserved functional engagement region. This observation does not by itself
show positive adaptation in the variable region. Relaxed constraint, changes
in regulatory interactions, and lineage-specific selection remain alternative
explanations.

## Functional interpretation

Experiments should first determine whether the substitutions have any
molecular consequence. DNMT3A1 wild type, DNMT3A2, and an empty-vector control
remain appropriate controls. Initial measurements should prioritize protein
abundance, chromatin association, localization, interaction profiles, and
methylation activity. A neuronal or developmental system is required before
making claims about brain function.

The known engagement region provides a mechanistic benchmark, not evidence that
the upstream candidates affect the same interface. DNMT3A1-tail deletion or a
published H2AK119ub-binding-defective construct can be added as a positive
control only when its exact published design is reproduced.

## Primary literature

- [Gu et al. (2022)](https://www.nature.com/articles/s41588-022-01063-6)
  showed that the disordered DNMT3A1 N terminus supports chromatin occupancy,
  cortical methylation, and postnatal development, and identified an
  H2AK119ub-recognition mechanism.
- [Wapenaar et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11624362/)
  mapped nucleosome engagement within the N-terminal region and provides the
  basis for treating residues 164–219 as a broad functional interval.
- [Molaro et al. (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7306680/)
  reported diversifying selection in part of the primate DNMT3A N terminus and
  proposed retrotransposon conflict as a hypothesis, not an established
  mechanism.

## Reproducible outputs

- `config/human_regulatory_subregions.tsv`
- `scripts/regulatory_tail_profiles.py`
- `results/04_selection/mammals/regulatory_tail/regulatory_subregion_summary.tsv`
- `results/04_selection/mammals/regulatory_tail/regulatory_subregion_comparisons.tsv`
- `results/04_selection/mammals/regulatory_tail/engagement_region_site_metrics.tsv`
- `results/04_selection/mammals/regulatory_tail/candidate_region_context.tsv`
- `results/04_selection/mammals/regulatory_tail/summary.json`
