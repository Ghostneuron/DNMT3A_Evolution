# S97 mammal-scope homology audit

## Why S97 was originally missing

S97 was absent from the production MACSE and MAFFT mammal subsets because
codon occupancy was calculated across all 475 vertebrates before the mammalian
subset was extracted. Its occupancy was 0.644 under MACSE and 0.577 under
MAFFT across that broad dataset, below the 0.70 filter.

This was a taxon-scope filtering artifact. When occupancy is recalculated for
the 250 mammals actually used in the selection analysis, S97 has 100%
occupancy under MACSE and MAFFT. PRANK and COBALT also give 100% occupancy.
All four alignments assign exactly the same amino acid to every mammal at this
coordinate.

Eight human sites are restored by correct mammal-scope filtering: D22, S97,
P106, A107, A116, E117, Q231, and G232. Only S97 reaches FUBAR posterior
≥0.90 under all four alignments. Thus the corrected filter does not introduce
a larger set of hidden positive-selection candidates.

## Selection sensitivity

| Alignment | Mammal occupancy | FUBAR posterior | MEME p |
|---|---:|---:|---:|
| mammal-scope MACSE | 1.000 | 0.9116 | 0.0225 |
| mammal-scope MAFFT | 1.000 | 0.9134 | 0.0226 |
| PRANK | 1.000 | 0.9212 | 0.0242 |
| COBALT | 1.000 | 0.9112 | 0.0212 |

The near-identical results are expected because the four methods provide the
same S97 character vector on the same tree. They demonstrate robustness to
alignment and filtering, but they are not four statistically independent
replications.

The mammalian states are S in 239 taxa, G in four, A in two, N in two, T in
two, and K in one. Several alternatives occur in related taxon pairs:
*Acomys minous*/*A. russatus* are G, *Artibeus jamaicensis*/*Sturnira
hondurensis* are G, *Mus pahari*/*Micromys minutus* are N, and both sampled
monotremes are A. This pattern is more consistent with genuine recurrent or
inherited substitutions than with isolated terminal annotation errors.

## Decision

S97 is promoted to the third-ranked, qualified coding-evolution hypothesis
behind G34 and T12. It is suitable for ancestral reconstruction and
comparative follow-up.

Ancestral reconstruction infers serine at the mammalian root, ten codon-change
events, and eight nonsynonymous events. Eight branches combine a reconstructed
amino-acid change with MEME positive-class posterior at least 0.90: the sampled
Stenodermatinae stem (S→G), hedgehog (S→K), ring-tailed lemur (S→T), sampled
*Acomys* stem (S→G), *Mus pahari* and *Micromys minutus* independently
(S→N), naked mole-rat (S→T), and the monotreme stem (S→A). These are
descriptive change-coupled leads, not branch-wide significance tests.

It is not a corrected positive-selection discovery. The completed 909-codon
mammal-scope scan gives MEME p=0.02245 but BH q=0.7206. The three additional
exact-site alignment sensitivities are also nominal and reuse the same taxa
and tree. No brain-development phenotype can yet be attributed to the site.

## Outputs

- `results/02_alignment/sensitivity/mammal_scope/`
- `results/04_selection/mammals/mammal_scope/`
- `results/04_selection/mammals/s97_homology_audit/summary.json`
- `results/04_selection/mammals/s97_homology_audit/s97_taxon_states.tsv`
- `results/04_selection/mammals/s97_homology_audit/s97_selection_sensitivity.tsv`
- `results/04_selection/mammals/s97_homology_audit/pairwise_s97_alignment_agreement.tsv`
- `results/04_selection/mammals/s97_homology_audit/lineages/summary.json`
- `results/04_selection/mammals/s97_homology_audit/lineages/prioritized_change_coupled_branches.tsv`
- `results/04_selection/mammals/mammal_scope_filter_audit/restored_human_sites.tsv`
- `results/04_selection/mammals/mammal_scope_filter_audit/summary.json`
