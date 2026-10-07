# DNMT3A candidate-lineage audit

## Rooting and reconstruction

The IQ-TREE production topology is unrooted. For ancestral-state
interpretation, the 250-tip HyPhy tree was rerooted on the branch separating
the sampled monotremes (*Ornithorhynchus anatinus* and *Tachyglossus
aculeatus*) from therians. Rerooting preserves all 247 informative unrooted
splits; the maximum change in any pairwise tip distance is
`1.01e-10` (rounding at Newick serialization). Validation details are in
`results/03_tree/mammals/rooting_audit.json`.

MEME was rerun on this rooted representation. The site-level conclusions are
unchanged:

| Human site | MEME p-value | Site-level conclusion |
|---|---:|---|
| T12 | 0.000920 | Targeted episodic signal |
| G34 | 0.000663 | Targeted episodic signal |
| A114 | 0.1246 | Not significant |

Sparse HyPhy ancestral codon reconstructions reproduce all 250 observed
terminal codons at each candidate site with zero mismatches.

## T12

The rooted reconstruction places an S-to-T change on the sampled eutherian
stem. Threonine is subsequently retained across most sampled eutherians, with
recurrent changes on several branches. Branches with positive-class posterior
probability at least 0.90 include:

- eutherian stem: S to T;
- cetacean stem: T to A;
- *Acomys* stem and *Mus caroli*: T to S;
- *Fukomys damarensis* and *Nannospalax galili*: independent T to N;
- *Callithrix jacchus*: T to I;
- *Dasypus novemcinctus*: T to A;
- *Puma concolor*: T to K;
- *Vombatus ursinus*: S to E;
- sampled didelphine stem: S to N.

This is a phylogenetically broad pattern rather than a change confined to one
ecological or neurodevelopmental category. The FUBAR result is also weakened
in the composition-pass sensitivity analysis, so T12 should remain secondary
to G34 until composition-aware and N-terminal alignment analyses are expanded.

## G34

Glycine is the reconstructed ancestral and predominant mammalian state. Three
branches have positive-class posterior probability at least 0.90:

- sampled lorisiform stem: G to S;
- *Carlito syrichta*: G to M;
- sampled *Rhinolophus* stem: G to A.

None of the descendant taxa for these three branches failed the IQ-TREE
nucleotide-composition screen. G34 is also the only site supported by the
primary, composition-pass, and MAFFT-projection FUBAR analyses. It is therefore
the strongest current coding-evolution candidate at the site level.

### Targeted G34 aBSREL follow-up

After excluding the alignment-sensitive *Carlito syrichta* terminal change,
the two remaining internal G34-bearing branches were audited as post hoc
aBSREL foregrounds. Both have positive branch lengths, 100/100 SH-aLRT/UFBoot
support, no composition-failing descendants, no exact sequence duplicates in
the HyPhy input, and concordant descendant G34 states in the PRANK alignment.
Neither branch is an order-constraint edge.

| Foreground | Reconstructed G34 change | aBSREL LRT | Uncorrected p | Holm p (two foregrounds) |
|---|---|---:|---:|---:|
| sampled *Rhinolophus* stem | G to A | 0.000 | 1.0000 | 1.0000 |
| sampled lorisiform stem | G to S | 3.400 | 0.06786 | 0.13572 |

Thus, aBSREL does not support branch-wide episodic diversifying selection on
either highlighted stem. The lorisiform fit assigns approximately 0.22% of
sites to a high-omega class, but its branch LRT is not significant even before
the explicit two-foreground correction. This negative result narrows the
claim: G34 retains cross-alignment site-level support, but the specific
lineages inferred from MEME branch posteriors are descriptive candidates, not
formally confirmed selected branches.

The foregrounds were chosen after observing the G34 MEME result, so this
aBSREL analysis is a post hoc concordance test rather than independent
replication. aBSREL tests selection across the whole branch and cannot
attribute a signal specifically to G34.

## A114

A114 has several reconstructed changes and some high branch posterior
probabilities, but the site-level MEME test is not significant and primary
FUBAR support is below 0.90. Branch posteriors at this site must not be treated
as evidence for episodic selection. A114 is retained as a sensitivity finding
only.

## Statistical boundary

The branch posterior threshold of 0.90 is a descriptive prioritization rule,
not a branch-wise significance test. MEME establishes a site-level episodic
signal; it does not independently prove that each highlighted branch underwent
positive selection. Formal lineage hypotheses require prespecified foreground
sets, appropriate branch or branch-site tests, and correction across tested
hypotheses.

No highlighted pattern is presently evidence of altered brain development.
The next biological test should use independently curated developmental
DNMT3A1/DNMT3A2 expression or neuronal mCH phenotypes rather than assigning
brain relevance from lineage identity alone.

## PRANK qualification

The independent PRANK codon alignment does not place the *Vombatus ursinus*
E12, *Puma concolor* K12, or *Carlito syrichta* M34 residues in the respective
human-mapped candidate columns. These three terminal substitutions are
alignment-sensitive and are excluded from robust lineage claims. The broader
T12 and G34 site signals survive PRANK; details are in
`docs/PRANK_SENSITIVITY.md`.

Machine-readable outputs:

- `results/04_selection/mammals/candidate_branch_evidence.tsv`
- `results/04_selection/mammals/candidate_substitution_events.tsv`
- `results/04_selection/mammals/candidate_high_posterior_branches.tsv`
- `results/04_selection/mammals/candidate_lineage_summary.json`
- `results/04_selection/mammals/g34_absrel/foreground_manifest.tsv`
- `results/04_selection/mammals/g34_absrel/branch_summary.tsv`
- `results/04_selection/mammals/g34_absrel/summary.json`

## Reliability-masked lineage update

Ancestral states and branch posteriors were re-extracted from the taxon-wide
N-terminal reliability-mask fits for T12 and G34. Informative terminal states
agree with the HyPhy reconstructions exactly; six T12 and four G34 terminal
states are missing after masking and are excluded from the consistency check.

For T12, the Puma and wombat terminal hypotheses disappear while nine other
change-coupled branches remain. The raw positive-class posterior exceeds 0.90
on 475 branches after masking, so posterior threshold alone is no longer
lineage-discriminating. The updated table additionally requires a reconstructed
nonsynonymous change and an informative terminal state or an internal branch.

For G34, the tarsier terminal hypothesis disappears. The *Rhinolophus* and
lorisiform stems remain, and four branches cross the descriptive threshold
only after refitting the masked model: Pecora, hippopotamus,
*R. ferrumequinum*, and Tasmanian devil. These new branches are post hoc
descriptive leads. The prior negative aBSREL results for the two retained
internal branches remain the only formal branch-wide follow-up.

See `docs/FINAL_CANDIDATE_EVIDENCE.md` for the consolidated decision table and
`results/04_selection/mammals/reliability_mask/lineages/` for machine-readable
outputs.
