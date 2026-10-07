# COBALT alignment sensitivity

## Scope

NCBI COBALT 3.0.0 was added as a fourth protein-alignment sensitivity method
for the same 250 nonredundant mammalian DNMT3A sequences used by HyPhy. COBALT
was run in `-norps T` mode, so the result uses COBALT's local-similarity
constraints but not a downloaded Conserved Domain Database. This is therefore
an alignment-algorithm sensitivity analysis, not a full CDD-constrained COBALT
analysis.

COBALT emitted 250 aligned proteins across 1,219 columns. Exact curated CDS
codons were threaded through that alignment, invalid codons were masked, and
human-mapped columns with occupancy at least 0.70 were retained. The resulting
HyPhy alignment contains 909 codons. All sequences retained their exact
ungapped amino-acid content, and human DNMT3A1 maps to 912 residues.

## Candidate results

| Site | Primary–COBALT AA agreement | COBALT FUBAR | COBALT MEME p | Holm p (8) | Interpretation |
|---|---:|---:|---:|---:|---|
| P5 | 0.992 | 0.0337 | 0.6667 | 0.6667 | unsupported |
| S6 | 0.992 | 0.8627 | 0.0831 | 0.3323 | unsupported by COBALT MEME |
| T12 | 0.992 | 0.9517 | 0.00402 | 0.0241 | COBALT-supported qualified candidate |
| A16 | 0.992 | 0.8121 | 0.3282 | 0.6563 | unsupported |
| E18 | 0.992 | 0.5250 | 2.87e-5 | 2.01e-4 | statistically supported but homology/annotation-sensitive |
| G34 | 0.996 | 0.9634 | 2.45e-6 | 1.96e-5 | strongest COBALT-supported candidate |
| S97 | 1.000 after mammal-scope refiltering | 0.9112 | 0.0212 | 0.1061 | qualified tertiary candidate |
| A114 | 1.000 | 0.8651 | 0.1159 | 0.3477 | unsupported |

The N-terminal primary–COBALT disagreements are concentrated in the already
audited Puma and wombat models. At G34, the sole disagreement is Puma: COBALT
places a proline where the primary alignment has a gap. These discrepancies
explain why COBALT cannot replace the taxon-wide reliability mask.

## Interpretation

COBALT strengthens the decision to carry G34 and T12 forward. Both exceed the
FUBAR 0.90 threshold and remain significant after Holm correction across the
eight exact-site COBALT follow-ups. Their authoritative status remains
“qualified hypotheses,” because their masked MEME results do not survive the
separate three-test correction and no branch or brain phenotype effect has
been established.

E18 is not promoted. Although COBALT gives a strong MEME result, its low FUBAR
posterior and prior annotation audit show that disputed first-exon homology
still affects the inference. P5, S6, A16, and A114 gain no robust support.

The subsequent S97 audit showed that its original MACSE/MAFFT omission was
caused by applying the occupancy filter across 475 vertebrates before mammal
subsetting. Mammal-scope MACSE and MAFFT both retain S97 at 100% occupancy and
give the same residue assignment for all 250 taxa as PRANK and COBALT. S97 is
therefore promoted to a qualified tertiary candidate, while its nominal
targeted MEME evidence remains explicitly uncorrected.

## Reproducible outputs

The alignment preparation and projection commands are:

```bash
python3 scripts/prepare_cobalt_sensitivity.py prepare
tools/cobalt-3.0.0/bin/cobalt \
  -i results/02_alignment/sensitivity/cobalt/DNMT3A_mammals_unique_proteins.faa \
  -norps T -outfmt mfasta \
  > results/02_alignment/sensitivity/cobalt/DNMT3A_COBALT_AA.fasta
python3 scripts/prepare_cobalt_sensitivity.py clean
python3 scripts/cobalt_sensitivity_qc.py
```

- `results/02_alignment/sensitivity/cobalt/COBALT_alignment_qc_summary.json`
- `results/02_alignment/sensitivity/cobalt/COBALT_human_coordinate_crosswalk.tsv`
- `results/02_alignment/sensitivity/cobalt/DNMT3A_COBALT_clean_mammals_HyPhy_unique.fasta`
- `results/04_selection/mammals/cobalt/DNMT3A.COBALT.FUBAR.json`
- `results/04_selection/mammals/cobalt/cobalt_alignment_site_agreement.tsv`
- `results/04_selection/mammals/cobalt/cobalt_candidate_sensitivity.tsv`
- `results/04_selection/mammals/cobalt/cobalt_sensitivity_summary.json`
