# Taxon-wide DNMT3A N-terminal reliability mask

## Purpose

The P5/S6/E18 audit showed that a very small number of uncertain N-terminal
placements can dominate a site-selection likelihood. The same rule was
therefore applied across the complete operational DNMT3A1 N-terminal region
(human residues 1–277), rather than only to sites already highlighted by
selection tests.

This is an alignment-placement sensitivity analysis. It does not independently
validate transcript expression and cannot detect every annotation error.

## Mask rule

The production MACSE, pre-refinement MAFFT, and phylogeny-aware PRANK codon
alignments contain the same 250 unique mammalian sequences. For every retained
N-terminal coordinate shared among their human-coordinate crosswalks:

1. retrieve each taxon's codon mapped to that human residue in all three
   alignments;
2. retain the primary codon only when all three codons are identical;
3. replace an informative discordant primary codon with `NNN`;
4. record discordant primary gaps or ambiguous codons without altering them.

The original production alignment remains unchanged. The masked copy retains
the same 250 taxa, 901 codons, taxon order, and sequence lengths.

## Scale of disagreement

| Quantity | Result |
|---|---:|
| Shared retained human N-terminal coordinates | 266 |
| Taxa | 250 |
| Taxon-site cells audited | 66,500 |
| Discordant cells | 138 (0.208%) |
| Newly masked informative codons | 113 (0.170%) |
| Discordant cells already missing in primary | 25 |

Of the 266 sites, 258 have at least 99% three-alignment agreement, six have
95–99%, and two have less than 95%. Human residues 108 and 109 are the two
lowest-agreement sites (237/250, 94.8%), primarily because PRANK places a
marsupial/monotreme segment differently.

The lowest-agreement taxa are:

| Taxon | Equal positions | Agreement | Discordant | Newly masked |
|---|---:|---:|---:|---:|
| *Puma concolor* | 228/266 | 85.7% | 38 | 24 |
| *Carlito syrichta* | 244/266 | 91.7% | 22 | 22 |
| *Vombatus ursinus* | 244/266 | 91.7% | 22 | 22 |

This confirms that Puma and wombat are broader N-terminal reliability outliers,
not isolated anomalies at P5/S6/E18. The tarsier annotation inspection is also
complete: all 22 discrepancies form a human 34–59 block in which MACSE/MAFFT
place unique Carlito model residues 1–22 and PRANK places gaps. Conserved
primate anchors indicate that the 2017 predicted model omits or replaces
approximately 34–37 upstream residues. The complete block is therefore
excluded from positional-homology inference.

## Candidate behavior after masking

FUBAR was rerun on the complete 901-codon masked alignment using the same
250-tip tree. Only three sites have posterior probability at least 0.90:

| Human site | Primary FUBAR | Masked FUBAR | MAFFT | PRANK | Status |
|---|---:|---:|---:|---:|---|
| T12 | 0.9807 | 0.9494 | 0.9480 | 0.9531 | retained after mask |
| A16 | 0.7442 | 0.9338 | 0.7992 | 0.9332 | mask-sensitivity-emergent |
| G34 | 0.9856 | 0.9616 | 0.9854 | 0.9648 | retained after mask |

P5, S6, and E18 have masked FUBAR posteriors of 0.0068, 0.8190, and
0.2690, respectively. A114 remains below threshold at 0.8464.

Exact-site MEME follow-ups were run for the two retained candidates and the
sensitivity-emergent A16:

| Human site | MEME LRT | MEME p | Holm p (3 tests) | Interpretation |
|---|---:|---:|---:|---|
| T12 | 4.667 | 0.0448 | 0.1277 | nominal only |
| A16 | 1.264 | 0.2732 | 0.2732 | unsupported |
| G34 | 4.765 | 0.0426 | 0.1277 | nominal only |

T12 and G34 therefore remain the most defensible coding-evolution hypotheses,
but the globally masked MEME evidence is nominal and does not survive the
three-test correction. A16 should not be promoted: it crosses the FUBAR
threshold only after masking and lacks MEME support.

G34's masked result is specifically independent of the tarsier artifact:
Carlito is already `NNN` at G34 in the masked alignment. Its retained FUBAR
posterior and nominal MEME result therefore reflect the remaining taxa.

## Interpretation

High site-wide agreement does not guarantee a robust selection result. P5,
S6, E18, and T12 each have 248/250 codon agreement, yet the two discordant taxa
have substantial likelihood leverage. Reliability must therefore be evaluated
by rerunning the evolutionary model after masking, not by imposing an arbitrary
agreement percentage alone.

The project can continue with G34 and T12 as explicitly qualified candidates.
Neither currently supports a claim of adaptive brain evolution. Candidate
lineages, close-relative annotations, and functional consequences remain to be
validated independently.

## Reproducible outputs

- `scripts/build_n_terminal_reliability_mask.py`
- `scripts/summarize_n_terminal_mask_fubar.py`
- `scripts/summarize_n_terminal_mask_meme.py`
- `results/02_alignment/reliability_mask/DNMT3A_mammals_HyPhy_unique_N_terminal_consensus_masked.fasta`
- `results/02_alignment/reliability_mask/n_terminal_mask_manifest.tsv`
- `results/02_alignment/reliability_mask/site_reliability.tsv`
- `results/02_alignment/reliability_mask/taxon_reliability.tsv`
- `results/02_alignment/reliability_mask/summary.json`
- `results/02_alignment/reliability_mask/candidate_fubar_comparison.tsv`
- `results/02_alignment/reliability_mask/fubar_summary.json`
- `results/02_alignment/reliability_mask/candidate_meme_comparison.tsv`
- `results/02_alignment/reliability_mask/meme_summary.json`
- `results/04_selection/mammals/DNMT3A.n_terminal_reliability_mask.FUBAR.json`
- `results/04_selection/mammals/reliability_mask/DNMT3A.T12.MEME.json`
- `results/04_selection/mammals/reliability_mask/DNMT3A.A16.MEME.json`
- `results/04_selection/mammals/reliability_mask/DNMT3A.G34.MEME.json`
