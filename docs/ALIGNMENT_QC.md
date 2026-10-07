# Production alignment QC

## Dataset and method

- 475 full-length, exact protein/CDS-matched DNMT3A1-like representatives.
- MAFFT v7.526 protein guide followed by exact codon projection.
- MACSE v2.07 `refineAlignment` with one conservative leaf-cut iteration
  (`optim=1`, `local_realign_init=0.1`).
- 107 of 475 evaluated local refinements were accepted by MACSE.

The initial projection is retained as an alignment-sensitivity dataset. The
MACSE-refined nucleotide alignment is the primary coordinate system.

## QC outcome

- 475 unique sequences; 2,015 codon columns before occupancy filtering.
- Human DNMT3A1 maps to all 912 expected residues before filtering.
- No frameshift, partial-gap, ambiguous, or internal-stop codons required
  masking in any representative.
- A 70% occupancy threshold retains 901 codon columns (2,703 nt).
- The retained HyPhy input has equal sequence lengths, valid `ACGT-`
  characters only, and no internal stop codons.
- All 901 retained columns map to a human DNMT3A1 residue.

Two representatives have sequence coverage below 85% of the retained columns:
*Podarcis vaucheri* (748/901) and *Ficedula albicollis* (753/901). Primary tree
diagnostics should inspect these taxa, and full-panel site/branch results must
be repeated without them as a sensitivity test. *Neopelma chrysocephalum*
(775/901) and *Chaetura pelagica* (784/901) are the next-lowest representatives
but remain above the prespecified 85% threshold.

Eleven human N-terminal residues occur in columns below 70% occupancy and are
excluded from the primary selection alignment: 1, 2, 3, 22, 97, 106, 107, 116,
117, 231, and 232. Their exclusion should be revisited in N-terminal-focused
sensitivity analyses because this regulatory region is intrinsically more
variable than the structured domains.

## Retained columns by operational human region

| Region | Retained codons |
|---|---:|
| N-terminal regulatory | 266 |
| PWWP extended | 150 |
| PWWP–ADD linker | 48 |
| ADD | 139 |
| ADD–methyltransferase linker | 19 |
| DNA methyltransferase | 279 |

Per-sequence details are in `results/02_alignment/sequence_codon_qc.tsv`; the
site-to-human mapping is in
`results/02_alignment/human_coordinate_crosswalk.tsv`.
