# PRANK codon-alignment sensitivity

## Design and QC

PRANK v170427 was run in codon mode on the 250 unique mammalian DNMT3A CDS
records using the monotreme-rooted production topology as its fixed guide tree.
This is algorithmically independent of the MAFFT-guided, MACSE-refined primary
alignment and explicitly models insertion/deletion history.

The input contains no internal stops or frameshifts. Ten ambiguous source
codons across four taxa were represented as unknown codons during alignment and
masked before selection analysis.

PRANK produced 1,657 codon columns. Applying the same 70% occupancy threshold
and requiring a mapped human DNMT3A1 residue retains 909 codons. All 912 human
residues are recovered before filtering. The clean alignment contains no
invalid codons.

Across the 901 human positions shared with the primary alignment:

- mean per-site amino-acid agreement is 99.91%;
- 767 sites have identical states in all 250 mammals;
- 893 sites have at least 99% agreement;
- the minimum per-site agreement is 94.8%.

PRANK additionally retains human positions 22, 97, 106, 107, 116, 117, 231,
and 232, which fall below the primary alignment's occupancy threshold.

## Original candidate results

| Human site | Primary FUBAR | PRANK FUBAR | Primary MEME p | PRANK MEME p | Primary–PRANK AA agreement |
|---|---:|---:|---:|---:|---:|
| T12 | 0.9807 | 0.9531 | 0.000920 | 0.0449 | 99.2% |
| G34 | 0.9856 | 0.9648 | 0.000663 | 0.0393 | 99.6% |
| A114 | 0.8498 | 0.8606 | 0.1246 | 0.1278 | 100% |

T12 and G34 therefore retain FUBAR support above 0.95 and nominal MEME support
at `p < 0.05` under PRANK. The PRANK MEME p-values do not pass a Bonferroni
threshold of 0.0167 for the three original candidate tests, so they should be
described as nominal alignment-sensitivity support rather than a second
corrected discovery.

A114 remains unsupported.

## Alignment-sensitive lineage assignments

PRANK does not place three extreme terminal residues in the same homologous
columns as the human candidates:

- T12: *Vombatus ursinus* E and *Puma concolor* K become gaps at the
  human-mapped PRANK column.
- G34: *Carlito syrichta* M becomes a gap at the human-mapped PRANK column.

Consequently, those three specific branch substitutions must not be presented
as robust evolutionary events. The broader T12 signal—including the eutherian
S-to-T reconstruction, cetacean T-to-A change, and several rodent changes—and
the G34 lorisiform and *Rhinolophus* changes remain aligned under PRANK.

## PRANK-only candidates

PRANK FUBAR also identifies A16 and S97:

| Site | Primary | Composition pass | MAFFT | PRANK | PRANK MEME p | Status |
|---|---:|---:|---:|---:|---:|---|
| A16 | 0.744 | 0.451 | 0.799 | 0.933 | 0.279 | PRANK-only FUBAR; unsupported |
| S97 | not retained (0.912 mammal-scope) | not retained | not retained (0.913 mammal-scope) | 0.921 | 0.0242 | Four-alignment-consistent tertiary candidate |

The later mammal-scope audit resolved S97. It fell below 70% only because
MACSE and MAFFT were filtered across 475 vertebrates before mammal subsetting.
Within the 250 mammals, all four aligners give 100% occupancy and identical
states. FUBAR is 0.911–0.921 and targeted MEME p=0.021–0.024 across methods.
S97 is now a qualified tertiary hypothesis, not a corrected discovery.

## Conclusion

PRANK strengthens the site-level prioritization of G34 and provides qualified
support for T12, while weakening several previously inferred terminal-branch
changes. It does not support A114. The correct interpretation is:

1. G34 is the most alignment-robust candidate.
2. T12 remains plausible but composition- and lineage-assignment-sensitive.
3. A114 should be deprioritized.
4. A16 remains alignment-sensitive; S97 is retained as a qualified tertiary
   hypothesis after mammal-scope validation.

COBALT was not run as a production sensitivity analysis. It is protein-only
and derives much of its value from conserved-domain and local-similarity
constraints, whereas the disputed sites lie in the low-complexity N terminus.
It would not provide a codon alignment suitable for direct FUBAR or MEME
comparison without an additional back-translation step.

Machine-readable outputs:

- `results/02_alignment/sensitivity/prank/PRANK_alignment_qc_summary.json`
- `results/04_selection/mammals/prank_candidate_sensitivity.tsv`
- `results/04_selection/mammals/prank_fubar_candidates.tsv`
- `results/04_selection/mammals/prank_alignment_site_agreement.tsv`
