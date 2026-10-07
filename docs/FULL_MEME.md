# Alignment-wide DNMT3A MEME scan

## Authoritative corrected mammal-scope analysis

The final scan recalculates occupancy within the 250-mammal analysis scope,
retaining 909 codons. It uses the rooted mammal tree, has no site filter, and
tests all 909 codons in one correction family. The run completed successfully
on 2026-07-24.

Of 909 sites, 31 have nominal p-values at or below 0.05. Three pass
Benjamini-Hochberg correction:

| Human site | Domain | MEME p | BH q | Bonferroni p |
|---|---|---:|---:|---:|
| P5 | N-terminal regulatory | 4.36e-7 | 0.000396 | 0.000396 |
| E18 | N-terminal regulatory | 1.67e-5 | 0.00758 | 0.0152 |
| S6 | N-terminal regulatory | 1.06e-4 | 0.0320 | 0.0961 |

The corrected scan therefore reproduces the same three primary-alignment
discoveries as the earlier 901-codon analysis. It does not change their
reliability classification: P5, S6, and E18 still fail the conservative
alignment/annotation checks described below.

The retained candidates have these corrected full-scan values:

| Human site | MEME p | BH q | Interpretation |
|---|---:|---:|---|
| T12 | 0.000958 | 0.0871 | strong nominal signal; not 5% FDR |
| G34 | 0.000922 | 0.0871 | strong nominal signal; not 5% FDR |
| S97 | 0.02245 | 0.7206 | nominal signal; not corrected |

Machine-readable corrected results are under
`results/04_selection/mammals/mammal_scope/full_meme/`.

## Earlier 901-codon analysis

MEME was run without a site filter on the rooted production tree, using 250
unique mammalian coding sequences and all 901 retained codons. The run reports
`limit-to-sites => null`, confirming that it is distinct from the deprecated
numeric-prefix-filtered experiment. HyPhy estimated a global omega of 0.0400.

Of 901 sites, 29 had nominal p-values at or below 0.05. Three passed
Benjamini-Hochberg correction across all 901 codons:

| Human site | Domain | MEME p | BH q | Bonferroni p |
|---|---|---:|---:|---:|
| P5 | N-terminal regulatory | 5.13e-7 | 0.000462 | 0.000462 |
| E18 | N-terminal regulatory | 1.47e-5 | 0.00662 | 0.0132 |
| S6 | N-terminal regulatory | 1.04e-4 | 0.0313 | 0.0938 |

P5 and E18 also pass the 901-site Bonferroni threshold. All three discoveries
fall in the low-complexity DNMT3A1-specific N-terminal region.

The earlier targeted candidates do not pass the alignment-wide 5% FDR
threshold:

| Human site | Full-scan p | BH q | Conclusion |
|---|---:|---:|---|
| T12 | 0.000921 | 0.0829 | targeted candidate, not full-scan discovery |
| G34 | 0.000668 | 0.0669 | targeted candidate, not full-scan discovery |
| A114 | 0.1248 | 0.721 | unsupported |

The full-scan correction family is the complete set of 901 codons. The earlier
three-site MEME jobs remain follow-up tests selected using FUBAR and must not be
described as alignment-wide discoveries.

## Independent-alignment validation

Each full-scan discovery was mapped by human DNMT3A1 residue into independently
filtered MAFFT and PRANK codon alignments. MEME was run as one exact-site job
per candidate and alignment. Validation p-values were Holm-corrected across
the three discovery candidates within each alternative alignment.

| Human site | MAFFT MEME p | MAFFT Holm p | PRANK MEME p | PRANK Holm p | Classification |
|---|---:|---:|---:|---:|---|
| P5 | 3.80e-5 | 7.60e-5 | 0.667 | 0.965 | primary+MAFFT only; PRANK-sensitive |
| S6 | 7.16e-7 | 2.15e-6 | 0.115 | 0.345 | primary+MAFFT only; PRANK-sensitive |
| E18 | 0.484 | 0.484 | 0.483 | 0.965 | primary only |

No discovery site is MEME-significant in all three alignments.

P5 and E18 are especially sensitive to the placement of *Puma concolor* and
*Vombatus ursinus* residues. At P5, the primary alignment contains
*Vombatus* K5 and *Puma* R5, whereas PRANK places both residues outside the
human-mapped column. At E18, PRANK likewise excludes the primary *Vombatus*
P18 and *Puma* E18 states; MAFFT also excludes the *Vombatus* state. S6 has
99.2% primary-versus-PRANK amino-acid agreement and elevated PRANK FUBAR
support (0.831), but its PRANK MEME p-value remains 0.115.

## Consensus-masked sensitivity

To isolate the disputed placements, three site-specific copies of the primary
alignment were made. At each candidate, a codon was replaced by `NNN` only
when its exact human-mapped placement differed among primary, MAFFT, and
PRANK. Exactly two taxa were masked at every site: *Puma concolor* and
*Vombatus ursinus*. The other 248 taxa and all other codons were unchanged.

| Human site | Masked MEME LRT | Masked p | Holm p (three candidates) |
|---|---:|---:|---:|
| P5 | 0.000 | 0.667 | 0.969 |
| S6 | 2.831 | 0.117 | 0.350 |
| E18 | 0.318 | 0.485 | 0.969 |

None remains nominally significant after the six disputed codon placements
are treated as missing. The close agreement with the PRANK results shows that
the primary full-scan discoveries depend on the *Puma*/*Vombatus* placements,
not on a signal distributed across the 248-taxon codon consensus.

## Annotation-aware validation

The disputed taxa were traced to their local NCBI gene, RNA, CDS, and protein
records. Both are single predicted protein-coding models on unplaced scaffolds.
Their predicted RNAs map perfectly back to their source loci with coherent
multi-exon structures, confirming that the model sequences are genome-encoded.
This is not independent transcript evidence and does not establish alignment
column homology.

Mapping the human-numbered columns through the production amino-acid alignment
shows that Puma contributes R2, P3, and E15 at human P5, S6, and E18. All are
in its unique predicted first exon. Puma residue 26 then begins a 144-aa exact
match to residue 57 in four close felid representatives, whose shared upstream
prefix is absent from the Puma model. The Puma placements are therefore
excluded from homologous-site inference. Puma has a gap, not a residue, at the
human G34 column.

Wombat contributes K2, G3, P25, and G47 at P5, S6, E18, and G34. These codons
are genome-encoded and do not cross splice junctions, but the model has an
unusual N-terminal extension and no independent transcript support in the
package. Its positional homology remains unresolved and should be masked in
conservative tests. See `docs/CANDIDATE_ANNOTATION_VALIDATION.md`.

A later region-wide sensitivity generalized this mask across all 266 retained
N-terminal coordinates. T12 and G34 retain FUBAR support, but their exact-site
masked MEME p-values (0.0448 and 0.0426) are nominal and become 0.1277 after
Holm correction across three follow-ups. See
`docs/N_TERMINAL_RELIABILITY_MASK.md`.

## Interpretation

The primary alignment contains three site-wide corrected MEME signals, but
none survives the prespecified standard of significance across both
alternative alignments. They are alignment-sensitive candidates, not robust
positive-selection discoveries. The consensus-masked analysis further shows
that none persists after the two aligner-discordant taxa are removed from the
site likelihood.

G34 remains the most consistently supported hypothesis across primary,
composition-pass, MAFFT, and PRANK FUBAR analyses. However, its full-scan
BH q-value is 0.0669, its PRANK MEME result is only nominal, and targeted
aBSREL does not confirm either highlighted G34-bearing branch. The defensible
current conclusion is therefore that DNMT3A's N-terminal region contains
candidate episodic-evolution signals whose exact sites and lineages remain
alignment-sensitive. No result yet supports altered brain development.

Machine-readable outputs:

- `results/04_selection/mammals/full_meme/all_sites.tsv`
- `results/04_selection/mammals/full_meme/significant_sites.tsv`
- `results/04_selection/mammals/full_meme/summary.json`
- `results/04_selection/mammals/full_meme/discovery_alignment_sensitivity.tsv`
- `results/04_selection/mammals/full_meme/discovery_alignment_sensitivity_summary.json`
- `results/04_selection/mammals/full_meme/sensitivity/run_manifest.tsv`
- `results/04_selection/mammals/full_meme/consensus_masked/mask_manifest.tsv`
- `results/04_selection/mammals/full_meme/consensus_masked/consensus_masked_summary.tsv`
- `results/04_selection/mammals/full_meme/consensus_masked/consensus_masked_summary.json`
- `results/04_selection/mammals/candidate_annotation_validation/candidate_site_exon_audit.tsv`
- `results/04_selection/mammals/candidate_annotation_validation/n_terminal_relative_anchor_audit.tsv`
- `results/04_selection/mammals/candidate_annotation_validation/candidate_annotation_summary.json`
