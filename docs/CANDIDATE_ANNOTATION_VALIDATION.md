# Annotation validation of disputed N-terminal candidates

## Question

The primary, MAFFT, and PRANK alignments disagree at the alignment-wide MEME
discoveries P5, S6, and E18 specifically for *Puma concolor* and *Vombatus
ursinus*. This audit asks two separate questions:

1. are the residues encoded by the source genome under the selected NCBI model?
2. are they defensible homologues of the human-numbered alignment positions?

Evidence for the first question does not establish the second.

## Source annotations

All three investigated records are old, single-product
`PROTEIN_CODING_MODEL` annotations on unplaced scaffolds, with no independent
transcript record in the local NCBI package.

| Species | Assembly | Annotation date | Model | Exons |
|---|---|---|---|---:|
| *Carlito syrichta* | GCF_000164805.1 | 2017-06-27 | XM_008060287.1 / XP_008058478.1 | 20 |
| *Puma concolor* | GCF_003327715.1 | 2018-07-24 | XM_025927162.1 / XP_025782947.1 | 23 |
| *Vombatus ursinus* | GCF_900497805.2 | 2019-01-18 | XM_027850369.1 / XP_027706170.1 | 22 |

Exonerate maps each predicted RNA back to its own extracted gene locus at 100%
identity. Thus the models are internally compatible with their source
assemblies. Because the predicted RNAs were generated from those same
assemblies, this is a model-consistency check, not independent expression
evidence.

## Alignment-column to exon mapping

Human residue numbers were first mapped to columns of the production MACSE
amino-acid alignment. Each species' actual ungapped protein position was then
mapped to its CDS and RNA exon. This avoids incorrectly assuming that human P5
means residue 5 in every species.

| Species | Human site | Model residue | Model position | RNA codon | Exon | Decision |
|---|---|---:|---:|---:|---:|---|
| *Puma* | P5 | R | 2 | 4–6 | 1 | exclude from homologous-site inference |
| *Puma* | S6 | P | 3 | 7–9 | 1 | exclude from homologous-site inference |
| *Puma* | E18 | E | 15 | 43–45 | 1 | exclude from homologous-site inference |
| *Puma* | G34 | gap | — | — | — | no Puma residue at this alignment column |
| *Vombatus* | P5 | K | 2 | 4–6 | 1 | unresolved model homology |
| *Vombatus* | S6 | G | 3 | 7–9 | 1 | unresolved model homology |
| *Vombatus* | E18 | P | 25 | 73–75 | 2 | unresolved model homology |
| *Vombatus* | G34 | G | 47 | 139–141 | 2 | unresolved model homology |

No mapped candidate codon crosses an exon junction.

## Close-relative protein evidence

Puma's predicted protein begins with a unique 25-aa segment. Beginning at Puma
residue 26, it matches residues 57–200 of four close felid representatives
(*Acinonyx*, *Herpailurus*, *Felis*, and *Prionailurus*) exactly for 144 aa.
The four relatives share the longer acidic and Arg/Gly-rich DNMT3A1 prefix.
The parsimonious annotation interpretation is that the Puma model replaces or
omits the conserved approximately 56-aa felid prefix. Puma R2, P3, and E15
therefore should not be treated as orthologues of human P5, S6, and E18.

Wombat also has an unusual leading sequence. Its exact downstream anchors begin
13–14 residues later than the corresponding anchors in koala, tammar wallaby,
and brushtail possum. Unlike Puma, the present records do not distinguish a
real wombat extension from an annotation error. Independent RNA evidence or a
new annotation is needed; conservative selection analyses should mask these
wombat placements meanwhile.

## Tarsier outlier follow-up

The taxon-wide mask identified *Carlito syrichta* as a third N-terminal outlier
with 22 discordant positions. All form one block spanning human residues 34–59.
MACSE and MAFFT place Carlito model residues 1–22 in these columns, while PRANK
places gaps. The codons lie within the model's first coding exon and map
perfectly to the source assembly.

The model protein is only 875 aa, compared with 908–909 aa in the primate
comparators. Conserved downstream anchors show the same offset:
`DLEKRSEPQPEEG` begins at Carlito residue 55 versus residues 89–92 in five
other primates, while `RGRLRGGLGWESS` begins at residue 130 versus 163–167.
The model therefore replaces or omits approximately 34–37 residues of the
conserved primate DNMT3A1 prefix. Its unique residues 1–22 are excluded as
positional homologues of human 34–59.

Carlito was already masked at these positions in the taxon-wide sensitivity
alignment. Consequently, the masked G34 results—FUBAR posterior 0.9616 and
MEME p=0.0426—do not depend on the tarsier placement.

## Consequence for candidate priority

| Candidate | Current status | Recommended use |
|---|---|---|
| P5 | driven by disputed Puma/wombat placements; fails PRANK and consensus masking | low-priority annotation-sensitive lead |
| S6 | strongest alternative-alignment profile, but loses MEME support after disputed placements are masked | exploratory lead only |
| E18 | fails both alternative alignments and annotation-aware masking | deprioritize |
| G34 | consistent FUBAR signal across alignments, but full-scan q=0.0669 and branch tests are negative | retain as the main coding-evolution hypothesis, not a discovery |
| T12 | targeted signal but full-scan q=0.0829 and low-complexity N-terminal context | retain as a secondary exploratory hypothesis |
| A114 | no MEME support | stop follow-up absent new evidence |

The project can continue, but the next selection step should use an
annotation-aware N-terminal reliability mask across all taxa rather than claim
P5/S6/E18 as selected sites. Functional or brain-development work should be
framed around hypotheses such as G34/T12 or regulatory DNMT3A1/DNMT3A2
evolution, not around the rejected full-scan discoveries.

That taxon-wide mask has now been completed. It retains T12 and G34 above the
FUBAR threshold, although their masked MEME results are nominal and not
significant after correction across three follow-ups. See
`docs/N_TERMINAL_RELIABILITY_MASK.md`.

## Reproducible outputs

- `results/04_selection/mammals/candidate_annotation_validation/record_manifest.tsv`
- `results/04_selection/mammals/candidate_annotation_validation/annotation_model_audit.tsv`
- `results/04_selection/mammals/candidate_annotation_validation/candidate_site_exon_audit.tsv`
- `results/04_selection/mammals/candidate_annotation_validation/n_terminal_relative_anchor_audit.tsv`
- `results/04_selection/mammals/candidate_annotation_validation/candidate_annotation_summary.json`
- `scripts/extract_candidate_annotation_records.py`
- `scripts/candidate_annotation_audit.py`
- `scripts/carlito_annotation_audit.py`
- `results/04_selection/mammals/candidate_annotation_validation/Carlito_syrichta/discordant_site_exon_audit.tsv`
- `results/04_selection/mammals/candidate_annotation_validation/Carlito_syrichta/primate_anchor_audit.tsv`
- `results/04_selection/mammals/candidate_annotation_validation/Carlito_syrichta/summary.json`
