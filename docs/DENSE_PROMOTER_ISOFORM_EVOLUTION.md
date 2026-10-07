# Dense mammalian analysis of DNMT3A promoter and isoform evolution

## Question

Does expanded taxon sampling identify a human-specific DNMT3A promoter change,
a broader mammalian gain or loss, or a phylogenetically informative origin of
the DNMT3A2-like internal transcript?

Two analyses address this question:

1. an annotation-aware screen of 253 mammalian species from 23 orders;
2. a human-branch substitution screen for the DNMT3A1 and DNMT3A2 promoter
   windows using the UCSC hg38 470-way alignment.

## Phylogenetic distribution of downstream-start products

Strict DNMT3A2-like products occur in 9 of 18 sampled placental orders.
Allowing leaderless and extended-leader downstream-core products expands the
distribution to 13 placental orders.

Two marsupial species also have downstream-core products:

- tammar wallaby, order Diprotodontia;
- fat-tailed dunnart, order Dasyuromorphia.

Their evidence is stronger than protein prediction alone. Targeted transcript
analyses recovered junction-spanning RNA evidence in both species, and
developing dunnart neocortex contains replicated internal-transcript signal.
The eutherian and marsupial products enter the shared DNMT3A coding structure
in the same narrow region.

This distribution is compatible with a DNMT3A internal-transcript architecture
in the therian ancestor. It does not prove a single ancestral promoter.
Convergent recruitment of the same coding neighborhood remains possible
because cross-therian promoter-sequence orthology has not been demonstrated.

## Why annotation non-detection cannot define losses

Downstream-core detection increases from 54.9% among species with fewer than
eight annotated DNMT3A protein products to 87.3% among species with at least
eight products. Apparent absence is therefore strongly associated with
annotation depth.

Neither of the two monotreme species in the package has a downstream-core
product. The targeted platypus RefSeq locus also lacks an annotated internal
transcript. These are annotation non-detections, not evidence that monotremes
biologically lack the isoform. No lineage-specific loss is inferred from this
matrix.

## Human-branch promoter screen

The [UCSC hg38 470-way alignment](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/)
was queried over each promoter and approximately 11 kb of local flanking
sequence. A base was callable only when chimpanzee, bonobo, gorilla, and
orangutan all carried the same nucleotide. A human difference from that
four-ape consensus was counted as a conservative fixed human-branch
substitution candidate.

### DNMT3A1 full-length TSS-cluster promoter

- callable promoter bases: 3,051;
- human-branch substitutions: 10;
- promoter substitution fraction: 0.00328;
- callable local-flank bases: 10,183;
- flank substitutions: 43;
- flank substitution fraction: 0.00422;
- promoter/flank rate ratio: 0.776;
- one-sided enrichment `p = 0.811`.

### DNMT3A2 internal promoter

- callable promoter bases: 2,337;
- human-branch substitutions: 5;
- promoter substitution fraction: 0.00214;
- callable local-flank bases: 10,962;
- flank substitutions: 38;
- flank substitution fraction: 0.00347;
- promoter/flank rate ratio: 0.617;
- one-sided enrichment `p = 0.896`.

Neither promoter has an elevated human-branch substitution fraction relative
to its aligned local flanks. This screen does not support human-specific
promoter acceleration.

The test is conservative and limited to fixed substitutions. It does not test
insertions, deletions, regulatory epigenetics, individual functional variants,
or changes outside the selected windows. It is not a formal branch-specific
phyloP analysis and is not a test of positive selection.

## Evolutionary interpretation

The expanded analysis supports mammalian regulatory evolution, but not a
human-specific DNMT3A innovation.

The strongest lineage-level hypothesis is that an internal DNMT3A transcript
architecture was present in the therian ancestor and was subsequently
remodeled. Exact promoter sequences and individual TSS positions have turned
over, whereas the coding-entry neighborhood and the separation between
full-length and internal transcript architectures have been retained.

The current data do not support:

- accelerated human evolution of either promoter;
- human-specific acquisition of DNMT3A1 or DNMT3A2 promoter activity;
- confident promoter gains or losses among mammalian orders;
- an association between promoter evolution and brain mCH.

## Manuscript decision

This result can support a comparative mammalian paper about conservation and
regulatory remodeling. It does not support a paper centered on human-specific
brain evolution.

The promoter analysis should be organized around the probable therian history
of the internal transcript architecture. Human and mouse provide direct TSS
validation, and marsupials extend the inferred ancestry. Monotreme
cap-enriched transcript data are the decisive remaining evidence for dating
the origin.

## Reproducible outputs

- `scripts/broad_isoform_phylogeny.py`
- `scripts/human_branch_promoter_screen.py`
- `config/human_promoter_intervals.tsv`
- `data/raw/ucsc_multiz470/source_manifest.tsv`
- `results/06_isoform_evolution/broad_isoform_phylogeny/`
- `results/06_isoform_evolution/human_branch_promoter_screen/`

