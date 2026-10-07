# DNMT3A evolutionary analysis plan

## Central question

How is evolutionary constraint partitioned among the isoform-specific
DNMT3A1 regulatory tail, its experimentally supported chromatin-engagement
region, and the conserved domains of mammalian DNMT3A? Neuronal methylation
and developmental expression provide functional context and future tests,
not a prespecified association claim for the current dataset.

## Phase 1: mammalian coding evolution

Primary dataset: one annotation-validated, full-length DNMT3A1-like ortholog per
species from the NCBI human-ortholog package. Because the package spans
sarcopterygian vertebrates, the primary reporting units will be (i) mammals and
(ii) the broader tetrapod/vertebrate panel. A trusted taxonomy table and species
tree are required before clade-specific tests.

Analyses:

- representative and CDS-translation audit;
- MAFFT protein guide with exact codon projection, followed by primary MACSE
  codon-aware refinement and pre-refinement sensitivity audit;
- alignment occupancy, frameshift, ambiguity, stop-codon, and long-branch QC;
- IQ-TREE nucleotide topology from the codon alignment, constrained only at
  the well-established mammalian-order level and compared with an unconstrained
  diagnostic tree;
- FUBAR pervasive-site, MEME episodic-site, and aBSREL branch tests;
- domain-level substitution/selection summaries;
- subregion-level comparison of residues 1–163, the experimentally supported
  nucleosome/H2AK119ub-engagement region at 164–219, and residues 220–277;
- mapping to human residues, known functional surfaces, and clinical variants.

Gate to interpretation: a candidate site must have reliable human-coordinate
mapping, adequate column occupancy, no obvious annotation/alignment artifact,
and support that survives at least one alignment or taxon-sampling sensitivity
analysis. Multiple-testing correction and the distinct null hypotheses of each
selection method must be reported.

Occupancy filtering must be calculated within the taxon set used by the
selection model. Filtering the 475-vertebrate alignment before extracting
mammals removed eight mammal-high-occupancy human sites and initially hid S97.
The mammalian production sensitivity therefore refilters raw MACSE and MAFFT
codon alignments across the 250 unique mammals before FUBAR/MEME. See
`docs/S97_HOMOLOGY_AUDIT.md`.

For low-complexity N-terminal candidates, this gate also requires
annotation-aware positional homology. Perfect mapping of a predicted RNA to
its source assembly does not suffice when a predicted first exon conflicts
with close-relative protein architecture. The Puma/wombat audit excludes
Puma's P5/S6/E18 states and leaves the corresponding wombat states unresolved;
see `docs/CANDIDATE_ANNOTATION_VALIDATION.md`.

## Phase 2: gene-family and isoform history

- add DNMT3B and non-mammalian DNMT3 sequences as outgroups;
- reconcile gene and species trees and validate ambiguous copies by synteny;
- analyze teleost `dnmt3aa`/`dnmt3ab` separately;
- reconstruct DNMT3A1/DNMT3A2 transcription-start-region evolution using
  whole-genome alignments and cross-species transcript evidence.

## Phase 3: brain-development integration

Preferred phenotypes are mechanistically proximal:

- neuronal or brain mCH/mCA abundance;
- developmental timing or slope of mCH accumulation;
- DNMT3A1 and DNMT3A2 developmental expression trajectories;
- neuronal maturation timing;
- MeCP2 abundance or mCA-reading context where comparable.

The template at `data/traits/brain_phenotypes.tsv` records values, units,
tissue/cell type, developmental stage, assay, accession, and citation. Values
from unmatched tissues, life stages, or assays must not be combined without an
explicit harmonization model.

Association tests will use phylogenetic regression or phylogenetic mixed
models, report effective species coverage, and include sensitivity analyses for
clade, body mass, lifespan, and data-source effects where relevant. Brain size,
encephalization, flight, aquatic adaptation, and longevity are secondary
exploratory traits—not substitutes for neuroepigenetic measurements.

## Phase 4: functional prioritization

Prioritize evolutionary candidates using convergent evidence:

- domain/structure location;
- biochemical or chromatin-binding relevance;
- overlap with pathogenic or experimentally tested residues;
- association with developmental expression/methylation;
- recurrence or convergent substitution after phylogenetic correction.

Selection scans alone will not be described as evidence of altered brain
development.

### Staged production execution

FUBAR is run first on the validated mammalian tree. MEME is enabled only after
site-coordinate and FUBAR output QC. An indiscriminate all-branch aBSREL run is
not a default production step because the tree contains near-zero branches,
exact sequence-identity groups, composition-test failures, and
constraint-enforced nodes. aBSREL foreground branches must be defined in an
auditable manifest after these cases are excluded or grouped.

The first targeted aBSREL follow-up uses the two alignment-robust internal
branches with reconstructed G34 changes: the sampled *Rhinolophus* and
lorisiform stems. The *Carlito syrichta* terminal change is excluded because
its G34 residue is not retained under PRANK. Foreground selection, branch
length, support, constraint status, composition status, exact-sequence
deduplication, and PRANK concordance are recorded in
`results/04_selection/mammals/g34_absrel/foreground_manifest.tsv`. Because the
foregrounds were chosen after the G34 MEME result, aBSREL is interpreted only
as post hoc branch-wide concordance and is Holm-corrected across the two
foreground tests.

The first FUBAR sensitivity analysis repeats the fit after exact-sequence
deduplication and exclusion of the 23 taxa failing IQ-TREE's nucleotide
composition chi-square screen. This exclusion is diagnostic, not a claim that
the flagged taxa are biologically invalid.

The second sensitivity analysis repeats FUBAR on the pre-refinement MAFFT
codon projection after occupancy is recalculated within the 250 mammals.
Candidate comparison is performed by human DNMT3A1 residue through an
independently generated coordinate crosswalk, not by assuming that raw
alignment columns are interchangeable.

The third alignment sensitivity uses PRANK codon mode with the rooted mammalian
guide topology. It begins from ungapped CDS records, independently maps human
coordinates, applies the same occupancy threshold, and repeats FUBAR plus
targeted MEME. Terminal residue assignments that become gaps under PRANK are
not retained as robust lineage changes.

The initial MEME execution is restricted to the union of FUBAR candidates
(filtered sites 9, 30, and 107). This is a hypothesis-follow-up analysis and
not a substitute for a future alignment-wide MEME scan. Each candidate is run
as a separate job and output because HyPhy 2.5.79's list matcher can also match
shorter numeric prefixes of multi-digit sites. Only the explicitly requested
row is extracted from each result.

The alignment-wide MEME scan has now been completed as a separate production
analysis with no site filter. Discovery statistics are corrected across all
901 codons. Sites passing 5% BH correction are mapped independently into
MAFFT and PRANK, tested one site per job, and Holm-corrected across the
discovery candidates within each alternative alignment. A primary-alignment
discovery is not called alignment-robust unless its human-mapped site remains
significant in both alternative alignments.
