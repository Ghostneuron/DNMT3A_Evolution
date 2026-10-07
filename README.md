# DNMT3A molecular evolution and developmental isoform specialization

This repository contains analysis code, configuration, validation tests,
documentation, figure files, and compact supplementary tables supporting the
manuscript *Regional constraint and developmental isoform specialization in
mammalian DNMT3A*.

The study integrates mammalian DNMT3A coding-sequence evolution, alternative
transcript architecture, developmental methylation, neuronal non-CG
methylation, gene expression, and chromatin context. The principal conclusions
concern regional constraint in the DNMT3A1 N terminus and developmental
allocation of methylation function between DNMT3A isoforms. The repository
does not treat individual recurrent substitutions as established adaptive or
causal variants.

## Repository contents

```text
config/                 project settings, sample tables, and track manifests
data/supplementary/     compact workbook and machine-readable manifest
docs/                   analysis and quality-control documentation
environment/            Python dependency snapshot
figures/manuscript/     manuscript figures in PDF, PNG, and SVG formats
scripts/                analysis and submission-preparation scripts
tests/                  validation tests
```

The sequence alignments, phylogenetic tree, promoter-sequence sets, and eight
larger site-, window-, and gene-level data tables have been prepared for a
companion Zenodo deposit. They are not duplicated in this GitHub repository.
The Zenodo DOI will be added before publication. All underlying biological
inputs were obtained from public resources, including GEO accessions
GSE164265, GSE295720, and GSE161274.

## Environment and external tools

The Python package snapshot is recorded in `environment/requirements.txt`.
The main sequence workflow also requires MAFFT, MACSE 2, IQ-TREE 2, and HyPhy
2.5 or later. Assay-specific tool versions are recorded in
`config/mch_tool_versions.tsv`.

To inspect the primary pipeline configuration:

```bash
python3 scripts/dnmt3a_pipeline.py doctor
```

Data-dependent stages require the public source datasets and the processed
Zenodo tables described in `data/supplementary/README.md`.

## Data availability statement

Analysis code, versioned configuration files, documentation, validation tests,
manuscript figures, and Supplementary Data Tables S1-S20 are available from
[GitHub](https://github.com/Ghostneuron/DNMT3A_Evolution). The sequence
alignments, phylogenetic tree, promoter-sequence sets, and eight larger derived
tables have been prepared for archival deposition in Zenodo. The DOI will be
inserted before publication: `ZENODO_DOI_TO_BE_ADDED_BEFORE_PUBLICATION`.
Source genomic and functional-genomic data remain available from the public
repositories and accessions documented in the manuscript and project
manifests.

## Citation and licensing

Citation metadata are provided in `CITATION.cff`. No public-use license has
yet been assigned; see `LICENSE_REVIEW_REQUIRED.md` before the version 1.0.0
release.
