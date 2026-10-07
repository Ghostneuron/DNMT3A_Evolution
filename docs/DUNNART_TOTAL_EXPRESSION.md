# Dunnart total DNMT3A-like expression

## Question

P12 has more exact internal-junction support than P20. This analysis asks how
much of that difference can be explained by a general increase in total
DNMT3A expression.

## Independent source-study input

The analysis uses the processed RSEM/Trinity gene-count matrix deposited with
GEO series `GSE161274`, rather than deriving total expression from the
competitive RefSeq Salmon model. The original GSE Trinity assembly and its
transcript-to-gene map were downloaded with the matrix and audited together.
This matters because Trinity identifiers from a later multi-tissue assembly
are not interchangeable with identifiers in the source-study matrix.

DNMT3A-like source transcripts were identified by exact sequence: a candidate
must contain at least 100 distinct 31-mers also present in one or more dunnart
RefSeq DNMT3A models. Forty-four original-assembly transcripts map to 16
source-matrix gene components at this threshold. Trinity fragmentation means
these are operational DNMT3A-like components, not evidence for 16 biological
DNMT3A loci.

Expected gene counts were summed across those components and divided by the
total assigned gene count in each source-study library to obtain CPM.

## Result

| stage | specimen | DNMT3A-like CPM | internal junction pairs/M per gene CPM |
|---|---:|---:|---:|
| P12 | 436 | 55.4104 | 0.010973 |
| P12 | 451 | 58.6110 | 0.012853 |
| P12 | 458 | 64.9498 | 0.013281 |
| P20 | 666 | 39.1864 | 0.005469 |
| P20 | 746 | 38.9173 | 0.002999 |
| P20 | 758 | 42.1687 | 0.006585 |

Pooled DNMT3A-like CPM is 59.806 at P12 and 40.475 at P20, a 1.48-fold
difference. Every P12 animal exceeds every P20 animal; the exact one-sided
animal-label permutation P value is 0.05.

The exact internal-junction rate differs by 3.40-fold. After normalizing that
rate to the source-study DNMT3A-like CPM, the pooled P12/P20 fold is 2.30.
Again, every P12 animal exceeds every P20 animal and the exact permutation
P value is 0.05. Thus higher total expression explains part, but not all, of
the internal-junction contrast.

## Candidate-threshold sensitivity

The total-expression result is stable when the exact-sequence requirement is
made progressively stricter:

| minimum distinct exact target 31-mers | candidate gene components | P12/P20 CPM fold |
|---:|---:|---:|
| 100 | 16 | 1.478 |
| 250 | 10 | 1.479 |
| 500 | 7 | 1.478 |
| 1000 | 4 | 1.614 |

All four thresholds retain complete P12-over-P20 separation and an exact
one-sided P value of 0.05.

## Interpretation boundary

This result supports two statements:

1. total DNMT3A-like expression is higher at P12;
2. internal-junction support rises more than total DNMT3A-like expression.

It does not establish a discrete internal-over-full-length switch. The direct
internal/full-length first-junction ratio remains overlapping and
non-significant. CPM is a transparent library-size normalization of the
deposited expected-count matrix, not a new edgeR differential-expression
analysis. With three animals per stage, all exact P values remain minimally
resolved. Missing sex and explicit batch metadata also remain limitations.

## Reproduction

```bash
python scripts/dunnart_source_dnmt3a_expression.py
```

Outputs:

- `dunnart_source_dnmt3a_candidate_audit.tsv`;
- `dunnart_source_dnmt3a_expression.tsv`;
- `dunnart_source_dnmt3a_expression_summary.json`.
