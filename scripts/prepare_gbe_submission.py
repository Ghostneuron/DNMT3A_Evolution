#!/usr/bin/env python3
"""Prepare GBE-oriented manuscript and supplementary text from the master draft."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript"
SOURCE = MANUSCRIPT / "DNMT3A_FULL_MANUSCRIPT_DRAFT.md"
MAIN_OUT = MANUSCRIPT / "GBE_SUBMISSION_DRAFT.md"
SUPP_OUT = MANUSCRIPT / "GBE_SUPPLEMENTARY_TEXT.md"


SIGNIFICANCE = """## Significance statement

DNMT3A produces two major protein forms, but it has been unclear how their
distinct regulatory architectures evolved or how they divide methylation work
during brain development. Across mammals, the variable end of the long form
lies beside a strongly conserved chromatin-interaction region, while an
internal transcript compatible with the short form is broadly distributed and
shows a stage-associated difference in marsupial neocortex. Mouse perturbation data
further indicate that the short form contributes more strongly during
embryonic development, whereas the long form supplies most neuronal mCA after
birth, linking ancient transcript architecture to stage-specific epigenetic
function without invoking human-specific adaptation.
"""

KEYWORDS = (
    "**Keywords:** non-CG methylation; alternative promoter; neuronal "
    "epigenome; purifying selection; transcript architecture; marsupial development"
)

TITLE_BLOCK = """
**Article type:** Article

**Author:** Jie Lu

**Affiliation:** Independent researcher, United States

**Author for correspondence:** Jie Lu; email: ghostneuron@gmail.com
"""

ENDMATTER = """## Acknowledgments

OpenAI Codex was used for code drafting and debugging, workflow organization,
literature-search support, figure generation, and language editing. The author
directed the analyses, evaluated source data and machine-readable outputs,
determined the scientific interpretation, and made all manuscript decisions.
No AI system was treated as an author. The author accepts full responsibility
for the accuracy, originality, and integrity of the work.

## Author contributions

J.L. conceived the study, curated data, performed the analyses, interpreted
the results, prepared the figures and tables, and wrote the manuscript.

## Funding

This work received no external funding.

## Conflict of interest

The author declares no competing interests.
"""

SUPPLEMENTARY_STRUCTURAL_ANALYSIS = """## Exploratory structural analysis of L188H

### Variant selection

At human-reference DNMT3A1 residue 188, leucine was present in 231 of the 250
unique mammalian sequences and histidine was present in 19, comprising 18
primates and the rodent *Micromys minutus*. This distribution establishes that
H188 occurs in multiple sampled species but does not specify the number of
independent L-to-H substitutions. Five natural amino-acid substitutions in the
conserved 164–219 engagement interval could be evaluated on at least one of
two nucleosome-bound templates (Chen et al. 2024; Wapenaar et al. 2024). L188H
was resolved in both 8U5H and 8QZM and
received same-direction EvoEF2 interface scores of −2.55 and −0.29. Histidine
ranked third and sixth, respectively, among the 19 non-reference amino acids
at position 188 (Huang et al. 2020). These comparative scores were used only to select an
exploratory structural case; they are not binding free energies or evidence
of biological function.

### Simulation design and interpretation

The H2A-K119–ubiquitin-G76 isopeptide linkage was represented explicitly.
Amber14 parameters (ff14SB and DNA.OL15) were used with OpenMM 8.4 (Eastman et
al. 2017). Two independently prepared WT–L188H pairs were simulated for 500 ps
per state under periodic PME/NPT conditions at 300 K and 1 atm. One pair used
variant-specific solvent boxes; the other used identical starting solvent
coordinates and box vectors. Twenty frames per system were retained at 25 ps
intervals.

The variant-specific-box pair placed H188 0.642 Å farther from H2A and 0.250 Å
closer to ubiquitin than WT L188. The common-solvent pair instead placed H188
1.331 Å closer to H3, with H2A and ubiquitin differences of −0.102 Å and
−0.031 Å. UDR backbone RMSD and residue-188 RMSF also changed in opposite
directions between preparations. No material distance or mobility effect was
reproduced across both comparisons. These short simulations therefore do not
support altered binding, a fixed partner switch, or a connection to neuronal
methylation, and they were not used to estimate binding free energy.
"""

SUPPLEMENTARY_DATA_TABLES = """## Supplementary data tables

The accompanying workbook `GBE_Supplementary_Data_Tables.xlsx` contains the
following machine-readable tables. Values are reported at the precision of the
validated analysis outputs; blank cells indicate unavailable or inapplicable
fields. Large site-, window-, and gene-level records are provided in
`GBE_Supplementary_Machine_Readable_Data.zip` with a checksum manifest.

- **Table S1. Mammalian DNMT3A sequence inventory and selection-model inclusion.**
- **Table S2. Regional mammalian variability and negative-selection summary.**
- **Table S3. Engagement-region depletion relative to the upstream DNMT3A1 tail.**
- **Table S4. Alignment and selection evidence for recurrent coding candidates.**
- **Table S5. Species-level evidence for downstream-start-compatible DNMT3A products.**
- **Table S6. Clade- and order-level distribution of downstream-start-compatible products.**
- **Table S7. Sensitivity of isoform detection to annotation depth.**
- **Table S8. Human-reference promoter differences relative to aligned local flanks.**
- **Table S9. Dunnart neocortex junction counts and sample metadata.**
- **Table S10. Developmental isoform-specific methylation contrasts.**
- **Table S11. P21 cortical-neuron CpG methylation effect summary.**
- **Table S12. P21 cortical-neuron CpG sample summaries.**
- **Table S13. P21 neuronal EM-seq library quality control and methylation estimates.**
- **Table S14. P21 neuronal mCA and mCH genotype effects.**
- **Table S15. Gene-body mCA contrast summary.**
- **Table S16. Gene-body mCA covariate models.**
- **Table S17. Associations between gene-body mCA and cortex expression effects.**
- **Table S18. MeCP2 gene-set enrichment and methylation-effect analyses.**
- **Table S19. Chromatin-context associations with P21 neuronal mCA effects.**
- **Table S20. DNMT3A occupancy-redistribution models.**
"""

ALT_TEXT = {
    "1": (
        "Three-panel figure showing the DNMT3A1 domain map, high mammalian "
        "variation in residues 1-163, strong constraint in residues 164-219, "
        "and alignment-sensitive support for three candidate sites."
    ),
    "2": (
        "Five-panel figure showing full-length and internal DNMT3A transcript "
        "architectures, their order-level distribution, higher internal-junction "
        "signal at P12 than P20 in dunnart neocortex, annotation-depth effects, "
        "and no promoter enrichment for human-reference differences from a "
        "concordant four-great-ape consensus."
    ),
    "3": (
        "Four-panel figure showing that *Dnmt3a2* knockout has the larger "
        "embryonic brain methylation effect, whereas *Dnmt3a1* knockout has the "
        "larger postnatal effect, including sex-stratified estimates for both "
        "isoforms at both stages."
    ),
    "4": (
        "Four-panel figure showing near-complete loss of P21 neuronal mCA after "
        "*Dnmt3a1* knockout, intermediate retention after N-terminal deletion, "
        "gene-length effects, and enrichment of MeCP2-repressed genes among "
        "large absolute mCA losses."
    ),
    "S1": (
        "Two-panel coefficient plot showing chromatin-context associations with "
        "mCA and DNMT3A occupancy effects; joint models do not isolate a direct "
        "H2AK119ub mediation pathway."
    ),
    "S2": (
        "Two-panel bar plot showing that L188H-minus-WT interface-distance and "
        "mobility changes reverse or differ between molecular-system preparations."
    ),
}

JOURNAL_ABBREVIATIONS = {
    "Bioinformatics": "Bioinformatics",
    "Cell": "Cell",
    "Communications Biology": "Commun Biol",
    "Cytogenetic and Genome Research": "Cytogenet Genome Res",
    "Development": "Development",
    "EMBO Journal": "EMBO J",
    "EMBO Reports": "EMBO Rep",
    "Genome Biology and Evolution": "Genome Biol Evol",
    "Journal of Biological Chemistry": "J Biol Chem",
    "Methods in Molecular Biology": "Methods Mol Biol",
    "Molecular Biology and Evolution": "Mol Biol Evol",
    "Nature Ecology & Evolution": "Nat Ecol Evol",
    "Nature Communications": "Nat Commun",
    "Nature Genetics": "Nat Genet",
    "Nature Neuroscience": "Nat Neurosci",
    "PLoS Computational Biology": "PLoS Comput Biol",
    "PLoS Genetics": "PLoS Genet",
    "Scientific Data": "Sci Data",
    "Scientific Reports": "Sci Rep",
}


def format_reference(entry: str) -> str:
    pattern = re.compile(
        r"^(?P<authors>.+?)\. (?P<title>.+?)\. \*(?P<journal>.+?)\*\. "
        r"(?P<year>\d{4});(?P<citation>.+?)\. "
        r"\[doi:(?P<doi>[^\]]+)\]\(.+\)\.?$"
    )
    match = pattern.match(entry)
    if not match:
        return entry
    fields = match.groupdict()
    authors = fields["authors"]
    if "et al." not in authors and len(authors.split(", ")) > 5:
        authors = authors.split(", ", 1)[0] + ", et al."
    journal = JOURNAL_ABBREVIATIONS.get(fields["journal"], fields["journal"])
    authors = authors.rstrip(".")
    return (
        f"{authors}. {fields['year']}. {fields['title']}. "
        f"*{journal}* {fields['citation']}. doi:{fields['doi']}."
    )


def sort_references(section: str) -> str:
    entries = []
    for line in section.strip().splitlines():
        match = re.match(r"^\d+\.\s+(.*)$", line.strip())
        if match:
            entries.append(match.group(1))
    if not entries:
        raise ValueError("No numbered references found")
    entries = [format_reference(entry) for entry in entries]
    entries.sort(key=lambda value: value.casefold())
    return "\n\n".join(entries)


def split_legends(legend_section: str) -> tuple[list[str], list[str]]:
    blocks = [block.strip() for block in legend_section.strip().split("\n\n") if block.strip()]
    main, supplementary = [], []
    for block in blocks:
        match = re.match(r"\*\*Figure (S?\d+)\.(.*?)\*\*(.*)", block, flags=re.S)
        if not match:
            continue
        number, title, body = match.groups()
        formatted = f"**Fig. {number}.{title}**{body}\n\n**Alt text:** {ALT_TEXT[number]}"
        (supplementary if number.startswith("S") else main).append(formatted)
    return main, supplementary


def main() -> None:
    text = SOURCE.read_text()
    title, rest = text.split("\n", 1)
    abstract_at = rest.index("## Abstract")
    intro_at = rest.index("## Introduction")
    legends_at = rest.index("## Figure legends")
    refs_at = rest.index("## References")

    abstract = rest[abstract_at:intro_at].strip()
    body = rest[intro_at:legends_at].strip()
    legend_text = rest[legends_at + len("## Figure legends"):refs_at]
    reference_text = rest[refs_at + len("## References"):]
    main_legends, supplementary_legends = split_legends(legend_text)

    body = body.replace("## Methods", "## Materials and Methods", 1)
    body = body.replace("### Data and code availability", "## Data availability", 1)

    prepared = "\n\n".join(
        [
            title,
            TITLE_BLOCK.strip(),
            abstract,
            KEYWORDS,
            SIGNIFICANCE.strip(),
            body,
            ENDMATTER.strip(),
            "## Figure legends\n\n" + "\n\n".join(main_legends),
            "## References\n\n" + sort_references(reference_text),
        ]
    )
    MAIN_OUT.write_text(prepared + "\n")

    supplementary = "\n\n".join(
        [
            f"# Supplementary information for: {title.removeprefix('# ')}",
            SUPPLEMENTARY_DATA_TABLES.strip(),
            SUPPLEMENTARY_STRUCTURAL_ANALYSIS.strip(),
            "## Supplementary figure legends\n\n" + "\n\n".join(supplementary_legends),
            "## Supplementary figure files\n\n"
            "- `figures/manuscript/Figure_S1_chromatin_context.pdf`\n"
            "- `figures/manuscript/Figure_S2_L188H_preparation_sensitivity.pdf`\n\n"
            "## Supplementary data files\n\n"
            "- `GBE_Supplementary_Data_Tables.xlsx`\n"
            "- `GBE_Supplementary_Machine_Readable_Data.zip`",
        ]
    )
    SUPP_OUT.write_text(supplementary + "\n")
    print(MAIN_OUT)
    print(SUPP_OUT)


if __name__ == "__main__":
    main()
