#!/usr/bin/env python3
"""Audit measurable GBE Article requirements in the prepared manuscript."""

from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "manuscript" / "GBE_SUBMISSION_DRAFT.md"
SUPPLEMENT_PATH = ROOT / "manuscript" / "GBE_SUPPLEMENTARY_TEXT.md"
PACKAGE = ROOT / "submission" / "GBE_initial_submission_2026-10-06"


def words(text: str) -> int:
    return len(re.findall(r"\b[\w’'-]+\b", text))


def between(text: str, start: str, end: str) -> str:
    return text.split(start, 1)[1].split(end, 1)[0]


def main() -> None:
    text = PATH.read_text()
    supplement = SUPPLEMENT_PATH.read_text()
    title = text.splitlines()[0].removeprefix("# ")
    abstract = between(text, "## Abstract", "**Keywords:**")
    significance = between(text, "## Significance statement", "## Introduction")
    text_without_refs = text.split("## References", 1)[0]
    keywords_line = next(line for line in text.splitlines() if line.startswith("**Keywords:**"))
    keywords = [item.strip() for item in keywords_line.split(":**", 1)[1].split(";")]
    reference_lines = [
        line for line in text.split("## References", 1)[1].splitlines()
        if re.search(r"\. \d{4}\.", line)
    ]
    reference_keys = [line.casefold() for line in reference_lines]
    headings = [
        "## Abstract",
        "## Significance statement",
        "## Introduction",
        "## Results",
        "## Discussion",
        "## Materials and Methods",
        "## Data availability",
        "## Figure legends",
        "## References",
    ]
    positions = [text.index(heading) for heading in headings]

    audit = {
        "title_characters": len(title),
        "title_within_150": len(title) <= 150,
        "abstract_words": words(abstract),
        "abstract_within_250": words(abstract) <= 250,
        "significance_words": words(significance),
        "significance_within_150": words(significance) <= 150,
        "significance_sentences": len(re.findall(r"[.!?](?:\s|$)", significance.strip())),
        "keyword_count": len(keywords),
        "keywords_within_6": len(keywords) <= 6,
        "text_words_excluding_references": words(text_without_refs),
        "text_within_10000": words(text_without_refs) <= 10000,
        "main_figure_legends": len(re.findall(r"^\*\*Fig\. [1-4]\.", text, flags=re.M)),
        "main_alt_text_entries": len(re.findall(r"^\*\*Alt text:\*\*", text, flags=re.M)),
        "supplementary_legends_absent": not re.search(r"^\*\*Fig\. S", text, flags=re.M),
        "reference_count": len(reference_lines),
        "references_alphabetical": reference_keys == sorted(reference_keys),
        "required_section_order": positions == sorted(positions),
        "em_dash_absent": "—" not in text,
        "supplementary_data_table_entries": len(
            re.findall(r"^- \*\*Table S\d+\.", supplement, flags=re.M)
        ),
        "supplementary_data_table_range_cited": "Supplementary Data Tables S1–S20" in text,
        "data_table_workbook_present": (PACKAGE / "GBE_Supplementary_Data_Tables.xlsx").exists(),
        "machine_readable_archive_present": (
            PACKAGE / "GBE_Supplementary_Machine_Readable_Data.zip"
        ).exists(),
    }
    audit["all_measurable_checks_pass"] = all(
        value for key, value in audit.items()
        if key.endswith(("_150", "_250", "_6", "_10000", "_absent", "_alphabetical", "_order"))
    ) and all(
        (
            audit["main_figure_legends"] == 4,
            audit["main_alt_text_entries"] == 4,
            audit["supplementary_data_table_entries"] == 20,
            audit["supplementary_data_table_range_cited"],
            audit["data_table_workbook_present"],
            audit["machine_readable_archive_present"],
        )
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
