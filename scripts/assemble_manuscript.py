#!/usr/bin/env python3
"""Assemble the DNMT3A manuscript from audited section drafts."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript"


def body_after_heading(text: str, heading: str) -> str:
    if not text.startswith(heading):
        raise ValueError(f"expected leading heading {heading!r}")
    return text[len(heading) :].strip()


def main() -> None:
    core = (MANUSCRIPT / "DNMT3A_RESULTS_DISCUSSION_DRAFT.md").read_text()
    introduction = body_after_heading(
        (MANUSCRIPT / "INTRODUCTION_DRAFT.md").read_text(), "# Introduction"
    )
    methods = body_after_heading(
        (MANUSCRIPT / "METHODS_DRAFT.md").read_text(), "# Methods"
    )
    if methods.startswith("## "):
        methods = "### " + methods[3:]
    methods = methods.replace("\n## ", "\n### ")
    references = (MANUSCRIPT / "WORKING_REFERENCES.md").read_text()

    results_at = core.index("\n## Results")
    legends_at = core.index("\n## Figure legends")
    front = core[:results_at].strip()
    results_discussion = core[results_at:legends_at].strip()
    legends = core[legends_at:].strip()

    reference_start = references.index("\n1. ") + 1
    reference_end_marker = "\nThe entries above"
    reference_end = references.find(reference_end_marker, reference_start)
    if reference_end == -1:
        reference_end = len(references)
    reference_list = references[reference_start:reference_end].strip()

    assembled = "\n\n".join(
        [
            front,
            "## Introduction\n\n" + introduction,
            results_discussion,
            "## Methods\n\n" + methods,
            legends,
            "## References\n\n" + reference_list,
        ]
    )
    output = MANUSCRIPT / "DNMT3A_FULL_MANUSCRIPT_DRAFT.md"
    output.write_text(assembled + "\n")
    print(output)


if __name__ == "__main__":
    main()
