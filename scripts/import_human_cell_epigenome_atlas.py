#!/usr/bin/env python3
"""Summarize the adult human single-cell methylome atlas for DNMT3A context."""

from __future__ import annotations

import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import fmean, median


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/human_cell_epigenome_atlas"
CELLS = RAW / "5kCG100k3C_summary.csv.gz"
SUBTYPES = RAW / "subtype_meta.tsv"
OUT = ROOT / "results/05_brain_integration"
MOTOR_CORTEX_CODE = "M1C"


def read_table(path: Path, delimiter: str) -> list[dict[str, str]]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", newline="") as handle:
            return list(csv.DictReader(handle, delimiter=delimiter))
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def quantile(values: list[float], probability: float) -> float:
    """Return a linearly interpolated quantile (NumPy method='linear')."""
    if not values:
        raise ValueError("quantile requires at least one value")
    ordered = sorted(values)
    index = (len(ordered) - 1) * probability
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (index - lower)


def cell_class(annotation: str) -> str:
    """Assign paper annotations to interpretable motor-cortex compartments."""
    if annotation.startswith("Neu Excitatory"):
        return "neuron_excitatory"
    if annotation.startswith("Neu Inhibitory"):
        return "neuron_inhibitory"
    if annotation.startswith("Glia Astrocyte"):
        return "glia_astrocyte"
    if annotation.startswith("Glia Oligodendrocyte Progenitor"):
        return "glia_oligodendrocyte_progenitor"
    if annotation.startswith("Glia Oligodendrocyte"):
        return "glia_oligodendrocyte"
    if annotation.startswith("Hema Myeloid Microglia"):
        return "microglia"
    if annotation.startswith("Endo "):
        return "vascular_endothelial"
    if annotation.startswith("Perivascular "):
        return "vascular_perivascular"
    if annotation.startswith("Fibro "):
        return "fibroblast_leptomeningeal"
    return "other_or_cross_tissue"


def summarize(
    rows: list[dict[str, object]], group_field: str
) -> list[dict[str, object]]:
    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[str(row[group_field])].append(row)

    output = []
    for label, members in sorted(groups.items()):
        mch = [float(row["mCHFrac"]) for row in members]
        mcg = [float(row["mCGFrac"]) for row in members]
        output.append({
            group_field: label,
            "n_cells": len(members),
            "n_donors": len({str(row["Donor"]) for row in members}),
            "n_subtypes": len({str(row["subtype"]) for row in members}),
            "mCH_mean": fmean(mch),
            "mCH_median": median(mch),
            "mCH_q25": quantile(mch, 0.25),
            "mCH_q75": quantile(mch, 0.75),
            "mCG_mean": fmean(mcg),
            "mCG_median": median(mcg),
        })
    return output


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    cells = read_table(CELLS, ",")
    subtype_rows = read_table(SUBTYPES, "\t")
    subtype_map = {row["celltype_L2_both"]: row for row in subtype_rows}

    missing = sorted({row["subtype"] for row in cells} - set(subtype_map))
    if missing:
        raise ValueError(f"{len(missing)} subtype labels lack metadata")

    enriched: list[dict[str, object]] = []
    for row in cells:
        metadata = subtype_map[row["subtype"]]
        enriched.append({
            **row,
            "subtype_abbreviation": metadata["celltype_L2_both_abbr"],
            "subtype_annotation": metadata["celltype_L2_both_annot"],
        })

    tissue_rows = summarize(enriched, "Tissue")
    cortex = [row for row in enriched if row["Tissue"] == MOTOR_CORTEX_CODE]
    for row in cortex:
        row["cell_class"] = cell_class(str(row["subtype_annotation"]))
        row["broad_class"] = (
            "neuronal"
            if str(row["cell_class"]).startswith("neuron_")
            else "non_neuronal"
        )
    class_rows = summarize(cortex, "cell_class")

    donor_groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(
        list
    )
    for row in cortex:
        donor_groups[(str(row["Donor"]), str(row["broad_class"]))].append(row)
    donor_class_rows = []
    for (donor, broad_class), members in sorted(donor_groups.items()):
        base = summarize(members, "broad_class")[0]
        donor_class_rows.append({
            "donor": donor,
            **base,
        })

    subtype_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in cortex:
        subtype_groups[str(row["subtype"])].append(row)
    subtype_summary = []
    for subtype, members in sorted(subtype_groups.items()):
        base = summarize(members, "subtype")[0]
        subtype_summary.append({
            "subtype": subtype,
            "subtype_abbreviation": members[0]["subtype_abbreviation"],
            "subtype_annotation": members[0]["subtype_annotation"],
            "cell_class": members[0]["cell_class"],
            **{key: value for key, value in base.items() if key != "subtype"},
        })
    subtype_summary.sort(
        key=lambda row: (-float(row["mCH_median"]), str(row["subtype"]))
    )

    write_tsv(OUT / "human_cell_atlas_tissue_summary.tsv", tissue_rows)
    write_tsv(
        OUT / "human_cell_atlas_motor_cortex_class_summary.tsv", class_rows
    )
    write_tsv(
        OUT / "human_cell_atlas_motor_cortex_donor_class_summary.tsv",
        donor_class_rows,
    )
    write_tsv(
        OUT / "human_cell_atlas_motor_cortex_subtype_summary.tsv",
        subtype_summary,
    )

    neuronal = [
        row for row in cortex if str(row["cell_class"]).startswith("neuron_")
    ]
    non_neuronal = [
        row for row in cortex if not str(row["cell_class"]).startswith("neuron_")
    ]
    donor_medians = {
        (str(row["donor"]), str(row["broad_class"])): float(row["mCH_median"])
        for row in donor_class_rows
    }
    donor_ratios = {
        donor: donor_medians[(donor, "neuronal")]
        / donor_medians[(donor, "non_neuronal")]
        for donor in sorted({key[0] for key in donor_medians})
    }
    summary = {
        "source_cells": len(cells),
        "tissues": len({row["Tissue"] for row in cells}),
        "major_types": len({row["majortype"] for row in cells}),
        "annotated_subtypes": len(subtype_map),
        "unmapped_subtypes": len(missing),
        "motor_cortex_code": MOTOR_CORTEX_CODE,
        "motor_cortex_cells": len(cortex),
        "motor_cortex_subtypes": len({row["subtype"] for row in cortex}),
        "motor_cortex_neuronal_cells": len(neuronal),
        "motor_cortex_non_neuronal_cells": len(non_neuronal),
        "motor_cortex_neuronal_mCH_median": median(
            float(row["mCHFrac"]) for row in neuronal
        ),
        "motor_cortex_non_neuronal_mCH_median": median(
            float(row["mCHFrac"]) for row in non_neuronal
        ),
        "motor_cortex_donor_neuronal_to_non_neuronal_mCH_median_ratio": (
            donor_ratios
        ),
        "interpretation": (
            "Adult human cell-type context only. This dataset cannot test "
            "cross-species candidate substitutions, developmental mCH "
            "accumulation, or DNMT3A1-versus-DNMT3A2 expression."
        ),
    }
    (OUT / "human_cell_atlas_import_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
