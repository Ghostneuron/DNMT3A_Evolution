import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_published_loss_controls_are_directionally_destabilizing():
    path = ROOT / "structural_modeling/8QZM/results/summary/combined_evoef_mutation_screen.tsv"
    rows = read_tsv(path)
    controls = [row for row in rows if row["class"] == "positive_control"]
    assert controls
    for row in controls:
        for structure in ("8U5H", "8QZM"):
            value = row[f"{structure}_delta_binding_score"]
            if value:
                assert float(value) > 0


def test_natural_variant_panel_is_represented():
    path = ROOT / "structural_modeling/8QZM/results/summary/combined_evoef_mutation_screen.tsv"
    observed = {row["mutation"] for row in read_tsv(path) if row["class"] == "natural"}
    assert {"S166A", "L174Q", "M185V", "P186Q", "L188H"} <= observed
