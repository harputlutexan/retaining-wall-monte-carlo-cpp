"""Compare the generated summary and per-design results to tracked references."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "sample_results"

expected_json = json.loads((REFERENCE / "retaining_wall_optimization_cpp.json").read_text())
actual_json = json.loads((ROOT / "retaining_wall_optimization_cpp.json").read_text())
assert expected_json.keys() == actual_json.keys(), "JSON summary field mismatch"
for key in expected_json:
    expected = expected_json[key]
    actual = actual_json[key]
    if isinstance(expected, float):
        assert math.isclose(expected, actual, abs_tol=0.011), (key, expected, actual)
    else:
        assert expected == actual, (key, expected, actual)

with (REFERENCE / "all_design_results_cpp.csv").open(newline="") as file:
    expected_rows = list(csv.DictReader(file))
with (ROOT / "all_design_results_cpp.csv").open(newline="") as file:
    actual_rows = list(csv.DictReader(file))
assert len(expected_rows) == len(actual_rows) == 2592, "Unexpected design row count"
for expected_row, actual_row in zip(expected_rows, actual_rows):
    assert expected_row.keys() == actual_row.keys(), "Design CSV columns differ"
    for field in expected_row:
        a, b = expected_row[field], actual_row[field]
        if not a and not b:
            continue
        if field in {"design_index", "success_count", "baseline", "robust"}:
            assert a == b, (expected_row["design_index"], field, a, b)
        else:
            assert math.isclose(float(a), float(b), rel_tol=1e-8, abs_tol=1e-7), (
                expected_row["design_index"], field, a, b
            )
print("PASS: all 10 JSON fields and all 2,592 design records match references")
