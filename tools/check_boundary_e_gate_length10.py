#!/usr/bin/env python3
"""Heavy length-ten spot-check for the boundary e-gate identity."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from check_all_disjoint_fusion_placements import crossing, solve_span
from check_boundary_e_gate_crossing import cut_ranks, tensor_e_layer
from check_quantum_braid_fusion_seam import rank
from quantum_braid_fusion_seam import FIELD

ROOT = Path(__file__).resolve().parents[1]
LENGTH = 10
LABELS = 34


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check(output):
    endpoint = crossing(LABELS).reshape(-1)
    tensor = tensor_e_layer(LENGTH, 0, (LENGTH - 1,))
    slices = [
        tensor[degree] for degree in range(tensor.shape[0]) if np.any(tensor[degree])
    ]
    matrix = np.array([piece.reshape(-1) for piece in slices], dtype=np.uint8)
    dimension = rank(matrix)
    augmented = rank(np.concatenate([matrix, endpoint.reshape(1, -1)]))
    assert augmented == dimension, "length-ten e-gate span misses the crossing"
    coefficients = solve_span(matrix, endpoint)
    rebuilt = np.zeros(LABELS**4, dtype=np.uint8)
    for coefficient, piece in zip(coefficients, matrix):
        if coefficient:
            rebuilt ^= FIELD[coefficient, piece]
    if not np.array_equal(rebuilt, endpoint):
        raise AssertionError("length-ten coefficients do not reconstruct the crossing")
    support = int(np.count_nonzero(coefficients))
    values = sorted(set(int(value) for value in coefficients if value))
    assert support == len(slices) == 10
    assert values == [11]
    b0 = list(cut_ranks(slices[0], LABELS))
    assert b0 == [LABELS**2, 1]
    receipt = {
        "field": "F16 with x^4+x+1",
        "status": "completed exact finite spot-check",
        "length": LENGTH,
        "labels": LABELS,
        "site": LENGTH - 1,
        "layer": 0,
        "terms": int(tensor.shape[0]),
        "nonzero_slices": len(slices),
        "span_dimension": dimension,
        "augmented_dimension": augmented,
        "contains_crossing": True,
        "support": support,
        "coefficients": [int(value) for value in coefficients],
        "coefficient_values": values,
        "b0_cut_ranks": b0,
        "scope": (
            "At length 10 (34 labels) the same central first-layer pure e "
            "insertion has ten nonzero degree slices, spans the crossing in "
            "dimension 10, and reconstructs it with coefficient 11 on every "
            "slice. The b=0 slice has cut ranks (1156,1), so the full-first-cut "
            "profile persists. This is a finite spanning check, not a rank-sum "
            "or cost certificate. The full per-slice cut-rank table is omitted "
            "because it is expensive and does not improve the bound."
        ),
        "source_sha256": {
            "tools/check_boundary_e_gate_length10.py": digest(
                ROOT / "tools/check_boundary_e_gate_length10.py"
            ),
            "tools/check_boundary_e_gate_crossing.py": digest(
                ROOT / "tools/check_boundary_e_gate_crossing.py"
            ),
            "tools/quantum_braid_graded_seam.py": digest(
                ROOT / "tools/quantum_braid_graded_seam.py"
            ),
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "evidence/boundary-e-gate-length10.json",
    )
    args = parser.parse_args()
    receipt = check(args.output)
    print(
        json.dumps(
            {
                "length": receipt["length"],
                "labels": receipt["labels"],
                "span_dimension": receipt["span_dimension"],
                "support": receipt["support"],
                "coefficients": receipt["coefficient_values"],
                "b0_cut_ranks": receipt["b0_cut_ranks"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
