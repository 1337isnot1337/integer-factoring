#!/usr/bin/env python3
"""Exact boundary e-gate identity for the graded crossing."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from check_all_disjoint_fusion_placements import crossing, solve_span
from check_quantum_braid_fusion_seam import rank
from paired_path_crossing import power
from quantum_braid_fusion_seam import FIELD, Q
from quantum_braid_graded_seam import _apply_e, _layout

ROOT = Path(__file__).resolve().parents[1]
LENGTHS = (2, 4, 6, 8)


def tensor_e_layer(length, layer, sites):
    labels, states_count, boundary, routes = _layout(length)
    states = np.zeros((states_count, labels * labels, 1), dtype=np.uint8)
    states[boundary, np.arange(labels * labels), 0] = 1
    for diagonal in range(2 * length - 1):
        if diagonal == layer:
            for site in sites:
                states = _apply_e(states, routes[site])
            continue
        for row in range(length):
            column = diagonal - row
            if not 0 <= column < length:
                continue
            action = _apply_e(states, routes[length + row - 1 - column])
            if row == column:
                following = np.zeros(
                    states.shape[:-1] + (states.shape[-1] + 1,), dtype=np.uint8
                )
                following[..., :-1] = FIELD[Q, states]
                following[..., 1:] ^= action
                states = following
            else:
                states = FIELD[Q, states] ^ action
    return FIELD[
        power(Q, length - length * length // 2),
        states[boundary].transpose(2, 1, 0).reshape(-1, labels, labels, labels, labels),
    ]


def cut_ranks(tensor, labels):
    first = rank(tensor.reshape(labels * labels, labels * labels))
    second = rank(tensor.transpose(0, 2, 1, 3).reshape(labels * labels, labels * labels))
    return first, second


def identity_record(length):
    labels, _, _, _ = _layout(length)
    endpoint = crossing(labels).reshape(-1)
    tensor = tensor_e_layer(length, 0, (length - 1,))
    slices = [tensor[degree] for degree in range(tensor.shape[0]) if np.any(tensor[degree])]
    matrix = np.array([piece.reshape(-1) for piece in slices], dtype=np.uint8)
    coefficients = solve_span(matrix, endpoint)
    rebuilt = np.zeros(labels**4, dtype=np.uint8)
    for coefficient, piece in zip(coefficients, matrix):
        if coefficient:
            rebuilt ^= FIELD[coefficient, piece]
    if not np.array_equal(rebuilt, endpoint):
        raise AssertionError("boundary e-gate identity does not reconstruct the crossing")
    support = int(np.count_nonzero(coefficients))
    values = sorted(set(int(value) for value in coefficients if value))
    ranks = [list(cut_ranks(piece, labels)) for piece in slices]
    return {
        "length": length,
        "labels": labels,
        "site": length - 1,
        "layer": 0,
        "terms": int(tensor.shape[0]),
        "nonzero_slices": len(slices),
        "support": support,
        "coefficients": [int(value) for value in coefficients],
        "coefficient_values": values,
        "cut_ranks": ranks,
        "matches_crossing": True,
    }


def site_census(length):
    labels, _, _, _ = _layout(length)
    endpoint = crossing(labels).reshape(-1)
    records = []
    for site in range(2 * length - 2):
        tensor = tensor_e_layer(length, 0, (site,))
        slices = [
            tensor[degree]
            for degree in range(tensor.shape[0])
            if np.any(tensor[degree])
        ]
        matrix = np.array([piece.reshape(-1) for piece in slices], dtype=np.uint8)
        dimension = rank(matrix)
        contains = rank(np.concatenate([matrix, endpoint.reshape(1, -1)])) == dimension
        records.append(
            {
                "site": site,
                "nonzero_slices": len(slices),
                "span_dimension": dimension,
                "contains_crossing": bool(contains),
            }
        )
    return records


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check(output):
    records = [identity_record(length) for length in LENGTHS]
    expected = {2: (1, 13), 4: (4, 7), 6: (6, 3), 8: (8, 13)}
    for record in records:
        support, coefficient = expected[record["length"]]
        assert record["support"] == support
        assert record["coefficient_values"] == [coefficient]
        assert record["matches_crossing"]
    census = site_census(8)
    central = [entry for entry in census if entry["contains_crossing"]]
    assert len(central) == 1 and central[0]["site"] == 7
    assert all(entry["span_dimension"] == 7 for entry in census if entry["site"] != 7)
    receipt = {
        "field": "F16 with x^4+x+1",
        "status": "completed exact finite checks",
        "records": records,
        "site_census_length8": census,
        "scope": (
            "At lengths 2,4,6,8 the graded crossing is an exact scalar multiple "
            "of the sum of the degree slices of one pure e insertion at the "
            "central site in the first layer. At length8 no other site spans "
            "the crossing. This is a finite spanning identity, not a rank-sum "
            "or cost certificate; every length8 slice has full first-cut rank."
        ),
        "source_sha256": {
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
        default=ROOT / "evidence/boundary-e-gate-crossing.json",
    )
    args = parser.parse_args()
    receipt = check(args.output)
    print(
        json.dumps(
            [
                {
                    "length": record["length"],
                    "labels": record["labels"],
                    "support": record["support"],
                    "coefficients": record["coefficient_values"],
                }
                for record in receipt["records"]
            ],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
