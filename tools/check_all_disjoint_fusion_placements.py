#!/usr/bin/env python3
"""Exhaustive disjoint f3 placements and exact graded crossing spans."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from check_quantum_braid_fusion_seam import rank
from paired_path_crossing import power
from quantum_braid_fusion_seam import FIELD, Q
from quantum_braid_graded_seam import (
    _layout,
    all_disjoint_projector_sites,
    coefficient_tensor,
    coefficient_tensor_sites,
    free_counts,
    projector_sites,
)

ROOT = Path(__file__).resolve().parents[1]
INVERSE = np.array([0] + [power(a, 14) for a in range(1, 16)], dtype=np.uint8)


def crossing(labels):
    tensor = np.zeros((labels, labels, labels, labels), dtype=np.uint8)
    for p in range(labels):
        for q in range(labels):
            tensor[p, q, q, p] = 1
    return tensor


def identity(labels):
    tensor = np.zeros((labels, labels, labels, labels), dtype=np.uint8)
    for p in range(labels):
        for q in range(labels):
            tensor[p, q, p, q] = 1
    return tensor


def add_piece(pieces, labels, tensor, placement, degree, orientation):
    if not np.any(tensor):
        return
    pieces.append((tensor.reshape(-1), placement, degree, orientation))
    pieces.append(
        (
            tensor.transpose(0, 2, 1, 3).reshape(-1),
            placement,
            degree,
            "column" if orientation == "row" else "row",
        )
    )


def solve_span(matrix, target):
    augmented = np.concatenate([matrix.T, target.reshape(-1, 1)], axis=1)
    rows = augmented.copy()
    pivots = []
    pivot_row = 0
    for column in range(matrix.shape[0]):
        candidates = np.flatnonzero(rows[pivot_row:, column])
        if not len(candidates):
            continue
        pivot = pivot_row + candidates[0]
        rows[[pivot_row, pivot]] = rows[[pivot, pivot_row]]
        rows[pivot_row] = FIELD[INVERSE[rows[pivot_row, column]], rows[pivot_row]]
        for row in range(len(rows)):
            if row != pivot_row and rows[row, column]:
                rows[row] ^= FIELD[rows[row, column], rows[pivot_row]]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    for row in range(pivot_row, len(rows)):
        if not rows[row, :-1].any() and rows[row, -1]:
            raise ValueError("crossing is not in the span")
    coefficients = np.zeros(matrix.shape[0], dtype=np.uint8)
    for row, column in enumerate(pivots):
        coefficients[column] = rows[row, -1]
    return coefficients


def span_record(length, placements):
    labels, _, _, _ = _layout(length)
    pieces = []
    for placement in placements:
        tensor = coefficient_tensor_sites(length, placement)
        for degree in range(length + 1):
            add_piece(pieces, labels, tensor[degree], placement, degree, "row")
    add_piece(pieces, labels, identity(labels), None, None, "identity")
    matrix = np.array([piece[0] for piece in pieces], dtype=np.uint8)
    dimension = rank(matrix)
    endpoint = crossing(labels)
    augmented = np.concatenate([matrix, endpoint.reshape(1, -1)], axis=0)
    membership = rank(augmented) == dimension
    recombination = []
    if membership:
        coefficients = solve_span(matrix, endpoint)
        for index, coefficient in enumerate(coefficients):
            if coefficient:
                _, placement, degree, orientation = pieces[index]
                recombination.append(
                    {
                        "placement": None if placement is None else list(placement),
                        "degree": degree,
                        "orientation": orientation,
                        "coefficient": int(coefficient),
                    }
                )
        rebuilt = np.zeros(matrix.shape[1], dtype=np.uint8)
        for index, coefficient in enumerate(coefficients):
            if coefficient:
                rebuilt ^= FIELD[coefficient, matrix[index]]
        assert np.array_equal(rebuilt, endpoint.reshape(-1))
    return {
        "length": length,
        "labels": labels,
        "placements": len(placements),
        "nonzero_pieces": len(pieces),
        "span_dimension": dimension,
        "crossing_in_span": membership,
        "recombination_support": len(recombination),
        "recombination": recombination,
    }


def prefix_record(length):
    placements = []
    for left in free_counts(length):
        for right in free_counts(length):
            placement = projector_sites(length, left, right)
            assert np.array_equal(
                coefficient_tensor_sites(length, placement),
                coefficient_tensor(length, left, right),
            )
            placements.append(placement)
    return span_record(length, tuple(placements))


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check(output):
    records = []
    for length in (4, 6):
        all_placements = all_disjoint_projector_sites(length)
        all_record = span_record(length, all_placements)
        prefix = prefix_record(length)
        expected = {
            4: {"placements": 13, "prefix_dimension": 13},
            6: {"placements": 60, "prefix_dimension": 80},
        }[length]
        assert all_record["placements"] == expected["placements"]
        assert prefix["span_dimension"] == expected["prefix_dimension"]
        assert all_record["crossing_in_span"]
        assert prefix["crossing_in_span"] == (length == 4)
        records.append({"all_placements": all_record, "prefix_control": prefix})

    receipt = {
        "field": "F16 with x^4+x+1",
        "status": "completed exact finite checks",
        "records": records,
        "scope": (
            "All disjoint three-strand f3 placements on 2*length central sites, "
            "both cut orientations and the two identity pairings, over F16. "
            "The span membership and recombination are exact. No rank-sum "
            "optimum, uniform family, or stronger factoring score is claimed."
        ),
        "source_sha256": {
            "research/invertible-fusion-graded-seams.md": digest(
                ROOT / "research/invertible-fusion-graded-seams.md"
            ),
            "tools/check_all_disjoint_fusion_placements.py": digest(
                ROOT / "tools/check_all_disjoint_fusion_placements.py"
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
        default=ROOT / "evidence/all-disjoint-fusion-placements.json",
    )
    args = parser.parse_args()
    receipt = check(args.output)
    print(
        json.dumps(
            [
                {
                    "length": record["all_placements"]["length"],
                    "placements": record["all_placements"]["placements"],
                    "span_dimension": record["all_placements"]["span_dimension"],
                    "crossing_in_span": record["all_placements"]["crossing_in_span"],
                    "prefix_dimension": record["prefix_control"]["span_dimension"],
                    "prefix_contains": record["prefix_control"]["crossing_in_span"],
                }
                for record in receipt["records"]
            ],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
