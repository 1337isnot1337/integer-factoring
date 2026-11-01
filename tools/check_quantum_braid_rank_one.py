#!/usr/bin/env python3
"""All-extension rank-one loci in two six-site braid coefficient spaces."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from random import Random

import numpy as np

from check_fixed_field_descent import Field
from check_quantum_braid_path_profiles import rank
from quantum_braid_matrix_aux import coefficient_tensor as auxiliary_tensor
from quantum_braid_path import EMBED, coefficient_tensor, crossing_weights
from quantum_braid_rank_one import (
    FIELD,
    add_row,
    coefficient_basis,
    combine,
    minor_equations,
    nullspace,
    pack,
    rank_one_envelope,
    square_roots,
)

ROOT = Path(__file__).resolve().parents[1]


def arithmetic_controls(random):
    field = Field(0b111)
    intersection_controls = 0
    for _ in range(48):
        k = random.randrange(1, 8)
        width = k * (k + 1) // 2
        rows = np.array(
            [
                [random.randrange(4) for _ in range(width)]
                for _ in range(random.randrange(1, 24))
            ],
            dtype=np.uint8,
        )
        packed = {}
        for row in rows:
            add_row(packed, *pack(row))
        reference_rank = field.rank(rows.tolist())
        assert len(packed) == reference_rank
        linears = square_roots(packed, k)
        units = np.eye(width, dtype=np.uint8)[:k]
        intersection = (
            reference_rank + k - field.rank(np.concatenate((rows, units)).tolist())
        )
        assert len(linears) == intersection
        kernel = nullspace(linears, k)
        assert len(kernel) == k - intersection
        for v in kernel:
            for row in linears:
                assert np.bitwise_xor.reduce(FIELD[row, v]) == 0
        intersection_controls += 1
    return intersection_controls


def assembly_controls(basis, random):
    """Literal F16 determinants at extension-valued coordinates, not F4 points."""
    field = Field(0b10011)
    k, size, _ = basis.shape
    mixed = list(combinations(range(k), 2))
    count = 0
    for _ in range(16):
        r, s = sorted(random.sample(range(size), 2))
        c, d = sorted(random.sample(range(size), 2))
        equation = minor_equations(basis, r, s, [(c, d)])[0]
        variables = [random.randrange(16) for _ in range(k)]
        assert any(v not in EMBED for v in variables) or k == 1
        monomials = [field.mul(v, v) for v in variables] + [
            field.mul(variables[i], variables[j]) for i, j in mixed
        ]
        predicted = 0
        for a, v in zip(equation, monomials):
            predicted ^= field.mul(EMBED[int(a)], v)

        def entry(row, column):
            value = 0
            for v, a in zip(variables, basis[:, row, column]):
                value ^= field.mul(v, EMBED[int(a)])
            return value

        actual = field.mul(entry(r, c), entry(s, d)) ^ field.mul(
            entry(r, d), entry(s, c)
        )
        assert actual == predicted
        count += 1
    return count


def analyze(tensor, name, column_cut, expected_dimension, witness_index, random):
    basis, indices = coefficient_basis(tensor, column_cut)
    original = basis.copy()
    coordinates = np.eye(len(basis), dtype=np.uint8)
    controls = assembly_controls(basis, random)
    stages = []
    for stage in range(4):
        # Column cuts first have a two-dimensional Frobenius envelope. It is
        # enough to restrict to that envelope and reassemble the minors.
        target = 2 if column_cut and stage == 0 else expected_dimension
        kernel, receipt = rank_one_envelope(basis, target)
        stages.append(receipt)
        assert len(kernel) < len(basis)
        if not len(kernel):
            basis = np.zeros((0,) + basis.shape[1:], dtype=np.uint8)
            coordinates = np.zeros((0, len(original)), dtype=np.uint8)
            break
        coordinates = np.array(
            [
                np.bitwise_xor.reduce(FIELD[v[:, None], coordinates], axis=0)
                for v in kernel
            ]
        )
        basis = np.array([combine(v, basis) for v in kernel])
        if len(basis) <= expected_dimension:
            break
    assert len(basis) == expected_dimension
    field = Field(0b111)
    if expected_dimension:
        assert expected_dimension == 1
        assert field.rank(basis[0].tolist()) == 1
        assert np.array_equal(basis[0], combine(coordinates[0], original))
        witness = tensor[witness_index]
        d = witness.shape[0]
        if column_cut:
            witness = witness.transpose(0, 2, 1, 3)
        witness = witness.reshape(d * d, d * d)
        assert rank(witness) == 1
        assert rank(np.array((witness.reshape(-1), basis[0].reshape(-1)))) == 1
    return {
        "name": name,
        "cut": "column" if column_cut else "row",
        "initial_dimension": len(original),
        "basis_coefficient_indices": indices,
        "stages": stages,
        "rank_one_locus_dimension": expected_dimension,
        "rank_one_witness_coefficient": witness_index,
        "final_original_coordinates": coordinates.tolist(),
    }, controls


def check():
    random = Random(931045)
    arithmetic = arithmetic_controls(random)
    raw = coefficient_tensor(6)
    scalar = FIELD[
        np.array(crossing_weights(6), dtype=np.uint8)[:, None, None, None, None], raw
    ]
    auxiliary = auxiliary_tensor(6)
    cases, controls = [], 0
    for tensor, name, column_cut, dimension, witness in (
        (scalar, "scalar-six-row", False, 1, 36),
        (scalar, "scalar-six-column", True, 1, 0),
        (auxiliary, "auxiliary-six-row", False, 0, None),
        (auxiliary, "auxiliary-six-column", True, 1, 0),
    ):
        result, count = analyze(tensor, name, column_cut, dimension, witness, random)
        cases.append(result)
        controls += count
    # A complete four-site auxiliary spans all 4x4 matrices in either cut.
    # Hence its rank-one locus spans dimension16: a negative control against
    # wrongly treating every span enlargement as a rank-one obstruction.
    small = auxiliary_tensor(4)
    for column_cut in (False, True):
        basis, _ = coefficient_basis(small, column_cut)
        assert len(basis) == 16
        envelope, _ = rank_one_envelope(basis, 16)
        assert len(envelope) == 16
    sources = list(
        json.loads((ROOT / "evidence/quantum-braid-matrix-auxiliary.json").read_text())[
            "source_sha256"
        ]
    )
    sources += [
        "tools/quantum_braid_rank_one.py",
        "tools/check_quantum_braid_rank_one.py",
        "research/quantum-braid-rank-one-loci.md",
    ]
    return {
        "status": "pass",
        "field": "F4; consequences hold over every extension field",
        "independent_rank_and_intersection_controls": arithmetic,
        "independent_F16_minor_assembly_controls": controls,
        "full_matrix_space_negative_controls": 2,
        "cases": cases,
        "scope": "Exact six-site all-extension rank-one loci; no exclusion of higher low ranks, longer paths, or stronger factoring algorithms.",
        "source_sha256": {
            f: sha256((ROOT / f).read_bytes()).hexdigest() for f in sources
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
