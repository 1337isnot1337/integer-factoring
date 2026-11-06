#!/usr/bin/env python3
"""Identify the one coefficient direction that completes the ten-site algebra."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from bounded_spin_batched_blocks import coefficient_tensor
from check_bounded_spin_growth import original_coefficients
from check_bounded_spin_six_vertex import operators
from check_bounded_spin_transposition_subalgebra import (
    frame, replica_swap, transposition_class_sum, transposition_matrices,
)
from native_gf5_flint import inverse, multiply, pivot_columns, rank

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / 'evidence/bounded-spin-ten-algebra-model.json'
OUTPUT = ROOT / 'evidence/bounded-spin-missing-direction.json'
LENGTH = 10


def main() -> None:
    model = json.loads(MODEL.read_text())
    monomials, tensors = coefficient_tensor(LENGTH)
    old = original_coefficients(monomials, tensors, LENGTH)
    operators_all = operators(old, LENGTH)[0]
    basis = operators_all[model['basis_indices']].reshape(8, -1)
    paths, words, vectors, inverse_norms = frame(LENGTH)
    transpositions = transposition_matrices(LENGTH, words, vectors, inverse_norms)
    class_sum = transposition_class_sum(transpositions)
    replica = replica_swap(len(paths))

    columns = np.array(pivot_columns(basis))
    assert len(columns) == 8
    inverse_frame = inverse(basis[:, columns])
    class_sum_coordinates = multiply(class_sum.reshape(-1)[columns][None, :],
                                     inverse_frame)[0]
    replica_coordinates = multiply(replica.reshape(-1)[columns][None, :],
                                   inverse_frame)[0]
    assert np.array_equal(multiply(class_sum_coordinates[None, :], basis)[0],
                          class_sum.reshape(-1))
    assert np.array_equal(multiply(replica_coordinates[None, :], basis)[0],
                          replica.reshape(-1))
    assert replica_coordinates.tolist() == [1, 0, 0, 0, 0, 0, 0, 0]
    assert class_sum_coordinates.tolist() == [0, 2, 0, 1, 0, 2, 2, 0]

    regular = np.array(model['regular'], dtype=np.uint8)
    unit = np.array(model['unit'], dtype=np.uint8)

    def left(v):
        return (np.tensordot(v, regular, axes=(0, 0)) % 5).astype(np.uint8)

    def product(v, w):
        return (left(v) @ w) % 5

    def generated(generators):
        space = [unit]
        frontier = [unit]
        while frontier:
            new = []
            for v in frontier:
                for g in generators:
                    w = product(v, g)
                    if rank(np.array(space + [w])) > len(space):
                        space.append(w)
                        new.append(w)
            frontier = new
        return np.array(space, dtype=np.uint8)

    subalgebra = generated([class_sum_coordinates, replica_coordinates])
    assert len(subalgebra) == 6
    missing = np.zeros(8, dtype=np.uint8)
    missing[4] = 1
    with_missing = generated([class_sum_coordinates, replica_coordinates, missing])
    assert len(with_missing) == 8
    assert model['basis_indices'][4] == 5
    # The other quotient direction alone does not complete the algebra.
    other = np.zeros(8, dtype=np.uint8)
    other[2] = 1
    assert len(generated([class_sum_coordinates, replica_coordinates, other])) == 7
    # The sum of the two quotient directions separates the three primaries in 335.
    separator = (other + missing) % 5
    characters = np.array(model['characters'], dtype=np.int64)
    values = [int(sum(int(separator[j]) * characters[i, j] for j in range(8)) % 5)
              for i in range(7)]
    dimensions = model['dimensions']
    target = [i for i, dimension in enumerate(dimensions) if dimension in (34, 1, 300)]
    target_values = [values[i] for i in target]
    assert target_values == [0, 1, 2]
    assert len(set(target_values)) == 3

    receipt = {
        'passed': True,
        'length': LENGTH,
        'prime': model['prime'],
        'basis_indices': model['basis_indices'],
        'replica_swap_coordinates': replica_coordinates.tolist(),
        'class_sum_coordinates': class_sum_coordinates.tolist(),
        'generated_dimension': len(subalgebra),
        'completing_coordinate_index': 4,
        'completing_polynomial_coefficient': model['basis_indices'][4],
        'dimension_with_completing_direction': len(with_missing),
        'dimension_with_other_quotient_direction': 7,
        'separator_coordinates': separator.tolist(),
        'separator_character_values': values,
        'coarse_335_primary_indices': target,
        'coarse_335_separator_values': target_values,
        'model_sha256': sha256(MODEL.read_bytes()).hexdigest(),
        'constructor_sha256': sha256(
            (ROOT / 'tools/bounded_spin_six_vertex.py').read_bytes()).hexdigest(),
        'scope': ('A ten-site structural diagnostic: the replica swap and transposition '
                  'class sum plus one coefficient direction generate the full eight-'
                  'dimensional algebra. No twelve-site model, rank bound, or score '
                  'change is claimed.'),
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'generated': len(subalgebra),
                      'completed_by_coefficient': model['basis_indices'][4],
                      'separator_values_34_1_300': target_values}))


if __name__ == '__main__':
    main()
