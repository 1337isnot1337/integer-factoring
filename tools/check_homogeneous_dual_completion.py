#!/usr/bin/env python3
"""Reconstruct and independently check the 88-dimensional crossing family."""

from pathlib import Path
import json
import random

from check_modular_crossing import determinant_operator_two, homogeneous_piece_two
from homogeneous_dual_completion import build_completion_seed, completion_basis, specialize


def binary_basis(rows):
    pivots, independent, relations = {}, [], []
    for index, row in enumerate(rows):
        combination = 1 << index
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (row, combination)
                independent.append(index)
                break
            row ^= pivots[pivot][0]
            combination ^= pivots[pivot][1]
        if not row:
            relations.append([j for j in range(len(rows)) if combination >> j & 1])
    return independent, relations


def verify():
    root = Path(__file__).resolve().parents[1]
    saved = json.loads((root / "evidence/homogeneous-dual-completion-q16-h32.json").read_text())
    fresh = build_completion_seed()
    # JSON normalization converts tuples to lists without changing coefficients.
    assert json.loads(json.dumps(fresh)) == saved
    seed = saved
    labels, reflected, kernels = completion_basis()
    dimension = len(labels)
    assert len(kernels) == 15 and dimension == 8
    assert seed["dual_parameter_dimension"] == 120
    assert seed["common_laurent_shift"] == -64
    assert seed["verified_exponent_step"] == 2
    kernel_entries = 0
    for terms in kernels:
        for label in labels:
            assert sum(determinant_operator_two(term, (), label) for term in terms) % 2 == 0
            kernel_entries += 1
    bit_width = max(value.bit_length() for slab in seed["kernel_slabs"] for value in slab)
    vectors = [sum(value << (bit_width * j) for j, value in enumerate(slab))
               for slab in seed["kernel_slabs"]]
    independent, relations = binary_basis(vectors)
    assert independent == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11]
    assert relations == [[6, 8, 9, 10], [6, 8, 11, 12],
                         [6, 9, 11, 13], [8, 9, 11, 14]]
    for relation in relations:
        total = 0
        for index in relation:
            total ^= vectors[index]
        assert not total

    # Independent original recurrence: complete y,u,v slices at selected modes.
    source_entries = 0
    for kernel_index in (2, 6, 8, 11):
        terms = kernels[kernel_index]
        duals = tuple(terms for _ in labels)
        for shift in (-64, -32, -16, 0):
            rows = homogeneous_piece_two(32, shift, labels, duals)
            packed_index = (shift + 64) // 2
            slab = seed["kernel_slabs"][kernel_index]
            for yi in range(dimension):
                for vi in range(dimension):
                    for ui in range(dimension):
                        expected = (rows[yi] >> (vi * dimension + ui)) & 1
                        actual = (slab[yi + dimension * ui + dimension ** 2 * vi]
                                  >> packed_index) & 1
                        assert actual == expected
                        source_entries += 1

    endpoint_entries = 0
    rng = random.Random(619830)
    for count in (0, 1, 3, 11, 40):
        changes = [(rng.randrange(dimension), rng.randrange(15), rng.randrange(1, 4))
                   for _ in range(count)]
        values = specialize(seed, changes)
        for xi in range(dimension):
            for yi in range(dimension):
                for ui in range(dimension):
                    for vi in range(dimension):
                        lo, hi = values[xi + dimension * yi + dimension ** 2 * ui
                                        + dimension ** 3 * vi]
                        assert (lo.bit_count() % 2, hi.bit_count() % 2) == (
                            int(xi == vi and yi == ui), 0)
                        endpoint_entries += 1
    return {"status": "exact controls passed; no pressure improvement asserted",
            "full_parameter_dimension": 120,
            "kernel_slab_rank": len(independent),
            "complete_tensor_direction_dimension": dimension * len(independent),
            "independent_kernel_indices": independent,
            "coefficientwise_kernel_relations": relations,
            "vacuum_kernel_entries": kernel_entries,
            "independent_source_piece_entries": source_entries,
            "F4_endpoint_entries": endpoint_entries,
            "full_family_reconstructed": True,
            "maximum_polynomial_degree": bit_width - 1}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
