#!/usr/bin/env python3
"""Exact finite checks of the far-reflection support reduction."""

from itertools import product
import argparse
import json

from check_coordinate_support_rank import local_polynomial, supported
from check_modular_crossing import column_bits, rank_two, regular_partitions
from fermionic_piece_polynomial import (
    monomial_piece_polynomial, piece_matrices, reflected_duals,
)


def incidence(label):
    return sum(1 << ((part - 1) // 2) for part in label)


def substitute_power(polynomial, power):
    result = 0
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        result ^= 1 << (power * (bit.bit_length() - 1))
    return result


def outside_support(rows, labels):
    dimension = len(labels)
    codes = [incidence(label) for label in labels]
    result = []
    for r, row in enumerate(rows):
        x, y = divmod(r, dimension)
        for column in range(dimension ** 2):
            if row >> column & 1:
                v, u = divmod(column, dimension)
                if not supported(codes[x], codes[y], codes[u], codes[v], 0):
                    result.append([x, y, u, v])
    return result


def verify():
    fixed_cases = [(8, 2), (12, 2), (16, 2), (20, 2), (15, 3), (21, 3), (24, 4)]
    polynomial_entries = 0
    fixed_records = []
    for weight, particles in fixed_cases:
        labels = tuple(label for label in regular_partitions(2, weight)
                       if len(label) == particles)[:6]
        largest = max(max(label) for label in labels)
        modes = 1 << largest.bit_length()
        codes = [incidence(label) for label in labels]
        duals = reflected_duals(labels, modes)
        for x, y, v, u in product(range(len(labels)), repeat=4):
            polynomial, offset = monomial_piece_polynomial(
                duals[x][0][0], labels[v], labels[u], labels[y])
            reference = local_polynomial(codes[x], codes[y], codes[u], codes[v],
                                         (largest + 1) // 2)
            assert offset == -2 * modes * particles
            assert polynomial == substitute_power(reference, modes)
            polynomial_entries += 1
        fixed_records.append({"weight": weight, "particles": particles,
                              "labels": labels, "mode_box": modes})

    projected_entries = mask_rank_checks = 0
    records = []
    saw_mask_necessity = False
    for weight in (9, 15, 16, 21):
        labels = regular_partitions(2, weight)
        dimension = len(labels)
        largest = max(max(label) for label in labels)
        modes = 1 << largest.bit_length()
        max_particles = max(map(len, labels))
        pieces = piece_matrices(labels, reflected_duals(labels, modes))
        mask = []
        for x, y in product(range(dimension), repeat=2):
            row = 0
            for v, u in product(range(dimension), repeat=2):
                if len(labels[x]) + len(labels[u]) == len(labels[v]) + len(labels[y]):
                    row |= 1 << (v * dimension + u)
            mask.append(row)
        mask_row_rank = rank_two(mask)
        mask_column_rank = rank_two(column_bits(mask, dimension))
        assert max(mask_row_rank, mask_column_rank) <= 2 * max_particles + 1
        total = [0] * dimension ** 2
        projected_shifts = []
        maximum = 0
        unsupported_before = 0
        for shift, rows in sorted(pieces.items()):
            unsupported_before += len(outside_support(rows, labels))
            projected = [a & b for a, b in zip(rows, mask)]
            assert not outside_support(projected, labels)
            if any(projected):
                assert shift % (2 * modes) == 0
                assert -2 * modes * max_particles <= shift <= 0
                projected_shifts.append(shift)
            row_rank = rank_two(rows)
            column_rank = rank_two(column_bits(rows, dimension))
            projected_row_rank = rank_two(projected)
            projected_column_rank = rank_two(column_bits(projected, dimension))
            assert projected_row_rank <= mask_row_rank * row_rank
            assert projected_column_rank <= mask_column_rank * column_rank
            diagonal = sum((row >> i) & 1 for i, row in enumerate(projected))
            assert dimension * projected_row_rank * projected_column_rank >= diagonal ** 2
            total = [a ^ b for a, b in zip(total, projected)]
            maximum = max(maximum, row_rank * column_rank)
            projected_entries += dimension ** 4
            mask_rank_checks += 2
        assert total == [1 << i for i in range(dimension ** 2)]
        assert len(projected_shifts) <= max_particles + 1
        assert maximum * (max_particles + 1) ** 2 * (2 * max_particles + 1) ** 2 >= dimension ** 3
        saw_mask_necessity |= bool(unsupported_before)
        records.append({"weight": weight, "dimension": dimension, "mode_box": modes,
                        "max_particles": max_particles, "projected_shifts": projected_shifts,
                        "mask_ranks": [mask_row_rank, mask_column_rank],
                        "unsupported_entries_before_projection": unsupported_before,
                        "maximum_original_rank_product": maximum})
    assert saw_mask_necessity

    # A box outside the theorem's range has a genuine third-pairing entry.
    labels = tuple(label for label in regular_partitions(2, 16) if len(label) == 2)
    small_box = piece_matrices(labels, reflected_duals(labels, 8))
    counterexample = outside_support(small_box[0], labels)
    assert counterexample and max(max(label) for label in labels) > 8
    return {
        "status": "pass",
        "scope": "exact finite checks of the proved far-reflection obstruction",
        "fixed_particle_polynomial_entries": polynomial_entries,
        "charge_projected_coefficient_entries": projected_entries,
        "hadamard_rank_inequalities": mask_rank_checks,
        "fixed_particle_cases": fixed_records,
        "varying_particle_cases": records,
        "small_box_negative_control": {"weight": 16, "mode_box": 8,
                                        "labels": labels, "shift": 0,
                                        "unsupported_index": counterexample[0]},
        "new_factoring_exponent_claimed": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    encoded = json.dumps(verify(), indent=2, sort_keys=True) + "\n"
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(encoded)
    print(encoded, end="")
