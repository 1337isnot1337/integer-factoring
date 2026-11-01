#!/usr/bin/env python3
"""Reproduce the finite permutation-label profiles of the intersection lift.

NumPy accelerates an exact XOR substitution. All coefficients remain bits;
matrix ranks use independent Python-integer Gaussian elimination.
"""

import argparse
from itertools import permutations, product
import json
from pathlib import Path
import random

import numpy as np

from check_permutation_homogeneous_labels import compose, edge_masks, inverse, rank_two
from fermionic_piece_polynomial import polynomial_product
from intersection_matching_lift import lifted_polynomial


def specialize_direct(polynomial, particles):
    """Independent sparse Kronecker substitution without NumPy."""
    radix = 4 * particles + 1
    result = 0
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        code = bit.bit_length() - 1
        exponent = code % radix
        code //= radix
        for vertex in range(2 * particles):
            exponent += (code % 5) << (vertex + 1)
            code //= 5
        assert code == 0
        result ^= 1 << exponent
    return result


def run(particles):
    ps = tuple(permutations(range(particles)))
    count = len(ps)
    index = {p: i for i, p in enumerate(ps)}
    multiplication = [[index[compose(p, q)] for q in ps] for p in ps]
    inverses = [index[inverse(p)] for p in ps]
    conjugations = [
        [multiplication[multiplication[g][p]][inverses[g]] for p in range(count)]
        for g in range(count)
    ]
    representatives = []
    orbit_of = {}
    for triple in product(range(count), repeat=3):
        if triple in orbit_of:
            continue
        orbit = len(representatives)
        representatives.append(triple)
        for action in conjugations:
            orbit_of[tuple(action[p] for p in triple)] = orbit

    labels = [edge_masks(p) for p in ps]
    radix = 4 * particles + 1
    height = particles + 2 * (2 ** (2 * particles) - 1)
    maximum_degree = 4 * height
    actions = []
    weights = []
    for left, right in product(range(count), repeat=2):
        actions.append([
            multiplication[multiplication[right][p]][inverses[left]]
            for p in range(count)
        ])
        physical = ps[left] + tuple(particles + j for j in ps[right])
        weights.append([1 << (j + 1) for j in physical])
    weights = np.array(weights, dtype=np.int32)

    rows = [[0] * (count * count) for _ in range(maximum_degree + 1)]
    columns = [[0] * (count * count) for _ in range(maximum_degree + 1)]
    central_values = []
    central_code = 2 * particles + radix * 2 * sum(
        5**j for j in range(2 * particles))
    seen = bytearray(count**4)
    monomials = 0
    for x, b, u in representatives:
        polynomial = lifted_polynomial(labels[x], labels[b], labels[u],
                                       labels[0], 2 * particles)
        assert polynomial.bit_count() % 2 == int(x == b and u == 0)
        central_values.append((polynomial >> central_code) & 1)
        raw = polynomial.to_bytes((polynomial.bit_length() + 7) // 8, "little")
        codes = np.flatnonzero(np.unpackbits(
            np.frombuffer(raw, dtype=np.uint8), bitorder="little"))
        lows = codes % radix
        codes //= radix
        charges = np.empty((len(codes), 2 * particles), dtype=np.int32)
        for vertex in range(2 * particles):
            charges[:, vertex] = codes % 5
            codes //= 5
        assert not np.any(codes)
        monomials += len(charges)
        if len(charges):
            exponents = weights @ charges.T + lows
            assert int(exponents.max()) <= maximum_degree
            coefficients = np.zeros((len(actions), maximum_degree + 1), dtype=np.uint8)
            np.bitwise_xor.at(coefficients,
                             (np.arange(len(actions))[:, None], exponents), 1)
        for action_index, action in enumerate(actions):
            a1, b1, u1, y1 = action[x], action[b], action[u], action[0]
            entry = ((a1 * count + y1) * count + b1) * count + u1
            if seen[entry]:
                continue
            seen[entry] = 1
            if not len(charges):
                continue
            row, column = a1 * count + y1, b1 * count + u1
            other_row, other_column = a1 * count + u1, b1 * count + y1
            for degree in np.flatnonzero(coefficients[action_index]):
                rows[degree][row] |= 1 << column
                columns[degree][other_row] |= 1 << other_column
        polynomial_product.cache_clear()
    assert all(seen)

    summed = [0] * (count * count)
    for matrix in rows:
        for row, value in enumerate(matrix):
            summed[row] ^= value
    assert summed == [1 << row for row in range(count * count)]

    rng = random.Random(82593 + particles)
    sample_count = 24
    for _ in range(sample_count):
        a, b, u, y = (rng.randrange(count) for _ in range(4))
        direct = specialize_direct(lifted_polynomial(
            labels[a], labels[b], labels[u], labels[y], 2 * particles), particles)
        assembled = sum(
            ((matrix[a * count + y] >> (b * count + u)) & 1) << degree
            for degree, matrix in enumerate(rows)
        )
        assert assembled == direct
        polynomial_product.cache_clear()

    central_rows = [0] * (count * count)
    central_columns = [0] * (count * count)
    for a, y, b, u in product(range(count), repeat=4):
        iy = inverses[y]
        orbit = orbit_of[(multiplication[iy][a], multiplication[iy][b],
                          multiplication[iy][u])]
        if central_values[orbit]:
            central_rows[a * count + y] |= 1 << (b * count + u)
            central_columns[a * count + u] |= 1 << (b * count + y)

    profile = []
    for degree, (row_matrix, column_matrix) in enumerate(zip(rows, columns)):
        if not any(row_matrix):
            continue
        row_rank = rank_two(row_matrix)
        column_rank = rank_two(column_matrix)
        profile.append([degree - 2 * height, row_rank, column_rank])
    maximum = max(r * c for _, r, c in profile)
    return {
        "particles": particles,
        "labels": count,
        "quartet_orbits": len(representatives),
        "multigraded_monomials": monomials,
        "complete_crossing_entries_checked": count**4,
        "independent_specializations": sample_count,
        "uniform_central_ranks": [rank_two(central_rows), rank_two(central_columns)],
        "nonzero_univariate_pieces": len(profile),
        "maximum_rank_product": maximum,
        "worst": [entry for entry in profile if entry[1] * entry[2] == maximum],
        "profile": profile,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum-particles", type=int, default=4, choices=(2, 3, 4))
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    result = {
        "coefficient_field": "F2",
        "kernel": "intersection",
        "scope": "complete finite profiles; no asymptotic extrapolation",
        "cases": [run(particles) for particles in range(2, options.maximum_particles + 1)],
    }
    expected = {2: (39, 8, 2), 3: (221, 792, 16), 4: (957, 207360, 364)}
    for case in result["cases"]:
        pieces, product_bound, rank = expected[case["particles"]]
        assert case["nonzero_univariate_pieces"] == pieces
        assert case["maximum_rank_product"] == product_bound
        assert case["uniform_central_ranks"] == [rank, rank]
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if options.output:
        options.output.write_text(encoded)
    print(json.dumps({
        "scope": result["scope"],
        "cases": [{key: value for key, value in case.items()
                   if key not in ("profile", "worst")} for case in result["cases"]],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
