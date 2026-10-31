#!/usr/bin/env python3
"""Four-copy endpoint circuits, literal gauge clearing, and full profile invariance."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from pathlib import Path
import random
from types import SimpleNamespace

from binary_typed_endpoint_gauges import (DIRECTIONS, CIRCUITS, circuit_certificate,
    eighth_power_potential, chi_power, labels, strict_laurent_gauge)
from full_alphabet_two_cut_lattices import cut_location
from weyl_lattice_filter_synthesis import polynomial_multiply
from native_matrix_smith_f4 import MatrixOracle
from check_full_alphabet_two_cut_lattices import independent_profiles


def balanced_multisets():
    result = []
    for length in range(1, 5):
        for indices in combinations_with_replacement(range(8), length):
            if all(sum(DIRECTIONS[i][coordinate] for i in indices) == 0 for coordinate in range(2)):
                result.append(indices)
    return tuple(result)


def circuit_controls():
    balanced = balanced_multisets()
    assert Counter(sum(weights) for indices, weights in CIRCUITS) == {2: 4, 3: 4, 4: 4}
    for indices, weights in CIRCUITS:
        assert all(sum(weight*DIRECTIONS[i][j] for i, weight in zip(indices, weights)) == 0 for j in range(2))
    complete = feasible = random_cases = 0
    certificates = []
    for values in product((-1, 0, 1), repeat=8):
        expected = all(sum(values[i] for i in indices) > 0 for indices in balanced)
        actual = circuit_certificate(dict(zip(DIRECTIONS, values)))['valid']
        assert actual == expected; complete += 1
        if actual:
            result = eighth_power_potential(dict(zip(DIRECTIONS, values)))
            assert result['valid']; feasible += 1
            certificates.append([values, result['integer_potential']])
    rng = random.Random(318181)
    for _ in range(600):
        support = sorted(rng.sample(range(8), rng.randrange(1, 9)))
        values = {DIRECTIONS[i]: rng.randrange(-5, 8) for i in support}
        expected = all(sum(values[DIRECTIONS[i]] for i in indices) > 0 for indices in balanced
                       if all(i in support for i in indices))
        result = eighth_power_potential(values)
        assert result['valid'] == expected; random_cases += 1
        if result['valid']:
            f, g = result['integer_potential']
            assert all(8*value+f*a+g*b >= 2 for (a, b), value in values.items())
    sharp = {(1, 1): 1, (-1, 1): 1, (0, -1): -1}
    available = [indices for indices in balanced if all(DIRECTIONS[i] in sharp for i in indices)]
    assert all(len(indices) == 4 for indices in available)
    assert not eighth_power_potential(sharp)['valid']
    # Literal rational entries in the sharp example use common denominator chi.
    sharp_values = [(chi_power(1) if (lambda t: t[0] == t[3] and t[1] == t[2])(labels(i)) else 0, 0)
                    for i in range(16)]
    sharp_values[3] = sharp_values[10] = (chi_power(2), 0)
    sharp_values[4] = (1, 0)
    assert not strict_laurent_gauge(sharp_values, 1)['valid']
    literal_product = 1
    for entry in (3, 10, 4, 4):
        literal_product = polynomial_multiply(literal_product, sharp_values[entry][0])
    assert literal_product == chi_power(4)
    return {'primitive_circuits': len(CIRCUITS),
        'primitive_circuit_lengths': dict(Counter(sum(weights) for indices, weights in CIRCUITS)),
        'all_balanced_multisets_through_four': len(balanced),
        'complete_valuation_assignments': complete, 'feasible_complete_assignments': feasible,
        'independent_partial_support_assignments': random_cases,
        'sharp_first_failure_length': 4,
        'complete_feasible_potentials_sha256': sha256(json.dumps(certificates, separators=(',', ':')).encode()).hexdigest()}


def profiles(values, pole, heights, independent=False):
    result = []
    for side in range(2):
        matrix = [[(0, 0)]*4 for _ in range(4)]
        for entry, pair in enumerate(values):
            i, j = cut_location(labels(entry), 2, side)
            matrix[i][j] = pair
        h = heights[4*side:4*side+4]
        rowlegs = (0, 2) if side == 0 else (0, 1)
        collegs = (1, 3) if side == 0 else (2, 3)
        row_weights = [h[rowlegs[0]][a]+h[rowlegs[1]][b] for a, b in product(range(2), repeat=2)]
        col_weights = [h[collegs[0]][a]+h[collegs[1]][b] for a, b in product(range(2), repeat=2)]
        raw = MatrixOracle(matrix, 4).profile(row_weights, col_weights, infinity=side == 0)['orders']
        assert all(abs(value-round(value)) < 1e-8 for value in raw)
        result.append(tuple(round(value)+(pole if side == 0 else 0) for value in raw))
    counts = 0
    if independent:
        degree = max(max(pair).bit_length()-1 for pair in values)
        scalars = [-sum(min(row) for row in heights[4*side:4*side+4])+(degree if side == 0 else 0)
                   for side in range(2)]
        model = SimpleNamespace(q=2, size=4, heights=heights, degree=degree,
            scalar_shifts=scalars, labels=[labels(i) for i in range(16)])
        alternate, counts = independent_profiles(model, values)
        alternate = tuple(tuple(value+(pole if side == 0 else 0) for value in row)
                          for side, row in enumerate(alternate))
        assert tuple(result) == alternate
    return tuple(result), counts


def tensor_controls():
    rng = random.Random(318182)
    rows = []; ranks = 0
    for trial in range(16):
        strict = []
        for entry in range(16):
            x, y, u, v = labels(entry)
            low = rng.randrange(16); high = rng.randrange(16)
            strict.append((int(x == v and y == u)^low^(low << 1), high^(high << 1)))
        a, b = rng.randrange(-3, 4), rng.randrange(-3, 4)
        pole = abs(a)+abs(b)
        rational = []
        for entry, pair in enumerate(strict):
            x, y, u, v = labels(entry)
            factor = chi_power(pole-a*(x-v)-b*(y-u))
            rational.append(tuple(polynomial_multiply(value, factor) for value in pair))
        gauged = strict_laurent_gauge(rational, pole)
        assert gauged['valid']
        f, g = gauged['integer_potential']
        heights = [(0, rng.randrange(-3, 4)) for _ in range(8)]
        new_heights = [tuple(8*value for value in row) for row in heights]
        for leg, shift in enumerate((f, g, -g, -f)):
            new_heights[leg] = (0, new_heights[leg][1]+shift)
        old, count = profiles(rational, pole, heights, independent=trial < 2); ranks += count
        new, count = profiles(gauged['values'], 0, new_heights, independent=trial < 2); ranks += count
        assert new == tuple(tuple(8*value for value in row) for row in old)
        for legs in ((0, 3), (1, 2)):
            old_energy = sum(heights[side+leg][1] for side in (0, 4) for leg in legs)
            new_energy = sum(new_heights[side+leg][1] for side in (0, 4) for leg in legs)
            assert new_energy == 8*old_energy
        rows.append({'original_integer_potential': [a, b], 'original_common_pole': pole,
            'clearing_integer_potential': [f, g], 'minimum_transformed_offmatch_order': gauged['minimum_transformed_margin'],
            'original_profiles': old, 'strict_profiles': new, 'strict_degree': gauged['degree']})
    return {'full_rational_to_strict_tensor_cases': len(rows),
        'complete_native_profiles': 4*len(rows), 'independent_toeplitz_ranks': ranks,
        'literal_strict_endpoint_entries': 16*len(rows), 'cases': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = {'status': 'Exact four-copy endpoint criterion and pressure-preserving strict gauges passed.',
        'circuits': circuit_controls(), 'tensors': tensor_controls(),
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in (Path(__file__).name, 'binary_typed_endpoint_gauges.py',
                         'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
                         'check_full_alphabet_two_cut_lattices.py')}}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'circuits': result['circuits'],
        'tensors': {key: value for key, value in result['tensors'].items() if key != 'cases'}}))


if __name__ == '__main__':
    main()
