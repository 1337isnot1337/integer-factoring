#!/usr/bin/env python3
"""Exact natural-kernel duality and downset Gaussian-factor controls."""
from __future__ import annotations
import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
from check_schur_gaussian_reduction import polynomial_multiply as mul
from reflection_dual_rail_gaussian import gauged_kernel, logical_heights, physical_bit_energies
from restrict_implicit_reflection import restricted_tensor


def pfaffian_function(kernel):
    @lru_cache(None)
    def pf(mask):
        if not mask:
            return 1
        if mask.bit_count() % 2:
            return 0
        first = mask & -mask
        i, rest = first.bit_length() - 1, mask ^ first
        remaining, value = rest, 0
        while remaining:
            chosen = remaining & -remaining
            remaining ^= chosen
            j = chosen.bit_length() - 1
            value ^= mul(kernel[i][j], pf(rest ^ chosen))
        return value
    return pf


def determinant(matrix):
    result = 0
    for permutation in permutations(range(len(matrix))):
        term = 1
        for row, column in enumerate(permutation):
            term = mul(term, matrix[row][column])
        result ^= term
    return result


def multiply(first, second):
    result = [[0] * len(second[0]) for _ in first]
    for i, row in enumerate(first):
        for k, value in enumerate(row):
            if value:
                for j, other in enumerate(second[k]):
                    if other:
                        result[i][j] ^= mul(value, other)
    return result


def restrict(matrix, rows, columns):
    return [[matrix[i][j] for j in columns] for i in rows]


def reverse(value, degree):
    result = 0
    while value:
        bit = value & -value
        value ^= bit
        exponent = bit.bit_length() - 1
        assert exponent <= degree
        result ^= 1 << (degree - exponent)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    kernel_entries = coefficient_entries = 0
    for interval in (4, 8, 16, 32, 64, 128):
        modes = interval // 2
        zero, infinity = gauged_kernel(interval, 'zero'), gauged_kernel(interval, 'infinity')
        permutation = [((0, 2, 1, 3)[i // modes]) * modes + i % modes for i in range(4 * modes)]
        for i in range(4 * modes):
            for j in range(4 * modes):
                assert infinity[i][j] == zero[permutation[i]][permutation[j]]
                kernel_entries += 1
    weights = physical_bit_energies(16)
    energy = lambda symbol: sum(value for bit, value in enumerate(weights) if symbol >> bit & 1)
    hot = sorted(range(256), key=lambda symbol: (energy(symbol), symbol))[:8]
    for interval, first, second in ((4, list(range(4)), list(range(4))),
                                    (8, [0, 1, 3], [0, 2]), (16, hot, hot)):
        seed = restricted_tensor(interval, first, second)
        modes, heights = interval // 2, logical_heights(interval, 'zero')
        state_height = lambda symbol: sum(heights[bit] for bit in range(modes) if symbol >> bit & 1)
        qx, qy = len(first), len(second)
        for index, (value, imaginary) in enumerate(seed['values']):
            assert imaginary == 0
            x, y = index % qx, index // qx % qy
            u, v = index // (qx * qy) % qy, index // (qx * qy * qy)
            swap = x + qx * (u + qy * (y + qy * v))
            bound = interval * interval + state_height(first[v]) - state_height(first[x])
            assert reverse(value, bound) == seed['values'][swap][0]
            coefficient_entries += 1
    # A separate polynomial Pfaffian expansion against permutation
    # determinants proves the complete small Gaussian factorization.
    kernel = gauged_kernel(4, 'zero')
    modes, half, masks = 2, 4, 16
    left = [row[:half] for row in kernel[:half]]
    right = [row[half:] for row in kernel[half:]]
    cross = [row[half:] for row in kernel[:half]]
    full_pf, left_pf, right_pf = map(pfaffian_function, (kernel, left, right))
    full = [[full_pf(row | (column << half)) for column in range(masks)] for row in range(masks)]
    U = [[left_pf(row ^ smaller) if smaller & ~row == 0 else 0
          for smaller in range(masks)] for row in range(masks)]
    V = [[right_pf(column ^ smaller) if smaller & ~column == 0 else 0
          for column in range(masks)] for smaller in range(masks)]
    E = [[0] * masks for _ in range(masks)]
    for row in range(masks):
        rows = [i for i in range(half) if row >> i & 1]
        for column in range(masks):
            columns = [i for i in range(half) if column >> i & 1]
            if len(rows) == len(columns):
                E[row][column] = determinant(restrict(cross, rows, columns))
    assert multiply(multiply(U, E), V) == full
    factor_entries = masks * masks
    endpoint_entries = charge_zero_entries = 0
    cases = []
    for first, second in (([0, 1, 2], [0, 1, 2]), ([0, 1, 2], [0, 1]),
                          ([0, 2], [0, 1, 2]), (list(range(4)), [0, 2])):
        assert all(smaller in code for code in (first, second) for state in code
                   for smaller in range(1 << modes) if smaller & ~state == 0)
        row_pairs = [(x, y) for x in first for y in second]
        column_pairs = [(u, v) for u in second for v in first]
        rows = [x | (y << modes) for x, y in row_pairs]
        columns = [u | (v << modes) for u, v in column_pairs]
        restricted = restrict(E, rows, columns)
        reconstructed = multiply(multiply(restrict(U, rows, rows), restricted), restrict(V, columns, columns))
        assert reconstructed == restrict(full, rows, columns)
        factor_entries += len(rows) * len(columns)
        for i, (x, y) in enumerate(row_pairs):
            for j, (u, v) in enumerate(column_pairs):
                assert restricted[i][j].bit_count() % 2 == int(x == v and y == u)
                endpoint_entries += 1
                if rows[i].bit_count() != columns[j].bit_count():
                    assert restricted[i][j] == 0
                    charge_zero_entries += 1
        cases.append({'first_code': first, 'second_code': second, 'cut_dimension': len(rows)})
    # The singleton non-downset misses a contribution through its vacuum.
    absent_vacuum = 1 | (1 << modes)
    assert full[absent_vacuum][absent_vacuum] != E[absent_vacuum][absent_vacuum]
    assert U[absent_vacuum][absent_vacuum] == V[absent_vacuum][absent_vacuum] == 1
    paths = [Path('evidence/implicit-reflection-code-controls.json'),
             Path('evidence/stable-hot-code-rational-pressure.json')]
    profile_pairs = []
    for path in paths:
        record = json.loads(path.read_bytes())
        orders = (record['large_code']['natural_profiles_infinity_then_zero']
                  if 'large_code' in record else record['integer_profiles_infinity_then_zero'])
        assert min(orders[0]) + min(orders[1]) == 0
        shifted = [[value - min(side) for value in side] for side in orders]
        assert shifted[0] == shifted[1]
        profile_pairs.append({'path': str(path), 'sha256': sha256(path.read_bytes()).hexdigest(),
                              'channels_per_cut': len(orders[0]), 'opposite_shifts': [min(side) for side in orders]})
    names = (Path(__file__).name, 'restrict_implicit_reflection.py', 'reflection_dual_rail_gaussian.py',
             'check_schur_gaussian_reduction.py')
    report = {'scope': 'Exact uniform natural-cut duality and downset compound reduction; no new factoring exponent.',
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'counts': {'kernel_entry_identities': kernel_entries, 'complete_coefficient_dualities': coefficient_entries,
                         'polynomial_factorization_entries': factor_entries, 'compound_endpoint_entries': endpoint_entries,
                         'off_charge_zero_entries': charge_zero_entries, 'non_downset_countercontrols': 1},
              'downset_cases': cases, 'full_profile_symmetry_controls': profile_pairs,
              'non_downset_difference': full[absent_vacuum][absent_vacuum] ^ E[absent_vacuum][absent_vacuum],
              'passed': True}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'counts': report['counts']}))


if __name__ == '__main__':
    main()
