#!/usr/bin/env python3
"""Independent exact endpoint and F2 restriction-of-scalars profile controls."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_block_smith_crossings import reverse_polynomial, smith_profile, toeplitz_rank
from physical_split_crossings import (
    F4_PRODUCT, SplitProfile, check_split, flatten, hyperplane_split,
    pair_split, retract,
)
from reflection_dual_rail_gaussian import tensor


def literal_retraction(seed, left_x, right_x, left_y, right_y):
    ambient_x, ambient_y = seed['q_x'], seed['q_u']
    maps = (left_x, left_y, right_y, right_x)
    dimensions = [len(matrix) for matrix in maps]
    result = []
    for reverse_output in product(*(range(d) for d in dimensions[::-1])):
        output = reverse_output[::-1]
        supports = [[(i, c) for i, c in enumerate(matrix[symbol]) if c]
                    for matrix, symbol in zip(maps, output)]
        low = high = 0
        for term in product(*supports):
            (x, a), (y, b), (u, c), (v, d) = term
            scalar = F4_PRODUCT[F4_PRODUCT[a][b]][F4_PRODUCT[c][d]]
            entry = seed['values'][x + ambient_x * y + ambient_x * ambient_y * u
                                   + ambient_x * ambient_y * ambient_y * v]
            # Multiply individual coefficients, independently of the bit-plane shortcut.
            for exponent in range(max(part.bit_length() for part in entry)):
                coefficient = (entry[0] >> exponent & 1) | ((entry[1] >> exponent & 1) << 1)
                value = F4_PRODUCT[scalar][coefficient]
                low ^= (value & 1) << exponent
                high ^= (value >> 1) << exponent
        result.append((low, high))
    return result


def scalar_restriction(matrix):
    size = len(matrix)
    result = [[0] * (2 * size) for _ in range(2 * size)]
    for i, row in enumerate(matrix):
        for j, (low, high) in enumerate(row):
            result[2 * i][2 * j] = low
            result[2 * i][2 * j + 1] = high
            result[2 * i + 1][2 * j] = high
            result[2 * i + 1][2 * j + 1] = low ^ high
    return result


def check_case(base, left_x, right_x, energy_x, left_y=None, right_y=None, energy_y=None):
    if left_y is None:
        left_y, right_y, energy_y = left_x, right_x, energy_x
    projected = retract(base, left_x, right_x, left_y, right_y)
    assert projected['values'] == literal_retraction(base, left_x, right_x, left_y, right_y)
    oracle = SplitProfile(projected)
    native = oracle.profiles(energy_x, energy_y)
    profiles = []
    for cut in (0, 1):
        matrix, rows, columns = flatten(projected, cut)
        degree = projected['degree']
        shifted = []
        for (a, b), row in zip(rows, matrix):
            values = []
            for (c, d), entry in zip(columns, row):
                if cut == 0:
                    shift = (energy_x[a] + energy_y[b]) // 2
                    values.append(tuple(part << shift for part in entry))
                else:
                    shift = (energy_y[b] + energy_x[d]) // 2
                    values.append(tuple(reverse_polynomial(part, degree) << shift for part in entry))
            shifted.append(values)
        binary = scalar_restriction(shifted)
        bound = max(value.bit_length() - 1 for row in binary for value in row)
        doubled_orders, precision = smith_profile(binary, bound)
        expected = tuple(sorted(round(value + (degree if cut else 0))
                                for value in native[cut] for _ in range(2)))
        assert doubled_orders == expected
        assert all(abs(value - round(value)) < 1e-9 for value in native[cut])
        checks = 0
        for length in sorted({1} | {value + 1 for value in doubled_orders}):
            assert toeplitz_rank(binary, length) == sum(max(0, length - value)
                                                      for value in doubled_orders)
            checks += 1
        profiles.append({'cut': cut, 'orders': dict(Counter(native[cut])),
                         'independent_precision': precision, 'toeplitz_checks': checks})
    return {'dimensions': projected['dimensions'], 'profiles': profiles,
            'endpoint_entries': len(projected['values']),
            'nonbinary_entries': sum(bool(value[1]) for value in projected['values'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = []
    base = tensor(4)
    examples = [
        ((0, 0, 0, 1), (0, 0, 0, 1)),
        ((1, 0, 0, 0), (1, 2, 3, 1)),
        ((1, 1, 0, 0), (2, 3, 1, 2)),
        ((1, 2, 3, 1), (1, 0, 0, 0)),
        ((1, 2, 3, 1), (0, 3, 0, 0)),
        ((0, 1, 3, 2), (2, 1, 0, 0)),
    ]
    natural = [0, 12, 4, 16]
    for functional, complement in examples:
        alphabet, left, right = hyperplane_split(functional, complement)
        cases.append(check_case(base, left, right, [natural[i] for i in alphabet]))
        cases.append(check_case(base, left, right, [0, 6, 18]))

    larger = tensor(8)
    for pair, lp, rp in [((6, 9), (2, 3), (1, 1)),
                         ((7, 15), (1, 0), (1, 1)),
                         ((0, 15), (3, 1), (1, 2))]:
        alphabet, left, right = pair_split(16, pair, lp, rp)
        keep = sorted({alphabet.index(pair[0])} | set(range(3)))[:3]
        # Include the mixed symbol even if it lies outside the first coordinates.
        if alphabet.index(pair[0]) not in keep:
            keep[-1] = alphabet.index(pair[0])
        left, right = [left[i] for i in keep], [right[i] for i in keep]
        cases.append(check_case(larger, left, right, [0, 8, 24]))
        cases.append(check_case(larger, left, right, [0, 8, 24],
                                left[:2], right[:2], [0, 16]))

    rejected = False
    try:
        check_split([[1, 1]], [[1, 1]])
    except ValueError:
        rejected = True
    assert rejected
    sources = [Path(__file__), Path(__file__).with_name('physical_split_crossings.py'),
               Path(__file__).with_name('reflection_dual_rail_gaussian.py'),
               Path(__file__).with_name('check_block_smith_crossings.py')]
    report = {
        'status': 'pass', 'cases': cases, 'invalid_split_rejected': rejected,
        'endpoint_entries': sum(case['endpoint_entries'] for case in cases),
        'independent_profiles': 2 * len(cases),
        'toeplitz_checks': sum(p['toeplitz_checks'] for case in cases for p in case['profiles']),
        'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in sources},
        'scope': 'exact split maps, literal four-leg contraction, complete endpoints, F2 scalar restriction Smith profiles and Toeplitz ranks',
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
