#!/usr/bin/env python3
"""Independent exact controls for polynomial physical boundary retractions."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from check_coherent_hot_boundary import logical, natural_profiles
from coherent_hot_boundary import compress as constant_compress
from moving_boundary_crossings import compress, evaluate_at_one, multiply, paired_maps
from physical_split_crossings import F4_PRODUCT
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor

ROOT = Path(__file__).resolve().parents[1]


def coefficients(value):
    return [(value[0] >> i & 1) | ((value[1] >> i & 1) << 1)
            for i in range(max(value[0].bit_length(), value[1].bit_length()))]


def reference_multiply(a, b):
    aa, bb = coefficients(a), coefficients(b)
    out = [0] * max(0, len(aa) + len(bb) - 1)
    for i, x in enumerate(aa):
        for j, y in enumerate(bb):
            out[i + j] ^= F4_PRODUCT[x][y]
    return tuple(sum(((value >> plane) & 1) << i for i, value in enumerate(out))
                 for plane in (0, 1))


def reference_compress(seed, low, a, b):
    """Literal four-leg contraction, with coefficientwise F4 convolution."""
    ambient, q = low + 2, low + 1
    ab = reference_multiply(a, b)
    maps = (((1, 0), a), ((1 ^ ab[0], ab[1]), b))
    out = []
    for reverse in product(range(q), repeat=4):
        logical_labels = reverse[::-1]
        support = [[(label, (1, 0))] if label < low else
                   [(low + i, value) for i, value in enumerate(maps[axis // 2])
                    if value != (0, 0)] for axis, label in enumerate(logical_labels)]
        answer = [0, 0]
        for term in product(*support):
            physical = [item[0] for item in term]
            value = seed['values'][sum(label * ambient**axis
                                       for axis, label in enumerate(physical))]
            for _, scalar in term:
                value = reference_multiply(value, scalar)
            answer[0] ^= value[0]
            answer[1] ^= value[1]
        out.append(tuple(answer))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    random = Random(1952001)
    arithmetic = 0
    for _ in range(500):
        a, b = ((random.randrange(1 << 12), random.randrange(1 << 12)) for _ in range(2))
        assert multiply(a, b) == reference_multiply(a, b)
        left, right = paired_maps(a, b)
        p, r = reference_multiply(left[0], right[0]), reference_multiply(left[1], right[1])
        assert (p[0] ^ r[0], p[1] ^ r[1]) == (1, 0)
        arithmetic += 1
    lower = [0, 1, 2, 3, 4, 5, 8]
    raw = restricted_tensor(16, list(map(logical, lower + [6, 9])))
    constant_entries = polynomial_entries = endpoints = 0
    profiles = []
    for av, bv in product(range(4), repeat=2):
        a, b = (av & 1, av >> 1), (bv & 1, bv >> 1)
        actual = compress(raw, 7, a, b)
        left, right = paired_maps(a, b)
        expected = constant_compress(raw, 7,
                                    [[evaluate_at_one(value) for value in left]],
                                    [[evaluate_at_one(value) for value in right]])
        assert actual['values'] == expected['values']
        constant_entries += len(actual['values'])
    for a, b in (((2, 0), (1, 2)), ((3, 2), (4, 3)), ((0, 2), (3, 0))):
        actual = compress(raw, 7, a, b)
        assert actual['values'] == reference_compress(raw, 7, a, b)
        assert actual['degree'] <= raw['degree'] + 4 * max(len(coefficients(a)) - 1, 0) + 2 * max(len(coefficients(b)) - 1, 0)
        for index, value in enumerate(actual['values']):
            x, y, u, v = ((index // 8**axis) % 8 for axis in range(4))
            assert evaluate_at_one(value) == int(x == v and y == u)
            endpoints += 1
        polynomial_entries += len(actual['values'])
        sides = natural_profiles(actual, [0, 4, 12, 16, 20, 24, 28, 32])
        assert all(len(side) == 64 for side in sides)
        profiles.append({'a': a, 'b': b, 'degree': actual['degree'],
                         'tensor_sha256': digest(actual['values']), 'orders_infinity_zero': sides})
    source_paths = ['tools/moving_boundary_crossings.py', 'tools/check_moving_boundary_crossings.py']
    report = {'scope': 'Exact polynomial contraction and complete natural-profile controls; no improved factoring bound.',
              'independent_arithmetic_and_split_controls': arithmetic,
              'constant_specialization_entries': constant_entries,
              'independent_full_polynomial_entries': polynomial_entries,
              'separately_checked_endpoint_entries': endpoints,
              'profiles': profiles,
              'source_sha256': {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in source_paths}}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in report.items() if isinstance(v, (str, int))}))


if __name__ == '__main__':
    main()
