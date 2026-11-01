#!/usr/bin/env python3
"""Independent dense polynomial controls for codimension-one Laurent frames."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from check_coherent_hot_boundary import logical
from check_moving_boundary_crossings import reference_multiply
from moving_boundary_crossings import compress
from moving_physical_frames import compress_frame
from restrict_implicit_reflection import restricted_tensor

ROOT = Path(__file__).resolve().parents[1]


def add(a, b):
    return a[0] ^ b[0], a[1] ^ b[1]


def shift(value, power):
    if value == (0, 0):
        return value
    assert power >= 0
    return value[0] << power, value[1] << power


def reference_maps(a, b, sa, sb):
    q = len(a)
    ml = max([0] + [-s for value, s in zip(a, sa) if value != (0, 0)])
    mp = max([0] + [ml - s for value, s in zip(b, sb) if value != (0, 0)])
    left = [[(int(i == j) << ml, 0) for j in range(q)] + [shift(a[i], ml + sa[i])]
            for i in range(q)]
    right = [[add((int(i == j) << mp, 0),
                  shift(reference_multiply(a[i], b[j]), mp + sa[i] + sb[j]))
              for i in range(q)] + [shift(b[j], mp + sb[j])] for j in range(q)]
    for i, j in product(range(q), repeat=2):
        value = (0, 0)
        for k in range(q + 1):
            value = add(value, reference_multiply(left[i][k], right[j][k]))
        assert value == (int(i == j) << (ml + mp), 0)
    return left, right, 2 * (ml + mp)


def reference_contract(seed, left, right):
    q, ambient = len(left), len(left) + 1
    result = []
    for reverse in product(range(q), repeat=4):
        labels = reverse[::-1]
        value = (0, 0)
        for physical in product(range(ambient), repeat=4):
            term = seed['values'][sum(label * ambient**axis for axis, label in enumerate(physical))]
            for axis, matrix in enumerate((left, left, right, right)):
                term = reference_multiply(term, matrix[labels[axis]][physical[axis]])
            value = add(value, term)
        result.append(value)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = Random(1952441)
    dense_entries = 0
    for _ in range(24):
        values = []
        for reverse in product(range(3), repeat=4):
            x, y, u, v = reverse[::-1]
            low, high = rng.randrange(16), rng.randrange(16)
            low ^= (low.bit_count() & 1) ^ int(x == v and y == u)
            high ^= high.bit_count() & 1
            values.append((low, high))
        seed = {'values': values, 'dimensions': [3] * 4, 'q_x': 3, 'q_u': 3, 'degree': 3}
        a, b = [[(rng.randrange(4), rng.randrange(4)) for _ in range(2)] for _ in range(2)]
        sa, sb = [[rng.randrange(-2, 3) for _ in range(2)] for _ in range(2)]
        left, right, clearing = reference_maps(a, b, sa, sb)
        actual = compress_frame(seed, a, b, sa, sb)
        assert actual['laurent_clearing_power'] == clearing
        assert actual['values'] == reference_contract(seed, left, right)
        dense_entries += len(actual['values'])
    raw = restricted_tensor(16, list(map(logical, [0, 1, 2, 3, 4, 5, 8, 6, 9])))
    specialization_entries = 0
    for av, bv in (((1, 0), (1, 0)), ((3, 0), (2, 0)), ((2, 1), (1, 2))):
        a, b = [(0, 0)] * 7 + [av], [(0, 0)] * 7 + [bv]
        actual = compress_frame(raw, a, b)
        assert actual['values'] == compress(raw, 7, av, bv)['values']
        specialization_entries += len(actual['values'])
    names = ['tools/check_moving_physical_frames.py', 'tools/moving_physical_frames.py',
             'tools/moving_boundary_crossings.py', 'tools/check_moving_boundary_crossings.py']
    report = {'scope': 'Literal dense F4 polynomial contractions and exact scaled split identities for Laurent moving frames; no improved factoring bound.',
              'dense_laurent_polynomial_entries': dense_entries,
              'full_reflection_specialization_entries': specialization_entries,
              'random_dense_frames': 24,
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in names}}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: value for key, value in report.items() if isinstance(value, (str, int))}))


if __name__ == '__main__':
    main()
