#!/usr/bin/env python3
"""Independent dense polynomial controls for codimension-one Laurent frames."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from check_coherent_hot_boundary import logical, natural_profiles
from check_moving_boundary_crossings import reference_multiply
from check_smith_profile_transfer import smith_orders
from moving_boundary_crossings import compress, evaluate_at_one
from moving_physical_frames import compress_frame
from physical_split_crossings import flatten
from restrict_implicit_reflection import restricted_tensor

ROOT = Path(__file__).resolve().parents[1]


def add(a, b):
    return a[0] ^ b[0], a[1] ^ b[1]


def shift(value, power):
    if value == (0, 0):
        return value
    assert power >= 0
    return value[0] << power, value[1] << power


def reference_maps(a, b, sa, sb, right_a=None):
    q = len(a)
    right_a = a if right_a is None else right_a
    ml = max([0] + [-s for value, s in zip(a, sa) if value != (0, 0)])
    mc = max([0] + [-s for value, s in zip(right_a, sa) if value != (0, 0)])
    mp = max([0] + [mc - s for value, s in zip(b, sb) if value != (0, 0)])
    left = [[(int(i == j) << ml, 0) for j in range(q)] + [shift(a[i], ml + sa[i])]
            for i in range(q)]
    right = [[add((int(i == j) << mp, 0),
                  shift(reference_multiply(right_a[i], b[j]), mp + sa[i] + sb[j]))
              for i in range(q)] + [shift(b[j], mp + sb[j])] for j in range(q)]
    for i, j in product(range(q), repeat=2):
        value = (0, 0)
        for k in range(q + 1):
            value = add(value, reference_multiply(left[i][k], right[j][k]))
        if right_a == a:
            assert value == (int(i == j) << (ml + mp), 0)
        else:
            assert evaluate_at_one(value) == int(i == j)
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


def independent_profiles(seed, energies):
    result = []
    for cut in (1, 0):
        matrix, rows, columns = flatten(seed, cut)
        weights = [[((energies[b] + energies[d]) // 2 if cut else
                     (energies[a] + energies[b]) // 2)
                    for c, d in columns] for a, b in rows]
        assert all(high == 0 for row in matrix for low, high in row)
        common = max(0, max(low.bit_length() - 1 - weights[i][j]
                            for i, row in enumerate(matrix) for j, (low, _) in enumerate(row))) if cut else 0
        binary = []
        for i, row in enumerate(matrix):
            output = []
            for j, (low, _) in enumerate(row):
                if cut:
                    degree = common + weights[i][j]
                    value = sum(1 << (degree - k) for k in range(low.bit_length()) if low >> k & 1)
                else:
                    value = low << weights[i][j]
                output.append(value)
            binary.append(output)
        degree = max(value.bit_length() - 1 for row in binary for value in row)
        result.append([order - common for order in smith_orders(binary, degree, precision=512)])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = Random(1952441)
    dense_entries = endpoint_entries = 0
    for trial in range(48):
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
        right_a = a
        if trial >= 24:
            sa, sb = [0, 0], [0, 0]
            scalars = [rng.randrange(4) for _ in range(2)]
            if not any(scalars):
                scalars[0] = 1
            right_a = [add(value, reference_multiply((3, 0), (c & 1, c >> 1)))
                       for value, c in zip(a, scalars)]
        left, right, clearing = reference_maps(a, b, sa, sb, right_a)
        actual = compress_frame(seed, a, b, sa, sb, right_a=right_a)
        assert actual['laurent_clearing_power'] == clearing
        assert actual['values'] == reference_contract(seed, left, right)
        if trial < 24:
            dense_entries += len(actual['values'])
        else:
            endpoint_entries += len(actual['values'])
    raw = restricted_tensor(16, list(map(logical, [0, 1, 2, 3, 4, 5, 8, 6, 9])))
    specialization_entries = 0
    for av, bv in (((1, 0), (1, 0)), ((3, 0), (2, 0)), ((2, 1), (1, 2))):
        a, b = [(0, 0)] * 7 + [av], [(0, 0)] * 7 + [bv]
        actual = compress_frame(raw, a, b)
        assert actual['values'] == compress(raw, 7, av, bv)['values']
        specialization_entries += len(actual['values'])
    energies = [0, 4, 12, 16, 20, 24, 28, 32]
    shifts = [(32 - value) // 4 for value in energies]
    profile_controls = []
    for am, bm, endpoint_only in ((128, 128, False), (170, 85, False),
                                 (255, 255, False), (170, 85, True)):
        a = [(int(am >> i & 1), 0) for i in range(8)]
        b = [(int(bm >> i & 1), 0) for i in range(8)]
        if endpoint_only:
            a = [(int(i == 7) | (2 if am >> i & 1 else 0), 0) for i in range(8)]
            r = [(2 if (am ^ 128) >> i & 1 else 0, 0) for i in range(8)]
            seed = compress_frame(raw, a, b, right_a=r)
        else:
            seed = compress_frame(raw, a, b, shifts, [-value for value in shifts])
        independent = independent_profiles(seed, energies)
        assert independent == natural_profiles(seed, energies)
        assert all(len(side) == 64 for side in independent)
        profile_controls.append({'a_mask': am, 'b_mask': bm, 'endpoint_only': endpoint_only,
                                 'orders_infinity_zero': independent})
    names = ['tools/check_moving_physical_frames.py', 'tools/moving_physical_frames.py',
             'tools/moving_boundary_crossings.py', 'tools/check_moving_boundary_crossings.py']
    report = {'scope': 'Literal dense F4 polynomial contractions, 24 exact Laurent splits, 24 endpoint-only splits, and independent complete binary Smith profiles; no improved factoring bound.',
              'dense_laurent_polynomial_entries': dense_entries,
              'dense_endpoint_only_polynomial_entries': endpoint_entries,
              'full_reflection_specialization_entries': specialization_entries,
              'random_dense_frames': 24,
              'random_endpoint_only_frames': 24,
              'independent_profile_orders': 512,
              'independent_profile_controls': profile_controls,
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in names}}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: value for key, value in report.items() if isinstance(value, (str, int))}))


if __name__ == '__main__':
    main()
