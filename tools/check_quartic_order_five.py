#!/usr/bin/env python3
"""Exact lowest-coefficient and rational controls for the quartic obstruction."""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import argparse
import json
import cube_homogeneous_gaussian as cube
from check_quartic_charge_fusion import ELL2, binary_rank


def coefficient_controls():
    kernel = [[cube.multiply(ELL2, p) for p in row] for row in cube.kernel()]
    pairs = [p for p in combinations(range(1, 8), 2) if p[0] ^ p[1] == 1]
    words = []
    for colors in pairs:
        for i, j in permutations(range(4), 2):
            w = [0] * 4
            w[i], w[j] = colors
            words.append(tuple(w))
    records = []
    for axis in range(3):
        left, right = cube.axis_parts(axis)
        y, z = [a for a in range(3) if a != axis]
        exponent = [0] * 3
        exponent[axis], exponent[y], exponent[z] = 5, 9, 2
        target = tuple(exponent)
        rows = []
        for a in words:
            aa = [(i, c) for i, c in enumerate(a) if c]
            row = 0
            for j, b in enumerate(words):
                bb = [(i, c) for i, c in enumerate(b) if c]
                same = sorted(c for _, c in aa) == sorted(c for _, c in bb)
                p = cube.ONE if same else cube.ZERO
                if same:
                    for i, c in aa:
                        k = next(k for k, d in bb if c == d)
                        p = cube.multiply(p, kernel[left[i]][right[k]])
                assert all(e[z] >= 2 for e in p if e[axis] == 5)
                expected = same and (sum((left[i] >> z) & 1 for i, _ in aa)
                                     + sum((right[i] >> z) & 1 for i, _ in bb)) % 2
                assert (target in p) == bool(expected)
                row |= int(target in p) << j
            rows.append(row)
        projected = [sum((((rows[i] >> j) ^ (rows[i] & 1)
                            ^ (rows[0] >> j) ^ (rows[0] & 1)) & 1) << j
                         for j in range(1, len(words)))
                     for i in range(1, len(words))]
        assert binary_rank(rows) == 6 and binary_rank(projected) == 4
        # Every adjustable quartic monomial has larger second-axis order.
        for extra in range(3):
            e = [5, 5, 5]
            e[extra] += 1
            assert e[axis] != 5 or e[z] >= 5
        records.append(dict(axis=axis,word_count=len(words),
                            coefficient=list(target),raw_rank=6,projected_rank=4))
    return records


def pressure_controls():
    f = Fraction
    q = 4096
    assert f((q - 1) * (q - 4), q*q) >= f(31, 32)
    assert f(q - 1, q - 4) <= f(257, 256)
    assert (f(28*257, 9*256))**7 * f(3, 4)**29 < 1
    assert f(6**6, 5**5) * f(31, 64)**5 / 2 > f(3, 4)**6
    assert 3**33 > 2**47
    assert 3 * 7**2 > 12**2
    # Dimension/degree tail allocation; monotonicity proves all larger q.
    m = q*q - 2*q + 2
    assert 33*m <= q**4 and 2*m <= q**4
    return dict(q_at_least=q,rho_lower_bound=[7,12],tail_order=33,
                split_scale=[3,4],quadratic_allocation=True,
                young_allocation=True,large_scale_allocation=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = dict(scope='Formal coefficient controls and rational allocations for the ordinary restricted-family proof; no complete Smith profile or stronger factoring bound.',
                  coefficient_controls=coefficient_controls(),
                  pressure=pressure_controls(),
                  source_sha256={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in ('check_quartic_order_five.py',
                                              'check_quartic_charge_fusion.py',
                                              'cube_homogeneous_gaussian.py')})
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('Quartic order five: lowest coefficients and restored 7/12 allocations PASS')


if __name__ == '__main__':
    main()
