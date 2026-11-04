#!/usr/bin/env python3
"""Exact fixed-curve controls for the ordinary colour-alternation lemma.

No full tensor profile, alphabet sweep or numerical pressure minimization.
Source construction: cube_homogeneous_gaussian and PR97 charge fusion.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb
from pathlib import Path

import cube_homogeneous_gaussian as cube
from check_quartic_charge_fusion import ELL2, QUARTIC, binary_rank

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/colour-alternation-joint-cold-curve.json'
LEFT, RIGHT = (0, 2, 4, 6), (7, 5, 3, 1)
WITNESSES = (((0,), (7,)), ((0, 4), (7, 3)),
             ((0, 2, 4), (7, 5, 3)), (LEFT, RIGHT))


def curve(poly):
    out = set()
    for a, b, c in poly:
        term = (2*a+c, 2*b+c, 0)
        if term in out:
            out.remove(term)
        else:
            out.add(term)
    return frozenset(out)


KERNEL = [[curve(cube.multiply(ELL2, entry)) for entry in row]
          for row in cube.kernel()]
H = curve(QUARTIC)


@lru_cache(None)
def gaussian(word):
    """Direct all-pairing expansion, independent of determinant expansion."""
    if not word:
        return cube.ONE
    (v, color), rest = word[0], word[1:]
    out = cube.ZERO
    for i, (w, other) in enumerate(rest):
        if color == other:
            out ^= cube.multiply(KERNEL[v][w],
                                 gaussian(rest[:i]+rest[i+1:]))
    return out


def valid_quartic(word):
    colors = [c for _, c in word]
    charge = 0
    for c in colors:
        charge ^= c
    return len(set(colors)) == 4 and charge == 0


def tensor_terms(word):
    """Separate Gaussian, one-quartic and two-quartic contributions."""
    out = [gaussian(word), cube.ZERO, cube.ZERO]
    for support in combinations(range(len(word)), 4):
        if valid_quartic(tuple(word[i] for i in support)):
            rest = tuple(item for i, item in enumerate(word) if i not in support)
            out[1] ^= cube.multiply(H, gaussian(rest))
    if len(word) == 8:
        # Each unordered split once, by requiring the first vertex in part one.
        for tail in combinations(range(1, 8), 3):
            support = (0,)+tail
            first = tuple(word[i] for i in support)
            rest = tuple(item for i, item in enumerate(word) if i not in support)
            if valid_quartic(first) and valid_quartic(rest):
                out[2] ^= cube.multiply(H, H)
    return out


def scalars_and_minors():
    assert all(x+y == 16 and z == 0 for row in KERNEL for p in row for x, y, z in p)
    assert all(x+y == 32 and z == 0 for x, y, z in H)
    a = frozenset(((6, 10, 0), (8, 8, 0), (10, 6, 0)))
    assert KERNEL[0][3] == KERNEL[4][7] == cube.ZERO
    assert KERNEL[0][7] == KERNEL[4][3] == a
    assert H == frozenset(((20, 12, 0), (17, 15, 0),
                           (15, 17, 0), (12, 20, 0)))
    one, b, c = cube.ONE, frozenset(((0, 0, 0), (3, 0, 0))), \
        frozenset(((3, 0, 0), (6, 0, 0)))
    n = [[one, b, cube.ZERO, c], [b, one, c, cube.ZERO],
         [cube.ZERO, c, one, b], [c, cube.ZERO, b, one]]
    unit = frozenset(((0, 0, 0), (2, 0, 0), (4, 0, 0)))
    for i, v in enumerate(LEFT):
        for j, w in enumerate(RIGHT):
            specialized = frozenset((x, 0, 0) for x, _, _ in KERNEL[v][w])
            assert specialized == cube.shift(cube.multiply(unit, n[i][j]), (6, 0, 0))
    assert cube.matrix_product(n, n) == [
        [cube.monomial(12, 0, 0) if i == j else cube.ZERO for j in range(4)]
        for i in range(4)]
    assert binary_rank([sum(((0, 0, 0) in p) << j for j, p in enumerate(row))
                        for row in n]) == 2
    records = []
    for axis in (0, 1):
        # y cut is obtained by swapping spatial bits zero and one.
        swap = lambda v: (v & 4) | ((v & 1) << 1) | ((v & 2) >> 1)
        left, right = (LEFT, RIGHT) if axis == 0 else \
            (tuple(map(swap, LEFT)), tuple(map(swap, RIGHT)))
        cross = [[KERNEL[i][j] for j in right] for i in left]
        minima = []
        for r in range(1, 5):
            values = [cube.minor(cross, rows, cols)
                      for rows in combinations(range(4), r)
                      for cols in combinations(range(4), r)]
            minima.append(min(cube.valuation(p, axis) for p in values if p))
        assert minima == [6, 12, 30, 48]
        records.append(dict(axis=axis, all_minors=69, minimum_valuations=minima))
    for (rows, cols), val in zip(WITNESSES, (6, 12, 30, 48)):
        assert cube.valuation(cube.minor(KERNEL, rows, cols), 0) == val
    return dict(cut_controls=records, one_particle_orders=[6, 6, 18, 18],
                kernel_square_order=12, constant_rank=2)


def alternation_controls():
    controls = [
        ((1,), (1,)), ((1,), (2,)),
        ((1, 2), (1, 2)), ((1, 2), (4, 7)),
        ((1, 2, 4), (1, 2, 4)), ((1, 2, 4), (1, 3, 5)),
        ((1, 2, 4, 7), (1, 2, 4, 7)),
        ((1, 2, 4, 7), (8, 9, 10, 11)),
    ]
    records = []
    for row_colors, col_colors in controls:
        rows, cols = WITNESSES[len(row_colors)-1]
        total = [cube.ZERO]*3
        raw_nonzero = [0]*3
        for p in permutations(col_colors):
            terms = tensor_terms(tuple(zip(rows, row_colors))+tuple(zip(cols, p)))
            for k, value in enumerate(terms):
                total[k] ^= value
                raw_nonzero[k] += bool(value)
        expected = cube.minor(KERNEL, rows, cols) if row_colors == col_colors else cube.ZERO
        assert total == [expected, cube.ZERO, cube.ZERO]
        records.append(dict(row_colors=list(row_colors), column_colors=list(col_colors),
                            raw_nonzero_by_quartic_count=raw_nonzero,
                            alternated_quartic_terms_zero=True))
    assert records[3]['raw_nonzero_by_quartic_count'][1] > 0
    assert records[5]['raw_nonzero_by_quartic_count'][1] > 0
    assert records[-1]['raw_nonzero_by_quartic_count'][2] > 0
    # Upper occupation blocks: a genuine raw quartic term must cancel.
    rows, cols = (0, 4), RIGHT
    total = [cube.ZERO]*3
    nonzero = 0
    for p in permutations((1, 4, 8, 14)):
        terms = tensor_terms(tuple(zip(rows, (1, 2)))+tuple(zip(cols, p)))
        nonzero += bool(terms[1])
        for k, value in enumerate(terms):
            total[k] ^= value
    assert nonzero and total == [cube.ZERO]*3
    return dict(equal_occupation=records, upper_block_raw_nonzero=nonzero,
                upper_block_alternation_zero=True)


def pressure_control():
    # Exact exponents for the three specified AM--GM comparisons.
    weights = [Fraction(2, 9), Fraction(20, 27), Fraction(4, 9)]
    orders = [(6, 12), (12, 48), (12, 48)]
    d_powers = [(1, 2), (2, 4), (2, 4)]
    assert [w*a+(1-w)*b for w, (a, b) in zip(weights, orders)] == \
        [Fraction(32, 3), Fraction(64, 3), Fraction(32)]
    assert [w*a+(1-w)*b for w, (a, b) in zip(weights, d_powers)] == \
        [Fraction(16, 9), Fraction(68, 27), Fraction(28, 9)]
    d0 = 1 << 90
    assert d0**7 >= 36**9 and d0**14 >= 1152**27 and d0 >= 384**9
    assert Fraction(9, 16) > Fraction(1, 2)
    # Endpoint d=6 validates the elementary coefficient comparison guards.
    assert comb(6, 2) >= Fraction(6**2, 4)
    assert comb(6, 4) >= Fraction(6**4, 192)
    return dict(sufficient_d_min=str(d0), rho_floor='3/4',
                method='Exact AM-GM and entropy inequalities; no parameter search')


def run():
    names = ['tools/check_colour_alternation.py', 'tools/cube_homogeneous_gaussian.py',
             'tools/check_quartic_charge_fusion.py']
    return dict(scope='Fixed-curve structural prefix bounds; no full tensor profile or complete factoring improvement.',
                scalars_and_minors=scalars_and_minors(),
                alternation=alternation_controls(), pressure=pressure_control(),
                source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(RECEIPT.read_text()), 'receipt differs; inspect before regeneration'
    print('Colour alternation: exact joint-cold minors, quartic cancellation and pressure PASS')
