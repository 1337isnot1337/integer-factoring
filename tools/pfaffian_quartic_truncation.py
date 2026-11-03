"""One fixed repeated-colour, degree-32 quartic completion on the q8 curve.

Kernel and exact rank scaffolding retain PR111/PR114 credit. This constructor
changes their quartic coefficients, not the Gaussian pair kernel or wiring.
The jet evaluator is complete through order eighteen only; two quartics
start at twenty-four and are not assumed to vanish.
"""
from collections import Counter
from functools import cache
from itertools import permutations, product

import check_fixed_curve_quartic_floor as base

assert base.PRECISION == 19


def multiply(a, b):
    """Untruncated binary polynomial product."""
    out = 0
    while b:
        bit = b & -b
        out ^= a << (bit.bit_length()-1)
        b ^= bit
    return out


@cache
def kernel(axis):
    rows = []
    for row in base.kernel:
        values = []
        for polynomial in row:
            value = 0
            for a, b, c in polynomial:
                value ^= 1 << (2*(a if axis == 0 else b)+c)
            values.append(value)
        rows.append(tuple(values))
    return tuple(rows)


@cache
def sigma_jet(word, axis):
    """Cross-Gaussian four-occupied coefficient through order fifteen."""
    left = tuple((v, c) for v, c in word if not (v >> axis & 1))
    right = tuple((v, c) for v, c in word if v >> axis & 1)
    if len(left) != 2:
        return 0
    out = 0
    for ordered in permutations(right):
        if any(c != d for (_, c), (_, d) in zip(left, ordered)):
            continue
        value = 1
        for (v, _), (w, _) in zip(left, ordered):
            value = multiply(value, kernel(axis)[v][w])
        out ^= value
    return out & ((1 << 16)-1)


@cache
def quartic(word):
    """Return the homogeneous binary coefficient as its y=1 polynomial.

    A bit j represents x**j*y**(32-j). The low jets at the two fixed cuts
    are disjoint; the central coefficient restores the endpoint.
    """
    if any(count % 2 for count in Counter(c for _, c in word).values()):
        return 0
    lowx, lowy = (sigma_jet(word, axis) for axis in (0, 1))
    out = lowx
    for j in range(16):
        if lowy >> j & 1:
            out ^= 1 << (32-j)
    if (lowx.bit_count()+lowy.bit_count()) % 2:
        out ^= 1 << 16
    assert out.bit_count() % 2 == 0
    assert all(12 <= j <= 20 for j in range(33) if out >> j & 1)
    assert out & ((1 << 16)-1) == lowx
    reverse = sum(((out >> j) & 1) << (32-j) for j in range(33))
    assert reverse & ((1 << 16)-1) == lowy
    return out


@cache
def quartic_at_axis(word, axis):
    polynomial = quartic(word)
    if axis:
        polynomial = sum(((polynomial >> j) & 1) << (32-j) for j in range(33))
    return polynomial


def jet_entry(aa, bb, left, right, pair, axis):
    """Cut-normalized coefficient through18; only cross creation gates remain."""
    total = len(aa)+len(bb)
    if not total:
        return 1
    if not aa or not bb or total > 6:
        return 0
    out = 0
    if len(aa) == len(bb):
        for ordered in permutations(bb):
            if any(c != d for (_, c), (_, d) in zip(aa, ordered)):
                continue
            value = 1
            for (i, _), (j, _) in zip(aa, ordered):
                value = base.pmul(value, pair[i][j])
            out ^= value
    if total == 4:
        word = tuple(sorted([(left[i], c) for i, c in aa]+[(right[j], c) for j, c in bb]))
        out ^= quartic_at_axis(word, axis) & base.MASK
    if total == 6:
        for i, c in aa:
            for j, d in bb:
                if c != d:
                    continue
                word = tuple(sorted([(left[k], v) for k, v in aa if k != i]
                                    + [(right[k], v) for k, v in bb if k != j]))
                out ^= base.pmul(quartic_at_axis(word, axis) & base.MASK, pair[i][j])
    return out


def charge_block(axis, target):
    left, right = base.cube.axis_parts(axis)
    pair = [[kernel(axis)[v][w] & base.MASK for w in right] for v in left]
    words = [w for w in product(range(8), repeat=4) if base.charge(w) == target]
    assert len(words) == 512
    active = [tuple((i, c) for i, c in enumerate(w) if c) for w in words]
    blocks = [[0]*512 for _ in range(19)]
    six_controls = 0
    for i, aa in enumerate(active):
        for j, bb in enumerate(active):
            polynomial = jet_entry(aa, bb, left, right, pair, axis)
            if len(aa) == len(bb) == 3:
                assert not polynomial & (1 << 18)
                six_controls += 1
            while polynomial:
                bit = polynomial & -polynomial
                blocks[bit.bit_length()-1][i] |= 1 << j
                polynomial ^= bit
    return blocks, six_controls
