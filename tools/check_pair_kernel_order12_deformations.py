#!/usr/bin/env python3
"""Finite q8 controls for order-twelve pair-kernel deformations.

The structural discussion is in research/pair-kernel-order12-gauge.md.
The canonical spatial-mode prefix is fixed.  This script replaces the
order-six pair kernel by

    P + q u v^T,

where q is a fixed order-twelve, endpoint-zero polynomial and the support
matrix u v^T is one of seven deterministic patterns.  It recomputes the
complete low Smith prefix for both cuts and all eight charges.  The
control is finite and does not prove a uniform gauge or pressure bound.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path

import check_fixed_curve_quartic_floor as fixed
import check_spatial_mode_colour_fusion as spatial

ROOT = Path(__file__).resolve().parents[1]

CANONICAL = Counter({0: 1, 6: 14, 12: 35, 15: 84, 16: 70, 17: 14, 18: 208})
SHIFTED = Counter({0: 1, 6: 14, 12: 35, 15: 91, 16: 56, 17: 21, 18: 208})
SINGLE = Counter({0: 1, 6: 14, 12: 42, 15: 84, 16: 56, 17: 14, 18: 248})

PATTERNS = {
    'canonical': (None, None),
    'rankone_all_all': ((1, 1, 1, 1), (1, 1, 1, 1)),
    'rankone_r_all': ((1, 1, 0, 0), (1, 1, 1, 1)),
    'rankone_all_c': ((1, 1, 1, 1), (0, 0, 1, 1)),
    'single_e0_e2': ((1, 0, 0, 0), (0, 0, 1, 0)),
    'rankone_e0_c': ((1, 0, 0, 0), (0, 0, 1, 1)),
    'rankone_r_e2': ((1, 1, 0, 0), (0, 0, 1, 0)),
}

EXPECTED = {
    'canonical': CANONICAL,
    'rankone_all_all': CANONICAL,
    'rankone_r_all': CANONICAL,
    'rankone_all_c': CANONICAL,
    'single_e0_e2': SINGLE,
    'rankone_e0_c': SHIFTED,
    'rankone_r_e2': SHIFTED,
}

def support_entries(u, v):
    if u is None or v is None:
        return []
    return [(i, j) for i in range(4) for j in range(4) if u[i] and v[j]]


def pair_matrix(axis, support):
    left, right = fixed.cube.axis_parts(axis)

    def evaluate(poly):
        value = 0
        for a, b, c in poly:
            value ^= 1 << (2 * (a if axis == 0 else b) + c)
        return value & fixed.MASK

    pair = [[evaluate(fixed.kernel[i][j]) for j in right] for i in left]
    if not support:
        return left, right, pair
    p6 = pair[0][2]
    q = fixed.pmul(p6, p6)
    q ^= (q << 2) & fixed.MASK
    assert fixed.val(p6) == 6
    assert fixed.val(q) == 12
    assert q & 1 == 0
    for i, j in support:
        pair[i][j] ^= q
    return left, right, pair


def prefix(axis, support):
    left, right, pair = pair_matrix(axis, support)
    combined = Counter()
    for target in range(8):
        words = [w for w in product(range(8), repeat=4) if fixed.charge(w) == target]
        active = [tuple((i, c) for i, c in enumerate(w) if c) for w in words]
        blocks = [[0] * len(words) for _ in range(19)]
        for i, aa in enumerate(active):
            for j, bb in enumerate(active):
                p = spatial.entry(aa, bb, left, right, pair, axis)
                while p:
                    bit = p & -p
                    blocks[bit.bit_length() - 1][i] |= 1 << j
                    p ^= bit
        combined.update(fixed.ranks(blocks)['smith_prefix'])
    return combined


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cases = []
    for axis in (0, 1):
        for name, (u, v) in PATTERNS.items():
            support = support_entries(u, v)
            result = prefix(axis, support)
            expected = EXPECTED[name]
            assert result == expected, (axis, name, result)
            cases.append({
                'axis': axis,
                'deformation': name,
                'u': list(u) if u is not None else None,
                'v': list(v) if v is not None else None,
                'support': support,
                'smith_prefix': {str(k): v for k, v in sorted(result.items())},
                'order_6': result.get(6, 0),
                'order_12': result.get(12, 0),
            })
    payload = {
        'scope': ('Finite q8 controls for seven rank-one order-twelve '
                  'pair-kernel deformation patterns; no uniform gauge or '
                  'pressure theorem'),
        'row_space': [[1, 1, 0, 0], [0, 0, 1, 1]],
        'canonical_prefix': {str(k): v for k, v in sorted(CANONICAL.items())},
        'minimum_order_12': min(c['order_12'] for c in cases),
        'cases': cases,
    }
    raw = json.dumps(payload, separators=(',', ':'), sort_keys=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(raw + '\n')
    print(json.dumps({
        'cases': len(cases),
        'minimum_order_12': payload['minimum_order_12'],
        'canonical_order_12': CANONICAL[12],
        'sha256': sha256(raw.encode()).hexdigest(),
    }, sort_keys=True))


if __name__ == '__main__':
    main()
