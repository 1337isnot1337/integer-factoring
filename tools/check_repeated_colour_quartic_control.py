#!/usr/bin/env python3
"""Finite q8 control for repeated-colour quartic gates.

A four-vertex gate on a support with two colours repeated twice can act on
the diagonal P=Q block, unlike a distinct-colour quartic.  This script tests
the strongest local cancellation form: replace the pair-product leading term
by its endpoint-safe shift pair_product*(1+x), on all total-four supports or
only on the repeated two-colour supports.  It recomputes the complete low
Smith prefix on both cuts.  The order-twelve channel never falls below the
canonical 35.  This is a finite control, not a repeated-colour theorem; the
q>=34 parity obstruction still applies to parity-preserving repeated gates.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from hashlib import sha256
from itertools import permutations
from pathlib import Path

import check_fixed_curve_quartic_floor as fixed
import check_spatial_mode_colour_fusion as spatial

ROOT = Path(__file__).resolve().parents[1]
MASK = fixed.MASK

CANONICAL = Counter({0: 1, 6: 14, 12: 35, 15: 84, 16: 70, 17: 14, 18: 208})
EXPECTED = {
    'baseline': CANONICAL,
    'repeated_cancel': Counter({0: 1, 6: 14, 12: 35, 14: 56, 16: 98, 17: 10, 18: 208}),
    'all_cancel': Counter({0: 1, 6: 14, 12: 91, 18: 308}),
}


def pair_product(aa, bb, pair):
    if not aa or not bb or len(aa) != len(bb):
        return 0
    value = 0
    for perm in permutations(bb):
        if any(c != d for (_, c), (_, d) in zip(aa, perm)):
            continue
        term = 1
        for (i, _), (j, _) in zip(aa, perm):
            term = fixed.pmul(term, pair[i][j])
        value ^= term
    return value


def make_entry(mode):
    original = spatial.entry

    def entry(aa, bb, left, right, pair, axis):
        if mode == 'baseline':
            return original(aa, bb, left, right, pair, axis)
        total = len(aa) + len(bb)
        if not total:
            return 1
        if mode == 'repeated_cancel' and total != 4:
            return original(aa, bb, left, right, pair, axis)
        out = pair_product(aa, bb, pair)
        if total != 4:
            return out
        colors = [c for _, c in aa + bb]
        repeated = len(set(colors)) == 2 and all(colors.count(c) == 2 for c in set(colors))
        if mode == 'repeated_cancel':
            out = original(aa, bb, left, right, pair, axis)
            if repeated:
                out ^= (pair_product(aa, bb, pair) << 2) & MASK
        elif mode == 'all_cancel':
            out ^= (pair_product(aa, bb, pair) << 2) & MASK
        return out

    return entry


def prefix(axis, mode):
    original = spatial.entry
    spatial.entry = make_entry(mode)
    try:
        combined = Counter()
        for target in range(8):
            combined.update(fixed.ranks(spatial.matrix(axis, target))['smith_prefix'])
    finally:
        spatial.entry = original
    return combined


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cases = []
    for axis in (0, 1):
        for name in ('baseline', 'repeated_cancel', 'all_cancel'):
            result = prefix(axis, name)
            expected = EXPECTED[name]
            assert result == expected, (axis, name, result)
            cases.append({
                'axis': axis,
                'control': name,
                'smith_prefix': {str(k): v for k, v in sorted(result.items())},
                'order_12': result.get(12, 0),
            })
    payload = {
        'scope': ('Finite q8 repeated-colour cancellation controls; no general '
                  'repeated-colour theorem, profile improvement or transfer claim'),
        'canonical_order_12': 35,
        'minimum_order_12': min(c['order_12'] for c in cases),
        'cases': cases,
    }
    raw = json.dumps(payload, separators=(',', ':'), sort_keys=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(raw + '\n')
    print(json.dumps({
        'cases': len(cases),
        'minimum_order_12': payload['minimum_order_12'],
        'canonical_order_12': 35,
        'sha256': sha256(raw.encode()).hexdigest(),
    }, sort_keys=True))


if __name__ == '__main__':
    main()
