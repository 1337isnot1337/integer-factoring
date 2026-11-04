#!/usr/bin/env python3
"""Finite q8 controls for lower-valuation determinant-shaped quartics.

The spatial-mode quartic uses exponents (12,20), (16,16), (20,12).
This script replaces the outer exponents by four lower positive
valuations while keeping degree 32, the middle exponent 16, the same
alpha support rule and the even endpoint cancellation.  It recomputes
the complete low Smith prefix on both cuts and checks that order twelve
stays 35.  This is a finite coefficient control, not a lower-valuation
theorem or a profile improvement.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path

import check_fixed_curve_quartic_floor as fixed
import check_spatial_mode_colour_fusion as spatial

ROOT = Path(__file__).resolve().parents[1]

CANONICAL = Counter({0: 1, 6: 14, 12: 35, 15: 84, 16: 70, 17: 14, 18: 208})
EXPECTED = {
    'canonical_12': CANONICAL,
    'valuation_6': Counter({0: 1, 6: 70, 12: 283, 15: 84, 16: 98, 18: 184}),
    'valuation_9': Counter({0: 1, 6: 14, 9: 56, 12: 35, 15: 248, 16: 98, 18: 86}),
    'valuation_10': Counter({0: 1, 6: 14, 10: 56, 12: 35, 16: 346, 18: 26}),
    'valuation_11': Counter({0: 1, 6: 14, 11: 56, 12: 35, 16: 98, 17: 248, 18: 26}),
}
PARAMETERS = {
    'canonical_12': (12, 20, 16),
    'valuation_6': (6, 26, 16),
    'valuation_9': (9, 23, 16),
    'valuation_10': (10, 22, 16),
    'valuation_11': (11, 21, 16),
}


def make_coefficient(a, b, middle):
    def coefficient(word, axis):
        ax, ay = spatial.alphas(word)
        terms = set()
        for enabled, exp in ((ax, (a, b)), (ax ^ ay, (middle, middle)), (ay, (b, a))):
            if enabled:
                terms.add(exp)
        assert all(x + y == 32 and x > 0 and y > 0 for x, y in terms)
        assert len(terms) % 2 == 0
        value = 0
        for exp in terms:
            value ^= 1 << exp[axis]
        return value & fixed.MASK
    return coefficient


def prefix(axis, coefficient):
    original = spatial.coefficient
    spatial.coefficient = coefficient
    try:
        combined = Counter()
        for target in range(8):
            result = spatial.matrix(axis, target)
            combined.update(fixed.ranks(result)['smith_prefix'])
    finally:
        spatial.coefficient = original
    return combined


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cases = []
    for axis in (0, 1):
        for name, (a, b, middle) in PARAMETERS.items():
            result = prefix(axis, make_coefficient(a, b, middle))
            expected = EXPECTED[name]
            assert result == expected, (axis, name, result)
            cases.append({
                'axis': axis,
                'coefficient': name,
                'exponents': [a, b, middle],
                'smith_prefix': {str(k): v for k, v in sorted(result.items())},
                'order_6': result.get(6, 0),
                'order_9': result.get(9, 0),
                'order_10': result.get(10, 0),
                'order_11': result.get(11, 0),
                'order_12': result.get(12, 0),
            })
    payload = {
        'scope': ('Finite q8 controls for four lower-valuation determinant-'
                  'shaped quartic coefficients; no general lower-valuation '
                  'theorem or pressure improvement'),
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
