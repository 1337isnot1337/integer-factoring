#!/usr/bin/env python3
"""Exact controls for the order-twelve floor under asymmetric quartic weights.

The proof is in research/colour-asymmetric-order12-floor.md.  This script
recomputes the q8 spatial-mode charge blocks for both cuts.  It compares
the canonical completion with the Gaussian-only block and deterministic
asymmetric subsets of the distinct-colour quartic gates.  The order-twelve
Smith increment must stay at least 4 on every nonzero charge and at least
7 on charge zero.
"""

from __future__ import annotations

import argparse
import json
from hashlib import sha256
from itertools import product
from pathlib import Path
import random

import check_fixed_curve_quartic_floor as s
import check_spatial_mode_colour_fusion as spatial

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/colour-asymmetric-order12-floor.json'


def prefix(axis, target, enabled=None):
    original = spatial.coefficient
    if enabled is not None:
        def coefficient(word, axis, _original=original, _enabled=enabled):
            return _original(word, axis) if _enabled(word, axis) else 0
        spatial.coefficient = coefficient
    try:
        blocks = spatial.matrix(axis, target)
    finally:
        spatial.coefficient = original
    return s.ranks(blocks)


def enabled_pattern(seed):
    rng = random.Random(seed)
    return lambda word, axis, _rng=rng: _rng.randrange(2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cases = []
    for axis in (0, 1):
        for charge in range(8):
            expected_six = 0 if charge == 0 else 2
            records = {}
            for name, enabled in (
                ('canonical', None),
                ('gaussian_only', lambda word, axis: 0),
                ('asymmetric_0', enabled_pattern(0x51554152 + 17*axis + charge)),
                ('asymmetric_1', enabled_pattern(0x54494353 + 29*axis + charge)),
            ):
                result = prefix(axis, charge, enabled)
                records[name] = {
                    'smith_prefix': {str(k): v for k, v in result['smith_prefix'].items()},
                    'order_6': result['smith_prefix'].get(6, 0),
                    'order_12': result['smith_prefix'].get(12, 0),
                }
            canonical = records['canonical']
            gaussian = records['gaussian_only']
            assert canonical['order_6'] == expected_six
            assert canonical['order_12'] == (7 if charge == 0 else 4)
            assert gaussian['order_6'] == expected_six
            assert gaussian['order_12'] == (7 if charge == 0 else 12)
            for name in ('asymmetric_0', 'asymmetric_1'):
                assert records[name]['order_6'] == expected_six
                assert records[name]['order_12'] >= (7 if charge == 0 else 4)
            cases.append({'axis': axis, 'charge': charge, 'records': records})
    payload = {
        'scope': ('Finite q8 controls for the order-twelve floor under '
                  'asymmetric distinct-colour quartic weights; no profile '
                  'or factoring exponent is recomputed'),
        'cases': cases,
        'nonzero_order_12_floor': 4,
        'zero_order_12': 7,
        'canonical_minimal_nonzero_order_12': True,
    }
    raw = json.dumps(payload, separators=(',', ':'), sort_keys=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(raw + '\n')
    print(json.dumps({
        'cases': len(cases),
        'nonzero_order_12_floor': 4,
        'zero_order_12': 7,
        'canonical_minimal': True,
        'sha256': sha256(raw.encode()).hexdigest(),
    }, sort_keys=True))


if __name__ == '__main__':
    main()
