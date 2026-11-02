#!/usr/bin/env python3
"""Show that the height-38 protected-core barrier does not extend to height 39."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from coherent_hot_compounds import height

A = 70
B = 99
Z_POINT = Fraction(899, 900)
CAPS = (37, 38, 39, 40)
EXPECTED_CAPACITIES = {37: 372, 38: 409, 39: 450, 40: 496}
EXPECTED_SIGNS = {37: True, 38: True, 39: False, 40: False}


def q_coefficients(cap: int) -> Counter[int]:
    q = Counter({0: 1})
    for odd in range(1, cap + 1, 2):
        for height, count in list(q.items()):
            if height + odd <= cap:
                q[height + odd] += count
    return q


def fraction_digest(value: Fraction) -> str:
    def packed(n: int) -> bytes:
        return n.to_bytes((n.bit_length() + 7) // 8, 'big')

    return sha256(packed(value.numerator) + b'\0' + packed(value.denominator)).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    core_path = Path('evidence/height24-core-low-orders.json')
    core = json.loads(core_path.read_text())
    assert core['completed'] is True
    assert core['words'] == [w for w in range(4096) if height(w) < 24]
    assert len(core['words']) == 81
    profile = core['profile']
    assert profile['q'] == 81 and profile['cutoff'] == 64
    orders = profile['cuts'][0]['retained_orders_t']
    assert orders == profile['cuts'][1]['retained_orders_t']
    assert len(orders) == 3359 and sum(orders) == 173663
    assert all(cut['tail_count'] == 3202 for cut in profile['cuts'])
    z_poly = Counter(orders)
    z_poly[82] += 3202
    assert sum(z_poly.values()) == 6561

    z = Z_POINT
    assert 0 < z < 1
    assert 2 * A * A < B * B
    rows = []
    for cap in CAPS:
        q = q_coefficients(cap)
        capacity = sum(q.values())
        assert capacity == EXPECTED_CAPACITIES[cap]
        left = sum(count * z ** (A * order) for order, count in z_poly.items())
        right = sum(count * z ** (B * mode_height) for mode_height, count in q.items())
        left **= B
        right **= 2 * A
        excludes = left > right
        assert excludes is EXPECTED_SIGNS[cap]
        rows.append({
            'maximum_height': cap,
            'ambient_height_inventory_q': capacity,
            'test_excludes': excludes,
            'left_bits': left.numerator.bit_length(),
            'right_bits': right.numerator.bit_length(),
            'left_sha256': fraction_digest(left),
            'right_sha256': fraction_digest(right),
        })

    source_names = sorted(core['source_sha256'])
    source_pins = {}
    for name in source_names:
        path = Path('tools') / name
        assert path.exists(), name
        assert sha256(path.read_bytes()).hexdigest() == core['source_sha256'][name], name
        source_pins[name] = core['source_sha256'][name]

    receipt = {
        'completed': True,
        'a': A,
        'b': B,
        'rho0': [A, B],
        'rho0_decimal': float(Fraction(A, B)),
        'rho0_below_inverse_sqrt_two': 2 * A * A < B * B,
        'z': [Z_POINT.numerator, Z_POINT.denominator],
        'test': 'Z(z^a)^b > Q_H(z^b)^(2a) excludes rho0; failure at one 0<z<1 does not exclude it.',
        'rows': rows,
        'height38_floor': [29, 41],
        'height38_floor_above_inverse_sqrt_two': 2 * 29 * 29 > 41 * 41,
        'core_receipt_sha256': sha256(core_path.read_bytes()).hexdigest(),
        'source_sha256': source_pins,
        'scope': (
            'Exact rational comparison at one interior point for the Jensen-tail '
            'necessary condition. It shows that the height-38 lower-bound argument '
            'does not exclude rho0=70/99 at height 39; it does not construct a split, '
            'prove a uniform pressure bound, or change the factoring exponent.'
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({
        'completed': True,
        'rho0': [A, B],
        'below_inverse_sqrt_two': receipt['rho0_below_inverse_sqrt_two'],
        'signs': {row['maximum_height']: 'excluded' if row['test_excludes'] else 'admitted' for row in rows},
    }))


if __name__ == '__main__':
    main()
