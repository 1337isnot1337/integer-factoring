#!/usr/bin/env python3
"""Five protected physical letters suffice for the natural 5/7 pressure floor."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from math import comb
import json
from pathlib import Path

from check_coherent_hot_boundary import logical, natural_profiles
from check_hot_boundary_core_barrier import integer_at, integer_hash
from check_moving_physical_frames import independent_profiles
from restrict_implicit_reflection import restricted_tensor

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    letters = [0, 1, 2, 3, 4]
    heights = [0, 1, 3, 4, 5, 6, 7, 8]
    raw = restricted_tensor(16, list(map(logical, letters)))
    profiles = independent_profiles(raw, [4 * h for h in heights[:5]])
    assert profiles == natural_profiles(raw, [4 * h for h in heights[:5]])
    orders = [(value - 128) // 2 for value in profiles[1]]
    assert profiles[1] == [128 + 2 * value for value in orders]
    assert profiles[0] == [-128 + 2 * value for value in orders]
    assert len(orders) == 25
    Z, Q = Counter(orders), Counter(heights)
    a, b, G = 5, 7, Counter()
    for i, ci in Z.items():
        for j, cj in Q.items():
            if i:
                G[a * (i - 1) + b * j] += i * ci * cj
            if j:
                G[a * i + b * (j - 1) + b - a] -= 2 * j * ci * cj
    G = {k: value for k, value in G.items() if value}
    degree = max(G)
    transformed = [sum(value * comb(degree - i, k - i)
                       for i, value in G.items() if i <= k) for k in range(degree + 1)]
    signs = [(value > 0) - (value < 0) for value in transformed]
    nonzero = [value for value in signs if value]
    variations = sum(x != y for x, y in zip(nonzero, nonzero[1:]))
    assert variations == 2 and G[0] == 1
    denominator, left, right = 10**8, 91056849, 91056850
    gl, gr = (integer_at(G, value, denominator) for value in (left, right))
    assert gl < 0 < gr
    zp, qp = {a * i: value for i, value in Z.items()}, {b * i: value for i, value in Q.items()}
    zn, qn = integer_at(zp, left, denominator), integer_at(qp, right, denominator)
    zd, qd = denominator ** max(zp), denominator ** max(qp)
    gap = zn**b * qd ** (2 * a) - qn ** (2 * a) * zd**b
    assert gap > 0 and 25**b > 8 ** (2 * a)
    names = ['tools/check_moving_five_core_floor.py', 'tools/check_moving_physical_frames.py',
             'tools/check_smith_profile_transfer.py', 'tools/check_hot_boundary_core_barrier.py',
             'tools/restrict_implicit_reflection.py']
    report = {'scope': 'Uniform natural common-scale 5/7 floor for arbitrary characteristic-two polynomial eight-letter tensors retaining the specified five-letter reflection subtensor. No degree limit, factoring lower bound, or exclusion of other physical gauges.',
              'status': 'passed', 'rho_floor': [a, b], 'protected_letters': letters,
              'candidate_heights': heights, 'core_orders_t': orders,
              'independent_binary_and_native_order_comparisons': 50,
              'derivative_degree': degree, 'derivative_coefficients': G,
              'mobius_signs': signs, 'mobius_sign_variations': variations,
              'minimum_box': [[left, denominator], [right, denominator]],
              'derivative_integer_hashes': [integer_hash(gl), integer_hash(gr)],
              'positive_box_gap_bits': gap.bit_length(), 'positive_box_gap_sha256': integer_hash(gap),
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in names}}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: report[key] for key in ('status', 'rho_floor', 'protected_letters', 'derivative_degree', 'mobius_sign_variations')}))


if __name__ == '__main__':
    main()
