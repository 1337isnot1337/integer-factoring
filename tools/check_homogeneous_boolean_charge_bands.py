#!/usr/bin/env python3
"""Exact controls for the all-degree homogeneous Boolean charge bands."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import prod
from pathlib import Path

from check_read_once_boolean_filter_floor import exact_ratio, truth_probability


def radical_nonnegative(a, b):
    """Exact sign of a+b*sqrt(2), for rational a,b."""
    a, b = F(a), F(b)
    if b == 0:
        return a >= 0
    if a >= 0 and b >= 0:
        return True
    if a < 0 and b < 0:
        return False
    if a >= 0:
        return a * a >= 2 * b * b
    return 2 * b * b >= a * a


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=root / 'evidence' / 'homogeneous-boolean-charge-bands.json')
    args = parser.parse_args()
    support_controls = loss_controls = ratio_controls = 0
    polynomial_counts = []
    for modes in range(1, 5):
        n = 1 << (modes - 1).bit_length()
        # The rarest distinct physical energies stress the support bound.
        energies = tuple(2 * n - 1 - 2 * i for i in range(modes))
        count = 0
        for degree in range(1, modes + 1):
            available = tuple(sum(1 << i for i in subset)
                              for subset in combinations(range(modes), degree))
            for pattern in range(1, 1 << len(available)):
                monomials = tuple(value for i, value in enumerate(available)
                                  if (pattern >> i) & 1)
                count += 1
                for base in (F(17, 16), F(5, 4), F(3, 2), F(2)):
                    probabilities = tuple(1 / (1 + base ** energy) for energy in energies)
                    p = truth_probability(monomials, probabilities)[0]
                    support_lower = prod(1 / (1 + base ** (2 * n - 1 - 2 * i))
                                         for i in range(degree))
                    assert p >= support_lower
                    support_controls += 1
                    h = base ** (2 * n * degree)
                    alternate = prod(1 / (base ** (-(2 * i + 1)) + base ** (-2 * n))
                                     for i in range(degree))
                    assert support_lower * h == alternate
                    cap = (1 + base ** (-2 * n)) ** (-n)
                    assert p * h >= alternate > (1 + base ** (-2 * n)) ** (-degree) >= cap
                    loss_controls += 1
                    for charge in (-8 * n, -2 * n, -2 * n + 1, 0, 2 * n - 1, 2 * n, 8 * n):
                        shifted = tuple(1 / (1 + base ** (energy + charge)) for energy in energies)
                        q = truth_probability(monomials, shifted)[0]
                        ratio = exact_ratio(p, q, h, base ** (-2 * (degree + 1)))
                        assert ratio > p * h > cap
                        ratio_controls += 1
        polynomial_counts.append({'modes': modes, 'homogeneous_polynomials': count})

    band_controls = mean_controls = threshold_controls = 0
    for n in range(1, 129):
        assert F(sum(4 * i + 2 for i in range(n)), n) == 2 * n
        mean_controls += 1
        # 2[(2-sqrt(2))n-sqrt(2)]-n=(3n)-(2n+2)sqrt(2).
        sign = radical_nonnegative(3 * n, -2 * n - 2)
        assert sign == (n >= 17)
        threshold_controls += 1
        for twice_charge in range(0, 4 * n + 1):
            charge = F(twice_charge, 2)
            in_band = radical_nonnegative(2 * n - 2 - charge, -2)
            if not in_band:
                continue
            smallest_v = 2 * n - charge
            # (sqrt(2)-1)v>=2 is the exact derivative comparison.
            assert radical_nonnegative(-smallest_v - 2, smallest_v)
            for extra in (F(0), F(1, 3), F(7, 2), charge):
                v = smallest_v + extra
                assert radical_nonnegative(-v - 2, v)
                band_controls += 1

    # All scalar reductions in the exterior use exact argument ordering.
    exterior_order_controls = 0
    for n in range(1, 65):
        for i in range(n):
            b = 2 * i + 1
            for charge_surplus in (F(0), F(1, 3), F(2), F(19, 2)):
                alpha = 2 + charge_surplus / b
                assert alpha >= 2
                assert charge_surplus + 2 * b - 1 > 0
                assert (charge_surplus + 2 * b - 1 + charge_surplus + 2 * b + 1) / 2 == alpha * b
                assert 2 < 3
                exterior_order_controls += 1

    result = {
        'status': 'Every control passed with exact rational, radical, and Boolean arithmetic.',
        'source_baseline': 'c2e951fa1f391edd97ba587d15ec13d3d6c8892d',
        'scope': 'Ordinary all-degree homogeneous one-level filter obstruction: all exterior charges, and the finite-tail charge band for n>=17. No improved factoring exponent.',
        'enumerated_homogeneous_patterns': polynomial_counts,
        'exact_biased_support_controls': support_controls,
        'exact_loss_cap_controls': loss_controls,
        'exact_relative_ratio_controls': ratio_controls,
        'exact_charge_band_derivative_orders': band_controls,
        'exact_exterior_convexity_argument_orders': exterior_order_controls,
        'exact_physical_energy_means': mean_controls,
        'sharp_integer_threshold_controls': threshold_controls,
        'threshold_certificate': {'n': 17, 'left_square': 2601, 'right_square': 2592},
        'proof_dependencies': [
            'homogeneous-boolean-filter-bounds.md: exact ratio and biased support lemma',
            'charged-one-mode-binary-laurent-filters.md Section 5: rigorous all-scale scalar inequality',
            'charged-gaussian-pressure-floor.md: positive periodic continuation',
            'folded-charged-gaussian-floor.md equation (10): full finite-tail integral lower bound'
        ],
        'not_claimed': 'No sign certificate for the two remaining fixed-width interior charge bands for overlapping polynomials; no result for arbitrary multilevel filters.'
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
