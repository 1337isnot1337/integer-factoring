#!/usr/bin/env python3
"""Exact Boolean support/transform controls and independent full-profile costs."""
import argparse
from itertools import combinations
import json
from math import exp, log, prod
from pathlib import Path
import random

from homogeneous_boolean_filter_bounds import (
    bounded_charge_small_scale_lower_bound, distinct_index_log_lower_bound,
    filter_data, relative_formula, scale_threshold, truth_mask,
)
from joint_identity_filter_pressure import pressure, relative_log_cost


def weighted_support_check(modes, mask, degree):
    denominators = tuple(range(2, modes + 2))
    common = prod(denominators)
    mass = sum(prod(1 if (state >> i) & 1 else denominator - 1
                    for i, denominator in enumerate(denominators))
               for state in range(1 << modes) if (mask >> state) & 1)
    assert mass * prod(sorted(denominators, reverse=True)[:degree]) >= common


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=root / 'evidence' / 'homogeneous-boolean-filter-bounds.json')
    args = parser.parse_args()
    support_controls = transform_controls = endpoint_controls = 0
    for modes in range(1, 5):
        monomial_masks = tuple(truth_mask(modes, (value,)) for value in range(1 << modes))
        for coefficients in range(1, 1 << (1 << modes)):
            mask = 0
            degree = 0
            remaining = coefficients
            while remaining:
                bit = remaining & -remaining
                monomial = bit.bit_length() - 1
                mask ^= monomial_masks[monomial]
                degree = max(degree, monomial.bit_count())
                remaining ^= bit
            weighted_support_check(modes, mask, degree)
            support_controls += 1
    for modes in range(1, 6):
        for degree in range(1, modes + 1):
            available = tuple(sum(1 << i for i in subset)
                              for subset in combinations(range(modes), degree))
            for coefficients in range(1, 1 << len(available)):
                monomials = tuple(value for i, value in enumerate(available) if (coefficients >> i) & 1)
                mask = truth_mask(modes, monomials)
                weighted_support_check(modes, mask, degree)
                support_controls += 1
                for m in (1, 3):
                    for increment in (1, 2):
                        final = m * degree + increment
                        values = [1 ^ (((1 << (m * degree)) ^ (1 << final))
                                       if (mask >> state) & 1 else 0)
                                  for state in range(1 << modes)]
                        assert all(value.bit_count() % 2 for value in values)
                        endpoint_controls += len(values)
                        for axis in range(modes):
                            bit = 1 << axis
                            for state in range(1 << modes):
                                if state & bit:
                                    source = values[state ^ bit]
                                    values[state] ^= source ^ (source << m)
                        expected = [1 << (m * state.bit_count() +
                                          (increment if (mask >> state) & 1 else 0))
                                    for state in range(1 << modes)]
                        assert values == expected
                        transform_controls += 1
    rng = random.Random(0xB001EA)
    discrepancies = []
    for trial in range(80):
        n, modes = rng.choice((8, 16, 32)), rng.randrange(1, 6)
        ramification = rng.choice((1, 3))
        degree = rng.randrange(1, modes + 1)
        available = [sum(1 << i for i in subset) for subset in combinations(range(modes), degree)]
        monomials = tuple(value for value in available if rng.randrange(2)) or (available[0],)
        final = n * ramification * degree + rng.choice((1, n, 3 * n))
        data = filter_data(modes, n, ramification, monomials, final)
        selected = tuple(sorted(rng.sample(range(n), modes)))
        charge = rng.choice((-3.0, 0.0, 1.0, 8.0, -2.0 * n))
        scale = rng.choice((.25, 1.0, 3.0, 6.0)) / n
        direct = relative_log_cost(2, n, ramification, charge, scale, selected,
                                   data['beta_orders'], data['degrees'])
        formula = relative_formula(n, ramification, selected, charge, scale,
                                   degree, final, data['truth_mask'])
        discrepancies.append(abs(direct - formula))
        assert abs(direct - formula) < 2e-10
        assert direct + 2e-10 >= distinct_index_log_lower_bound(n, degree, scale)
        if charge == 0 or scale >= scale_threshold(n):
            assert direct > 0
        if scale <= scale_threshold(n):
            full = pressure(2, n, ramification, charge, scale, selected,
                            data['beta_orders'], data['degrees'])
            assert full + 2e-10 >= bounded_charge_small_scale_lower_bound(n, degree, abs(charge))
    data = filter_data(1, 32, 1, (1,), 266)
    example_relative = relative_log_cost(2, 32, 1, -1000, .01, (0,),
                                        data['beta_orders'], data['degrees'])
    example_pressure = pressure(2, 32, 1, -1000, .01, (0,),
                                data['beta_orders'], data['degrees'])
    assert example_relative < 0 < example_pressure
    result = {
        'status': 'All exact Boolean support and polynomial transform controls passed.',
        'exact_biased_support_inequalities': support_controls,
        'complete_polynomial_transform_identities': transform_controls,
        'endpoint_equalities': endpoint_controls,
        'independent_full_profile_cost_comparisons': len(discrepancies),
        'maximum_formula_discrepancy': max(discrepancies),
        'small_scale_relative_gain_example': {
            'n': 32, 'ramification': 1, 'selected': [0], 'charge': -1000,
            'scale': .01, 'Boolean_degree': 1, 'final_polynomial_degree': 266,
            'relative_log_cost': example_relative, 'critical_pressure': example_pressure,
            'interpretation': 'A genuine relative gain at a very costly charged Gaussian, not a threshold witness.',
        },
        'bounded_charge_examples': [
            {'n': n, 'Boolean_degree': degree, 'charge_bound': charge,
             'proved_expression_evaluated_numerically': bounded_charge_small_scale_lower_bound(n, degree, charge)}
            for n, degree, charge in ((64, 1, 1), (128, 1, 3), (1024, 2, 3), (8192, 3, 3))
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('bounded_charge_examples', 'small_scale_relative_gain_example')}))


if __name__ == '__main__':
    main()
