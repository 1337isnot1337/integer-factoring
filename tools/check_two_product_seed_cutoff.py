#!/usr/bin/env python3
"""Exact 53/54-label rank-sum exclusion using the two-product theorem."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from check_weighted_cut_rank_budget import (
    encoded_fraction, exact_above_sqrt2_power, exact_below_sqrt2_power, weights,
)
from check_weighted_cut_residual_budget import (
    case_bound as previous_bound, case_parameters as previous_parameters,
    chord_coefficients,
)


def parameters(d, k, b):
    data = previous_parameters(d, k, b)
    if b == 2:
        complement = d * (d - 1) // 2
        data['row_rest'] = max(data['row_rest'], complement)
        data['extra_rank_budget'] = max(
            data['extra_rank_budget'], data['a_lower'] + complement - d*d)
    return data


def bound(d, k, b):
    if b != 2:
        return previous_bound(d, k, b)
    data = parameters(d, k, b)
    a, rr, cc = data['a_lower'], data['row_rest'], data['col_rest']
    result = Fraction(2*d*d + data['extra_rank_budget'], 2*k-1)
    for s, weight in enumerate(weights(k), 1):
        mu, nu = chord_coefficients(k, s)
        result += weight * (
            max(Fraction(s*d), min(a, s*s*b) + mu*rr + nu*cc)
            + max(Fraction(s*d), b + nu*rr + mu*cc))
    return result


def certify(d):
    best, selected, cases, improved = Fraction(d*d), {'k': 1}, 1, 0
    assert exact_above_sqrt2_power(best, d)
    by_k = []
    for k in range(2, d+1):
        local, selected_b = None, None
        upper = (d*d-1)//((k-1)**2)
        for b in range(1, upper+1):
            value = bound(d, k, b)
            assert exact_above_sqrt2_power(value, d), (d, k, b)
            if b == 2:
                assert value >= previous_bound(d, k, b)
                improved += value > previous_bound(d, k, b)
            cases += 1
            if local is None or value < local:
                local, selected_b = value, b
        by_k.append({'k': k, 'b': selected_b, 'bound': encoded_fraction(local),
                     'classes': upper})
        if local < best:
            best, selected = local, {'k': k, 'b': selected_b,
                                    **parameters(d, k, selected_b)}
    return {'d': d, 'lower_bound': encoded_fraction(best),
            'strict_above_d_to_sqrt2': True, 'scalar_classes': cases,
            'strictly_improved_b2_classes': improved,
            'minimizing_class': selected, 'minimum_by_k': by_k}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    assert 577**2 - 2*408**2 == 1
    assert 1393**2 - 2*985**2 == -1
    dimensions = [certify(d) for d in (53, 54)]
    control = bound(55, 51, 1)
    assert exact_below_sqrt2_power(control, 55)
    names = ['check_two_product_seed_cutoff.py',
             'check_weighted_cut_rank_budget.py',
             'check_weighted_cut_residual_budget.py',
             'check_two_product_complement.py']
    result = {'status': 'PASS', 'scope': 'All-field rank-sum exclusion through 54; no new factoring exponent',
              'previous_cutoff': 52, 'last_certified_dimension': 54,
              'arithmetic': 'Integers and rational numbers only',
              'dimensions': dimensions,
              'dimension_55_control': {'k': 51, 'b': 1, 'bound': encoded_fraction(control),
                  'below_target': True, 'meaning': 'This certificate alone fails; no realizable seed is asserted'},
              'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'dimensions'}, indent=2))
    for row in dimensions:
        print(json.dumps({k: v for k, v in row.items() if k != 'minimum_by_k'}))


if __name__ == '__main__':
    main()
