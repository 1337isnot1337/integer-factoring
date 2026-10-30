#!/usr/bin/env python3
"""Exact all-field geometric rank-sum exclusion through 95 labels."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from check_geometric_subspace_compression import geometric_complement
from check_flagged_cut_rank_budget import flag_complement, flag_budget
from check_weighted_cut_rank_budget import (
    encoded_fraction, exact_above_sqrt2_power, exact_below_sqrt2_power, weights,
)
from check_weighted_cut_residual_budget import chord_coefficients


@lru_cache(None)
def exact_geometric_complement(d: int, m: int) -> int:
    assert d >= 2 and m >= 0
    if m == 0:
        return d*d
    if m <= 2:
        return d*(d-1)//2
    return geometric_complement(d, m)


@lru_cache(None)
def geometric_budget(d: int, s: int, minimum_products: int = 0) -> int:
    assert 0 <= s <= d and minimum_products >= 0
    # Geometric bounds can be positive beyond the old flag cutoff d+1.
    end = max(minimum_products, 2*d-1)
    return min(exact_geometric_complement(d, m)+s*s*m
               for m in range(minimum_products, end+1))


def parameters(d: int, k: int, b: int) -> dict[str, int]:
    n = d*d
    assert 2 <= k <= d and 1 <= b <= (n-1)//((k-1)**2)
    a0, a1 = (k-1)**2*b+1, min(k*k*b, n)
    residual = exact_geometric_complement(d, b)
    rr, cc = max(n-a1, residual), n-b
    if b == 1:
        rr = max(rr, n-k*(k+1)//2)
        extra = (k*k-3*k+4)//2
    else:
        extra = max(0, a0+residual-n)
    return {'a_lower': a0, 'a_upper': a1, 'row_rest': rr,
            'col_rest': cc, 'extra_rank_budget': extra}


def class_bound(d: int, k: int, b: int) -> Fraction:
    p = parameters(d, k, b)
    a, rr, cc = p['a_lower'], p['row_rest'], p['col_rest']
    result = Fraction(2*d*d+p['extra_rank_budget'], 2*k-1)
    for s, weight in enumerate(weights(k), 1):
        mu, nu = chord_coefficients(k, s)
        first = max(Fraction(geometric_budget(d, s, b)),
                    min(a, s*s*b)+mu*rr+nu*cc)
        second = max(Fraction(geometric_budget(d, s)), b+nu*rr+mu*cc)
        result += weight*(first+second)
    return result


def dimension_certificate(d: int) -> dict[str, object]:
    best, selected, classes = Fraction(d*d), {'k': 1}, 1
    assert exact_above_sqrt2_power(best, d)
    for k in range(2, d+1):
        for b in range(1, (d*d-1)//((k-1)**2)+1):
            lower = class_bound(d, k, b)
            assert exact_above_sqrt2_power(lower, d), (d, k, b)
            classes += 1
            if lower < best:
                best, selected = lower, {'k': k, 'b': b, **parameters(d, k, b)}
    assert classes == 1+sum((d*d-1)//((k-1)**2) for k in range(2, d+1))
    return {'d': d, 'bound': encoded_fraction(best), 'above_target': True,
            'scalar_classes': classes, 'minimizing_class': selected,
            'complement_table': [exact_geometric_complement(d, m) for m in range(2*d)],
            'unconditioned_budgets': [geometric_budget(d, s) for s in range(d+1)]}


def check_closure_and_monotonicity() -> dict[str, int]:
    complement_checks, budget_checks, closure_checks = 0, 0, 0
    for d in range(2, 97):
        values = [exact_geometric_complement(d, m) for m in range(2*d+4)]
        for m, value in enumerate(values):
            assert value >= flag_complement(d, m)
            assert value == 0 if m >= 2*d-1 else value >= 0
            if m:
                assert values[m-1] >= value
            complement_checks += 1
        for s in range(d+1):
            for b in (0, 1, 2, 3, d, d+1, 2*d-1, 2*d+1):
                result = geometric_budget(d, s, b)
                assert result >= flag_budget(d, s, b)
                budget_checks += 1
                expanded = min(exact_geometric_complement(d, m)+s*s*m
                               for m in range(b, max(b, 3*d)+1))
                assert result == expanded
                closure_checks += 1
    return {'complement_comparisons': complement_checks,
            'budget_comparisons': budget_checks,
            'expanded_minimum_closures': closure_checks}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    assert 577**2-2*408**2 == 1
    assert 1393**2-2*985**2 == -1
    checks = check_closure_and_monotonicity()
    dimensions = [dimension_certificate(d) for d in range(83, 96)]
    control = class_bound(96, 52, 3)
    assert exact_below_sqrt2_power(control, 96)
    here = Path(__file__)
    names = [here.name, 'check_geometric_subspace_compression.py',
             'check_flagged_cut_rank_budget.py', 'check_weighted_cut_rank_budget.py',
             'check_weighted_cut_residual_budget.py']
    result = {'status': 'PASS', 'scope': 'All-field exact finite rank-sum exclusion through 95',
              'previous_cutoff': 82, 'last_certified_dimension': 95,
              'improved_factoring_exponent': False,
              'arithmetic': 'Integers and fractions; no numerical optimizer',
              'finite_minimum_cap': 'max(minimum_products, 2*d-1)',
              'scalar_classes': sum(row['scalar_classes'] for row in dimensions),
              'structural_checks': checks, 'dimensions': dimensions,
              'dimension_96_control': {'k': 52, 'b': 3, 'bound': encoded_fraction(control),
                  'below_target': True, 'meaning': 'Certificate limit; no realizable seed is asserted'},
              'source_sha256': {name: hashlib.sha256(here.with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'dimensions'}, indent=2))
    for row in dimensions:
        print(json.dumps({key: value for key, value in row.items()
                          if key not in ('complement_table', 'unconditioned_budgets')}))


if __name__ == '__main__':
    main()
