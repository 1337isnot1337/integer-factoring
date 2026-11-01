#!/usr/bin/env python3
"""Exact all-field flag budgets and finite rank-sum exclusion through 82.

The theorem is proved in research/flagged-cut-rank-budget.md. All certificate
decisions use integers and fractions. No numerical optimizer is used.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
from random import Random

from check_weighted_cut_rank_budget import (
    column_flattening, crossing, encoded_fraction, exact_above_sqrt2_power,
    exact_below_sqrt2_power, lower_interpolant, matrix_rank,
    product_operator, random_matrix, weights,
)
from check_weighted_cut_residual_budget import chord_coefficients


@lru_cache(None)
def flag_complement(d: int, m: int) -> int:
    """Lower bound for r(E-B) when B is a sum of at most m products."""
    assert d >= 2 and m >= 0
    if m == 0:
        return d*d
    if m <= 2:
        return d*(d-1)//2
    q = m-1
    t = d//(2*q)
    return max(t*(d-q*t), (t+1)*(d-q*(t+1)))


@lru_cache(None)
def flag_budget(d: int, s: int, minimum_products: int = 0) -> int:
    assert 0 <= s <= d and minimum_products >= 0
    # Beyond d+1 the complement bound is zero and s^2*m is nondecreasing.
    end = max(minimum_products, d+1)
    return min(flag_complement(d, m)+s*s*m
               for m in range(minimum_products, end+1))


def harmonic_bound(d: int) -> Fraction:
    return (Fraction(2*d*d, 2*d-1)
            + 2*sum((w*flag_budget(d, s)
                     for s, w in enumerate(weights(d), 1)), Fraction()))


def class_parameters(d: int, k: int, b: int) -> dict[str, int]:
    n = d*d
    assert 2 <= k <= d and 1 <= b <= (n-1)//((k-1)**2)
    a0, a1 = (k-1)**2*b+1, min(k*k*b, n)
    f = flag_complement(d, b)
    rr, cc = max(n-a1, f), n-b
    if b == 1:
        rr = max(rr, n-k*(k+1)//2)
        extra = (k*k-3*k+4)//2
    else:
        extra = max(0, a0+f-n)
    return {'a_lower': a0, 'a_upper': a1, 'row_rest': rr,
            'col_rest': cc, 'extra_rank_budget': extra}


def class_bound(d: int, k: int, b: int) -> Fraction:
    p = class_parameters(d, k, b)
    a0, rr, cc = p['a_lower'], p['row_rest'], p['col_rest']
    result = Fraction(2*d*d+p['extra_rank_budget'], 2*k-1)
    for s, weight in enumerate(weights(k), 1):
        mu, nu = chord_coefficients(k, s)
        left = max(Fraction(flag_budget(d, s, b)),
                   min(a0, s*s*b)+mu*rr+nu*cc)
        right = max(Fraction(flag_budget(d, s)), b+nu*rr+mu*cc)
        result += weight*(left+right)
    return result


def dimension_certificate(d: int) -> dict[str, object]:
    best, selected, count = Fraction(d*d), {'k': 1}, 1
    assert exact_above_sqrt2_power(best, d)
    local_rows = []
    for k in range(2, d+1):
        local, local_b = None, None
        for b in range(1, (d*d-1)//((k-1)**2)+1):
            value = class_bound(d, k, b)
            assert exact_above_sqrt2_power(value, d), (d, k, b)
            count += 1
            if local is None or value < local:
                local, local_b = value, b
            if value < best:
                best, selected = value, {'k': k, 'b': b,
                                        **class_parameters(d, k, b)}
        local_rows.append({'k': k, 'b': local_b,
                           'bound': encoded_fraction(local)})
    return {'d': d, 'bound': encoded_fraction(best), 'above_target': True,
            'scalar_classes': count, 'minimizing_class': selected,
            'minimum_by_k': local_rows}


def check_scalar_bounds() -> dict[str, int]:
    maxima, closures, strengthened = 0, 0, 0
    for d in range(2, 91):
        for m in range(3, 2*d+5):
            assert flag_complement(d, m) == max(
                t*(d-(m-1)*t) for t in range(d+1))
            maxima += 1
        for s in range(d+1):
            for b in (0, 1, 2, 3, d, d+1, d+3):
                expected = min(flag_complement(d, m)+s*s*m
                               for m in range(b, max(b, 2*d+4)+1))
                assert flag_budget(d, s, b) == expected
                closures += 1
            assert flag_budget(d, s) >= s*d
            if 1 <= s <= (d-1)//2:
                assert flag_budget(d, s) >= s*d+s*s
                strengthened += 1
    return {'quadratic_maxima': maxima, 'finite_minimum_closures': closures,
            'strictly_strengthened_small_s_budgets': strengthened}


def check_actual_profile(d: int, profile: list[tuple[int, int]]) -> int:
    profile = [(r, c) for r, c in profile if r or c]
    assert all(r > 0 and c > 0 for r, c in profile)
    checks = 0
    for mask in range(1 << len(profile)):
        m = sum(c for i, (_, c) in enumerate(profile) if mask >> i & 1)
        rest = sum(r for i, (r, _) in enumerate(profile) if not mask >> i & 1)
        assert rest >= flag_complement(d, m)
        m = sum(r for i, (r, _) in enumerate(profile) if mask >> i & 1)
        rest = sum(c for i, (_, c) in enumerate(profile) if not mask >> i & 1)
        assert rest >= flag_complement(d, m)
        checks += 2
    for s in range(d+1):
        assert sum(min(r, s*s*c) for r, c in profile) >= flag_budget(d, s)
        assert sum(min(s*s*r, c) for r, c in profile) >= flag_budget(d, s)
        checks += 2
    ratio, index, exchange = max(
        [(Fraction(r, c), i, False) for i, (r, c) in enumerate(profile)]
        + [(Fraction(c, r), i, True) for i, (r, c) in enumerate(profile)])
    if ratio == 1:
        assert sum(r for r, _ in profile) >= d*d
        return checks
    if exchange:
        profile = [(c, r) for r, c in profile]
    a, b = profile[index]
    k = isqrt((a-1)//b)+1
    for s in range(1, k):
        selected = [i for i, (r, c) in enumerate(profile) if r > s*s*c]
        assert index in selected
        assert sum(profile[i][1] for i in selected) >= b
        assert sum(min(r, s*s*c) for r, c in profile) >= flag_budget(d, s, b)
        checks += 1
    total = sum((lower_interpolant(k, r, c) for r, c in profile), Fraction())
    assert total >= class_bound(d, k, b)
    return checks+1


def check_finite_tensors() -> dict[str, int | list[int]]:
    rng = Random(8210302026)
    profiles, inequalities, complements = 0, 0, 0
    for prime in (2, 3, 5):
        for d in range(2, 6):
            n, swap = d*d, crossing(d)
            for count in (1, 2, 3, 4):
                for _ in range(5):
                    terms = [product_operator(random_matrix(d, prime, rng),
                                              random_matrix(d, prime, rng), prime)
                             for _ in range(count)]
                    combined = [[sum(t[i][j] for t in terms) % prime
                                 for j in range(n)] for i in range(n)]
                    rest = [[(swap[i][j]-combined[i][j]) % prime
                             for j in range(n)] for i in range(n)]
                    assert matrix_rank(rest, prime) >= flag_complement(d, count)
                    complements += 1
                    matrices = terms+[rest]
                    profile = [(matrix_rank(t, prime),
                                matrix_rank(column_flattening(t, d), prime))
                               for t in matrices]
                    inequalities += check_actual_profile(d, profile)
                    profiles += 1
    return {'primes': [2, 3, 5], 'dimensions': [2, 3, 4, 5],
            'actual_decompositions': profiles, 'complement_checks': complements,
            'subset_and_weighted_inequalities': inequalities,
            'random_seed': 8210302026}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    assert 577**2-2*408**2 == 1
    assert 1393**2-2*985**2 == -1
    scalar = check_scalar_bounds()
    finite = check_finite_tensors()
    simple = []
    for d in range(2, 78):
        value = harmonic_bound(d)
        assert exact_above_sqrt2_power(value, d)
        simple.append({'d': d, 'bound': encoded_fraction(value)})
    simple_control = harmonic_bound(78)
    assert exact_below_sqrt2_power(simple_control, 78)
    dimensions = [dimension_certificate(d) for d in range(78, 83)]
    control = class_bound(83, 45, 3)
    assert exact_below_sqrt2_power(control, 83)
    here = Path(__file__)
    names = [here.name, 'check_weighted_cut_rank_budget.py',
             'check_weighted_cut_residual_budget.py', 'check_two_product_complement.py']
    result = {'status': 'PASS', 'scope': 'All-field exact rank-sum exclusion through 82; no factoring improvement',
              'last_certified_dimension': 82, 'simple_bound_cutoff': 77,
              'arithmetic': 'Integers and rational numbers only',
              'scalar_checks': scalar, 'finite_field_checks': finite,
              'simple_dimensions': simple,
              'simple_dimension_78_control': encoded_fraction(simple_control),
              'dimensions': dimensions,
              'dimension_83_control': {'k': 45, 'b': 3, 'below_target': True,
                  'bound': encoded_fraction(control),
                  'meaning': 'This certificate is insufficient; no realizable seed is asserted'},
              'source_sha256': {name: hashlib.sha256(here.with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: val for key, val in result.items()
                     if key not in ('simple_dimensions', 'dimensions')}, indent=2))
    for row in dimensions:
        print(json.dumps({key: val for key, val in row.items() if key != 'minimum_by_k'}))


if __name__ == '__main__':
    main()
