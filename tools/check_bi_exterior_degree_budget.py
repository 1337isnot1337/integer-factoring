#!/usr/bin/env python3
"""Exact degree and quadratic-rank certificates for the bi-exterior model."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from check_diagonal_word_rank_sum import square_zero_rank, complement_coefficients
from check_weighted_cut_rank_budget import (
    encoded_fraction, exact_above_sqrt2_power, exact_below_sqrt2_power,
)


def support_degrees(variables: int, polynomial: int) -> list[int]:
    return [a.bit_count() for a in range(1 << variables) if polynomial >> a & 1]


def quadratic_ranks(variables: int, half_rank: int) -> tuple[int, int]:
    assert variables >= 2 and 1 <= half_rank <= variables//2
    return (1 << (variables-1))-(1 << (variables-half_rank-1)), 2*half_rank+2


def scalar_catalog(width: int):
    """A safe superset of all degree-endpoint and quadratic-rank classes."""
    n, dimension = 2*width, 1 << (2*width)
    for p in range(n+1):
        for q in range(p, n+1):
            aa = range(width+1) if p <= 2 <= q else (0,)
            bb = range(width+1) if p <= n-2 <= q else (0,)
            for a in aa:
                if (p == 2 or q == 2) and a == 0:
                    continue
                for b in bb:
                    if (p == n-2 or q == n-2) and b == 0:
                        continue
                    r, c = 1 << (n-p), 1 << q
                    if p == 2:
                        r = max(r, quadratic_ranks(n, a)[0])
                    if q == n-2:
                        c = max(c, quadratic_ranks(n, b)[0])
                    if q == 2:
                        c = max(c, quadratic_ranks(n, a)[1])
                    if p == n-2:
                        r = max(r, quadratic_ranks(n, b)[1])
                    assert 1 <= r <= dimension and 1 <= c <= dimension
                    coverage = [int(p <= 2*j <= q) for j in range(width+1)]
                    if a == 0:
                        coverage[1] = 0
                    if b == 0:
                        coverage[-2] = 0
                    yield p, q, a, b, r, c, coverage


def rational_certificate(width: int) -> dict[str, object]:
    assert width in (8, 9)
    d = 1 << width
    grade_weights = [int(j not in (1, width-1)) for j in range(width+1)]
    classes = 0
    minimum_slack = None
    for p, q, a, b, r, c, coverage in scalar_catalog(width):
        weighted_degrees = sum(x*y for x, y in zip(grade_weights, coverage))
        scaled_dual = 5*weighted_degrees+a+b
        slack = 25*r*c-d*d*scaled_dual*scaled_dual
        assert slack >= 0, (width, p, q, a, b, r, c)
        classes += 1
        minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
    lower = d*Fraction(sum(grade_weights)*5+2*width, 5)
    return {'width': width, 'labels': d, 'scalar_classes': classes,
            'grade_weights': grade_weights, 'quadratic_half_rank_weight': '1/5',
            'bound': encoded_fraction(lower), 'minimum_squared_slack': minimum_slack}


def finite_controls() -> dict[str, int]:
    polynomials = quadratics = degree_counts = 0
    for variables in (2, 4):
        d = 1 << (variables//2)
        for f in range(1, 1 << (1 << variables)):
            degrees = support_degrees(variables, f)
            p, q = min(degrees), max(degrees)
            r = square_zero_rank(variables, f)
            c = square_zero_rank(variables, complement_coefficients(variables, f))
            assert r >= 1 << (variables-p)
            assert c >= 1 << q
            assert r*c >= d*d*(1 << (q-p))
            even_count = len({x for x in degrees if x % 2 == 0})
            assert r*c >= d*d*even_count*even_count
            polynomials += 1
    for variables in range(2, 11, 2):
        for a in range(1, variables//2+1):
            f = sum(1 << (3 << (2*j)) for j in range(a))
            assert square_zero_rank(variables, f) == quadratic_ranks(variables, a)[0]
            dual = complement_coefficients(variables, f)
            assert square_zero_rank(variables, dual) == quadratic_ranks(variables, a)[1]
            quadratics += 1
    for k in range(1, 257):
        assert 1 << (k-1) >= k
        degree_counts += 1
    return {'all_nonzero_small_polynomials': polynomials,
            'quadratic_normal_forms': quadratics,
            'degree_count_inequalities': degree_counts}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    controls = finite_controls()
    elementary = []
    for width in range(1, 8):
        d = 1 << width
        lower = Fraction(d*(width+1))
        assert exact_above_sqrt2_power(lower, d)
        elementary.append({'width': width, 'labels': d, 'bound': int(lower)})
    final = rational_certificate(8)
    lower = Fraction(int(final['bound']['numerator']), int(final['bound']['denominator']))
    assert exact_above_sqrt2_power(lower, 256)
    control = rational_certificate(9)
    below = Fraction(int(control['bound']['numerator']), int(control['bound']['denominator']))
    assert exact_below_sqrt2_power(below, 512)
    here = Path(__file__)
    result = {'status': 'PASS',
        'scope': 'Full characteristic-two bi-exterior space; exact and supplied formal limits',
        'universal_bound': 'S >= d*(lambda+1), d=2^lambda',
        'last_certified_binary_width': 8, 'last_certified_labels': 256,
        'improved_factoring_exponent': False, 'controls': controls,
        'degree_only_dimensions': elementary, 'dimension_256_certificate': final,
        'dimension_512_control': {**control,
            'meaning': 'This scalar certificate is below target; no realizable seed is asserted'},
        'source_sha256': {name: hashlib.sha256(here.with_name(name).read_bytes()).hexdigest()
                         for name in (here.name, 'check_diagonal_word_rank_sum.py',
                                      'check_weighted_cut_rank_budget.py')}}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
