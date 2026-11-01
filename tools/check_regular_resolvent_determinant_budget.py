#!/usr/bin/env python3
"""Exact weighted AM-GM certificates for regular resolvent profiles.

Floating point is used only with --generate to choose integer witnesses.
Every asserted inequality is then checked with integers or Fractions.
"""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, exp, gcd, log
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'evidence/regular-resolvent-determinant-budget.json'
LABEL_SIZES = tuple(range(2, 19)) + (21,)


def threshold(n):
    """Least positive K with K^408 > n^577; 577/408 > sqrt(2)."""
    target, lo, hi = n**577, 0, n*n
    while hi - lo > 1:
        middle = (lo + hi)//2
        if middle**408 > target:
            hi = middle
        else:
            lo = middle
    return hi


def row_budget(n, r):
    return max(0, (n + r*r - 1)//(r*r) - 2*r)


def weights(n, r, k, a, b, c):
    ell, degree = n-r-k, n-1
    assert ell > 0 and min(a, b) > 0 and (k == 0 or c > 0)
    row = [r*b, k*c, ell*a]
    total = sum(row)
    numerator = degree*r*(b-a) + (degree-1)*k*(c-a)
    assert numerator % n == 0
    v = numerator//n
    assert 0 < v < total
    # Cancellation of the two weighted exponent means is exact.
    row_mean = (Fraction(row[0]*degree + row[1]*(degree-1), total)
                - Fraction(row[2]*(r*degree + k*(degree-1)), ell*total))
    assert row_mean == Fraction(n*v, total)
    return row, total, v


def verify_amgm(item):
    n, r, k = (item[key] for key in ('n', 'r', 'k'))
    a, b, c = (item[key] for key in ('a', 'b', 'c'))
    row, total, v = weights(n, r, k, a, b, c)
    bound = item['integer_bound']
    numerator = total**(2*total) * (n-1)**v
    denominator = (b**row[0] * c**row[1] * a**row[2]
                   * (total-v)**(total-v) * v**v)
    assert numerator > denominator * bound**total
    return {'weight_denominator': total, 'column_second_weight': v,
            'numerator_bits': numerator.bit_length(),
            'denominator_bits': denominator.bit_length(),
            'numerator_sha256': sha256(numerator.to_bytes((numerator.bit_length()+7)//8, 'big')).hexdigest(),
            'denominator_sha256': sha256(denominator.to_bytes((denominator.bit_length()+7)//8, 'big')).hexdigest()}


def log_sum(terms):
    top = max(terms)
    return top + log(sum(exp(value-top) for value in terms))


def choose_witness(n, r, k, bound):
    ell, degree = n-r-k, n-1
    last = -Fraction(r*degree + k*(degree-1), ell)
    def objective(y):
        terms = [log(r)+degree*y, log(ell)+float(last)*y]
        if k:
            terms.append(log(k)+(degree-1)*y)
        return log_sum(terms) + log_sum([0., log(n-1)-n*y])
    left, right = 0., 4.
    for _ in range(100):
        x, y = (2*left+right)/3, (left+2*right)/3
        if objective(x) < objective(y):
            right = y
        else:
            left = x
    scale = (left+right)/2
    ratio_b = exp((degree-float(last))*scale)
    ratio_c = exp((degree-1-float(last))*scale)
    choices = []
    common, step = gcd(r, n), n//gcd(r, n)
    inverse = pow(r//common, -1, step) if step > 1 else 0
    for a in range(1, 41):
        c_values = [1] if k == 0 else range(max(1, round(a*ratio_c)-n), round(a*ratio_c)+n+1)
        for c in c_values:
            rhs = r*a - 2*k*(c-a)
            if rhs % common:
                continue
            residue = ((rhs//common)*inverse) % step if step > 1 else 0
            nearest = round((a*ratio_b-residue)/step)
            for offset in (-1, 0, 1):
                b = residue + (nearest+offset)*step
                if b <= 0:
                    continue
                try:
                    row, total, v = weights(n, r, k, a, b, c)
                except AssertionError:
                    continue
                logarithm = (2*total*log(total) + v*log(n-1)
                             - row[0]*log(b) - row[1]*log(c) - row[2]*log(a)
                             - (total-v)*log(total-v) - v*log(v))/total
                if logarithm > log(bound) + 1e-7:
                    choices.append((total, a, b, c))
    if not choices:
        raise ValueError(('No exact-certificate candidate', n, r, k, bound, exp(objective(scale))))
    _, a, b, c = min(choices)
    item = {'n': n, 'r': r, 'k': k, 'a': a, 'b': b, 'c': c,
            'integer_bound': bound, 'method': 'weighted-amgm'}
    item['certificate'] = verify_amgm(item)
    return item


def generate():
    records = []
    for q in LABEL_SIZES:
        n, r0 = q*q, 1
        while 2*r0**3 < n:
            r0 += 1
        bound = threshold(n)
        cases = []
        for r in range(1, r0+1):
            k = row_budget(n, r)
            raw = n*r+(n-1)*k
            if raw**408 > n**577:
                cases.append({'n': n, 'r': r, 'k': k, 'method': 'early-pivots',
                              'integer_bound': raw})
            else:
                cases.append(choose_witness(n, r, k, bound))
        records.append({'q': q, 'n': n, 'zero_only_rank_threshold': r0, 'cases': cases})
    # These two short types also give the convenient all-size rank-two corollary.
    short_types = [dict(n=9, r=2, k=0, a=6, b=51, c=1, integer_bound=27, method='weighted-amgm'),
                   dict(n=16, r=2, k=0, a=9, b=121, c=1, integer_bound=64, method='weighted-amgm')]
    for item in short_types:
        item['certificate'] = verify_amgm(item)
    return {'label_sizes': list(LABEL_SIZES), 'all_rank_certificates': records,
            'short_rank_two_certificates': short_types}


def relaxation_witness():
    n, scale = 361, Fraction(73, 72)
    profile = [(0, 6), (366, 325), (367, 30)]
    assert sum(count for _, count in profile) == n
    assert sum(order*count for order, count in profile) == n*(n-1)
    assert row_budget(n, 6) == 0
    cost = sum(count*scale**(n-1-order) for order, count in profile) * (1+(n-1)*scale**(-n))
    assert cost.numerator**169 < n**239 * cost.denominator**169
    assert 239**2 < 2*169**2
    return {'q': 19, 'profile_order_multiplicity': profile, 'scale': [73, 72],
            'comparison': 'B^169 < 361^239, and 239/169 < sqrt(2)',
            'cost_numerator': str(cost.numerator), 'cost_denominator': str(cost.denominator),
            'actual_crossing_constructed': False}


def mobius_controls():
    """Check the cleared polynomial identity in F_p[S]/(S^n)."""
    count = 0
    for prime in (3, 5, 7):
        for gamma in range(2, prime):
            for n in (4, 9, 16):
                degree = n-1
                transformed = [0] + [(gamma-1)*pow(gamma, -j, prime) % prime
                                     for j in range(1, n)]
                def multiply(a, b):
                    return [sum(a[j]*b[i-j] for j in range(i+1)) % prime for i in range(n)]
                one_minus = [1] + [(-x) % prime for x in transformed[1:]]
                right, power = [[0]*n for _ in range(n)], [1]+[0]*(n-1)
                for j in range(n):
                    right[degree-j] = multiply(power, one_minus)
                    power = multiply(power, transformed)
                assert power == [0]*n
                assert pow(transformed[1], n-1, prime) != 0
                left = [[0]*n for _ in range(n)]
                for j in range(n):
                    for h in range(j+1):
                        coefficient = (pow(gamma, -j, prime)*comb(j, h)
                                       * pow(gamma-1, j-h, prime)) % prime
                        x_degree = degree-j+h
                        left[x_degree][j] = (left[x_degree][j]+coefficient) % prime
                        if j+1 < n:
                            left[x_degree][j+1] = (left[x_degree][j+1]-coefficient) % prime
                assert left == right
                count += 1
    return count


def verify(report):
    assert 577**2 > 2*408**2
    assert report['label_sizes'] == list(LABEL_SIZES)
    total = amgm = early = 0
    for q, record in zip(LABEL_SIZES, report['all_rank_certificates'], strict=True):
        n, r0 = q*q, record['zero_only_rank_threshold']
        assert record['q'] == q and record['n'] == n
        assert 2*(r0-1)**3 < n <= 2*r0**3
        assert len(record['cases']) == r0
        for r, item in enumerate(record['cases'], 1):
            k = row_budget(n, r)
            assert (item['n'], item['r'], item['k']) == (n, r, k)
            assert item['integer_bound']**408 > n**577
            if item['method'] == 'early-pivots':
                assert item['integer_bound'] == n*r+(n-1)*k
                early += 1
            else:
                assert item['certificate'] == verify_amgm(item)
                amgm += 1
            total += 1
        assert record['cases'][-1]['k'] == 0
    for item in report['short_rank_two_certificates']:
        assert item['certificate'] == verify_amgm(item)
        assert item['r'] == 2 and item['integer_bound']**2 == item['n']**3
    return {'status': 'passed', 'label_sizes': len(LABEL_SIZES),
            'rank_cases': total, 'early_pivot_certificates': early,
            'weighted_amgm_certificates': amgm, 'short_rank_two_certificates': 2,
            'mobius_polynomial_identities': mobius_controls(),
            'q19_relaxation_witness': relaxation_witness()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = generate() if args.generate else json.loads(args.evidence.read_text())
    report['verification'] = verify(report)
    report['source_sha256'] = {name: sha256((ROOT/name).read_bytes()).hexdigest()
                               for name in ('tools/check_regular_resolvent_determinant_budget.py',
                                            'research/partial-transpose-power-rank.md')}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report['verification'].items() if k != 'q19_relaxation_witness'}))


if __name__ == '__main__':
    main()
