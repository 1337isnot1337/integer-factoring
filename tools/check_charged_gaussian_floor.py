#!/usr/bin/env python3
"""Replay a rational rectangle cover and exact charged Gaussian identities."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

from charged_gaussian_floor import Dyadic, split_rectangle


def interval_controls(arithmetic):
    count = 0
    for first in (Fraction(-7, 3), Fraction(0), Fraction(1, 9), Fraction(5, 2)):
        for second in (Fraction(-3, 5), Fraction(2, 7), Fraction(4)):
            a = arithmetic.rational(first.numerator, first.denominator)
            b = arithmetic.rational(second.numerator, second.denominator)
            for actual, exact in ((arithmetic.add(a, b), first + second),
                                  (arithmetic.subtract(a, b), first - second),
                                  (arithmetic.multiply(a, b), first * second)):
                assert Fraction(actual[0], arithmetic.scale) <= exact <= Fraction(actual[1], arithmetic.scale)
                count += 1
    assert sum(Fraction(3 ** k, factorial(k)) for k in range(9)) > 20
    assert sum(Fraction(12, 5) ** k / factorial(k) for k in range(12)) > 11
    for value in (0, 1, 2, 7):
        positive = arithmetic.softplus_point(value * arithmetic.scale)
        negative = arithmetic.softplus_point(-value * arithmetic.scale)
        difference = arithmetic.subtract(negative, positive)
        assert difference[0] <= value * arithmetic.scale <= difference[1]
        count += 1
    return count


def charge_identity_controls():
    count = 0
    for n in range(1, 65):
        for k in range(-n, n + 1):
            charge = 2 * k
            physical = Counter(4 * i - 2 * (n - 1) + sign * charge
                               for i in range(n) for sign in (-1, 1))
            comparison = Counter(4 * i - 2 * (length - 1)
                                 for length in (n + k, n - k) for i in range(length))
            assert physical == comparison
            zero = [2 * i + 1 + sign * charge for i in range(n) for sign in (-1, 1)]
            infinity = [2 * i + 1 - 2 * n for i in range(2 * n)]
            assert -sum(min(value, 0) for value in zero + infinity) == n * n + k * k
            thermal = Counter(abs(value) for value in zero + infinity)
            expected = Counter(2 * i + 1 for length in (n + k, n - k, n, n) for i in range(length))
            assert thermal == expected
            count += 1
    return count


def parity_identity_controls():
    """Check the complete finite physical reductions for symbolic real J.

    A term (a,b,c) denotes c*sqrt(2)*f((a+b*J)*t/sqrt(2)).
    Normalizing its argument by f(-x)=f(x)+x leaves an exact
    polynomial ground term and a Counter of formal thermal terms.
    No choice of sign or numerical value of J is needed.
    """
    def normalize(terms, ground):
        constant, linear = ground
        thermal = Counter()
        for a, b, coefficient in terms:
            if a < 0 or (a == 0 and b < 0):
                constant -= coefficient * a
                linear -= coefficient * b
                a, b = -a, -b
            thermal[a, b] += coefficient
        return constant, linear, thermal

    for n in range(1, 129):
        raw = [(4 * i - 2 * (n - 1), sign, 1)
               for i in range(n) for sign in (-1, 1)]
        if n % 2 == 0:
            reduced = [(4 * i + 2, sign, 2)
                       for i in range(n // 2) for sign in (-1, 1)]
            ground = n * n, 0
        else:
            reduced = [(0, -1, 2)] + [
                (4 * i, sign, 2)
                for i in range(1, (n + 1) // 2) for sign in (-1, 1)]
            ground = n * n - 1, -1
        assert normalize(raw, (0, 0)) == normalize(reduced, ground)
    return 128


def dyadic_rectangle(arithmetic, rectangle):
    values = [arithmetic.rational(value.numerator, value.denominator) for value in rectangle]
    return (values[0][0], values[1][1]), (values[2][0], values[3][1])


def generate(arithmetic, terms):
    start = (Fraction(1), Fraction(6), Fraction(0), Fraction(2))
    pending = [('', start)]
    leaves = []
    evaluated = 0
    while pending:
        path, rectangle = pending.pop()
        t, J = dyadic_rectangle(arithmetic, rectangle)
        lower = arithmetic.continuation(t, J, terms)
        evaluated += 1
        if lower > 0:
            leaves.append([path, str(lower)])
        else:
            assert len(path) < 30, ('Unresolved rectangle', path, rectangle, lower)
            children = split_rectangle(rectangle)
            pending.extend((path + symbol, child) for symbol, child in reversed(children))
        if evaluated % 1000 == 0:
            print('RECTANGLES', evaluated, len(leaves), len(pending), flush=True)
        assert evaluated < 200000
    return leaves, evaluated


def verify_cover(arithmetic, terms, leaves):
    lookup = {path: int(lower) for path, lower in leaves}
    assert len(lookup) == len(leaves)
    prefixes = {path[:i] for path in lookup for i in range(len(path))}
    seen = set()
    pending = [('', (Fraction(1), Fraction(6), Fraction(0), Fraction(2)))]
    while pending:
        path, rectangle = pending.pop()
        if path in lookup:
            assert path not in prefixes
            t, J = dyadic_rectangle(arithmetic, rectangle)
            lower = arithmetic.continuation(t, J, terms)
            assert lower == lookup[path] and lower > 0
            seen.add(path)
        else:
            assert path in prefixes
            for symbol, child in split_rectangle(rectangle):
                pending.append((path + symbol, child))
    assert seen == set(lookup)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--certificate', type=Path,
                        default=root / 'evidence' / 'charged-gaussian-floor.json')
    args = parser.parse_args()
    arithmetic = Dyadic(80)
    exact_controls = interval_controls(arithmetic)
    identities = charge_identity_controls()
    parity_identities = parity_identity_controls()
    natural = arithmetic.finite_pressure(32, arithmetic.rational(3), arithmetic.zero)
    exterior = arithmetic.finite_pressure(32, arithmetic.rational(3), arithmetic.rational(323, 5))
    difference = arithmetic.subtract(exterior, natural)
    assert exterior[0] > 0 and natural[0] > 0 and difference[1] < 0
    minimization_counterexample = {
        'n': 32, 'scale': '3', 'charge': '323/5',
        'natural_pressure_numerators': list(map(str, natural)),
        'charged_pressure_numerators': list(map(str, exterior)),
        'charged_minus_natural_numerators': list(map(str, difference)),
        'interpretation': 'Both pressures are positive, but the exterior charge is strictly better than J=0.',
    }
    if args.generate:
        leaves, evaluated = generate(arithmetic, 20)
        record = {
            'status': 'Every rectangle has a strictly positive outward-rational continuation bound.',
            'dyadic_precision': 80, 'exponential_Taylor_degree': 24,
            'logarithm_atanh_terms': 32, 'paired_sum_terms': 20,
            'domain': {'scale': ['1', '6'], 'charge': ['0', '2']},
            'exact_interval_controls': exact_controls,
            'exact_even_charge_multiset_identities': identities,
            'exact_symbolic_charge_parity_identities': parity_identities,
            'rectangles_evaluated_during_generation': evaluated,
            'positive_cover_rectangles': len(leaves),
            'minimum_lower_bound_numerator': str(min(int(lower) for path, lower in leaves)),
            'common_denominator': str(arithmetic.scale),
            'counterexample_to_global_J0_minimization': minimization_counterexample,
            'cover_sha256': sha256(json.dumps(leaves, separators=(',', ':')).encode()).hexdigest(),
            'cover': leaves,
        }
        args.certificate.write_text(json.dumps(record, separators=(',', ':')) + '\n')
    else:
        record = json.loads(args.certificate.read_text())
        assert record['dyadic_precision'] == arithmetic.precision
        assert record['exponential_Taylor_degree'] == 24 and record['logarithm_atanh_terms'] == 32
        assert record['exact_interval_controls'] == exact_controls
        assert record['exact_even_charge_multiset_identities'] == identities
        assert record['exact_symbolic_charge_parity_identities'] == parity_identities
        assert record['domain'] == {'scale': ['1', '6'], 'charge': ['0', '2']}
        assert sha256(json.dumps(record['cover'], separators=(',', ':')).encode()).hexdigest() == record['cover_sha256']
        assert len(record['cover']) == record['positive_cover_rectangles']
        assert record['common_denominator'] == str(arithmetic.scale)
        assert record['minimum_lower_bound_numerator'] == str(min(int(lower) for path, lower in record['cover']))
        assert record['counterexample_to_global_J0_minimization'] == minimization_counterexample
        verify_cover(arithmetic, record['paired_sum_terms'], record['cover'])
    print(json.dumps({key: value for key, value in record.items() if key != 'cover'}))


if __name__ == '__main__':
    main()
