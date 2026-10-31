#!/usr/bin/env python3
"""Integer outward certificates for the three residual-charge Gaussian bases."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

from charged_gaussian_floor import Dyadic, split_rectangle


def absolute(interval):
    low, high = interval
    if low >= 0:
        return interval
    if high <= 0:
        return -high, -low
    return 0, max(-low, high)


def positive(interval):
    return max(interval[0], 0), max(interval[1], 0)


def cancelled_pressure(arithmetic, n, scale, charge, side):
    """Exact complete thermal expression with a piecewise linear ground."""
    a = arithmetic
    value = a.zero
    for index in range(n):
        odd = a.rational(2 * index + 1)
        for energy in (a.add(odd, charge), a.subtract(odd, charge)):
            value = a.add(value, a.softplus(a.multiply(absolute(energy), scale)))
        value = a.add(value, a.multiply(a.rational(2),
                                       a.softplus(a.multiply(odd, scale))))
    if n == 2 and side == 0:
        ground = positive(a.subtract(charge, a.one))
        physical = (a.subtract(a.rational(2), charge), a.add(a.rational(2), charge))
    elif n == 2:
        ground = positive(a.subtract(a.rational(3), charge))
        physical = (a.subtract(charge, a.rational(2)), a.add(a.rational(2), charge))
    else:
        ground = positive(a.subtract(a.one, charge))
        physical = (charge,) if n == 1 else (
            charge, a.subtract(a.rational(4), charge), a.add(a.rational(4), charge))
    value = a.add(value, a.multiply(ground, scale))
    negative = a.zero
    for energy in physical:
        assert energy[0] >= 0
        negative = a.add(negative, a.softplus(a.divide(a.multiply(energy, scale), a.root_two)))
    return a.subtract(value, a.multiply(a.multiply(a.rational(2), a.root_two), negative))


def cover_base(arithmetic, n, left_j, right_j, side):
    pending = [('', (Fraction(1, 32), Fraction(9), Fraction(left_j), Fraction(right_j)))]
    leaves = []
    evaluated = 0
    while pending:
        path, rectangle = pending.pop()
        endpoints = [arithmetic.rational(x.numerator, x.denominator) for x in rectangle]
        scale = endpoints[0][0], endpoints[1][1]
        charge = endpoints[2][0], endpoints[3][1]
        lower = cancelled_pressure(arithmetic, n, scale, charge, side)[0]
        evaluated += 1
        if lower > 0:
            leaves.append([path, str(lower)])
        else:
            assert len(path) < 28, (n, path, rectangle, lower)
            pending.extend((path + symbol, child)
                           for symbol, child in reversed(split_rectangle(rectangle)))
        assert evaluated < 100000
    return {'n': n, 'charge_interval': [left_j, right_j], 'physical_branch': side,
            'evaluated_rectangles': evaluated, 'positive_cover_rectangles': len(leaves),
            'minimum_lower_numerator': str(min(int(value) for path, value in leaves)),
            'cover': leaves}


def exact_profile_controls():
    from collections import Counter
    count = 0
    for n, maximum in ((1, 2), (2, 4), (3, 2)):
        for delta in (Fraction(k, 4) for k in range(4 * maximum + 1)):
            zero = [2 * i + 1 + sign * delta for i in range(n) for sign in (-1, 1)]
            infinity = [2 * i + 1 - 2 * n for i in range(2 * n)]
            physical = [4 * i - 2 * (n - 1) + sign * delta
                        for i in range(n) for sign in (-1, 1)]
            ground = -sum(min(x, 0) for x in zero + infinity) + sum(min(x, 0) for x in physical)
            if n == 2 and delta <= 2:
                expected_ground = max(delta - 1, 0)
                expected_physical = (2 - delta, 2 + delta)
            elif n == 2:
                expected_ground = max(3 - delta, 0)
                expected_physical = (delta - 2, 2 + delta)
            else:
                expected_ground = max(1 - delta, 0)
                expected_physical = (delta,) if n == 1 else (delta, 4 - delta, 4 + delta)
            assert ground == expected_ground
            assert Counter(abs(x) for x in physical) == Counter(x for x in expected_physical for _ in range(2))
            assert sum(max(x, 0) for x in zero + infinity) - sum(min(x, 0) for x in physical) <= 4 * n * n + 8
            count += 1
    for n in (1, 2, 3):
        assert Fraction(16 * n, 21) - Fraction(4 * n * n + 8, 32) > 0
    return count


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=root / 'evidence' / 'small-charged-gaussian-bases.json')
    args = parser.parse_args()
    arithmetic = Dyadic(80)
    bases = [cover_base(arithmetic, n, left, right, side)
             for n, left, right, side in ((1, 0, 2, 0), (2, 0, 2, 0),
                                          (2, 2, 4, 1), (3, 0, 2, 0))]
    record = {'status': 'All complete small-base pressure rectangles have positive integer outward bounds.',
              'dyadic_precision': 80, 'common_denominator': str(arithmetic.scale),
              'scale_interval': ['1/32', '9'],
              'exact_complete_profile_controls': exact_profile_controls(),
              'positive_cover_rectangles': sum(x['positive_cover_rectangles'] for x in bases),
              'bases': bases}
    args.output.write_text(json.dumps(record, separators=(',', ':')) + '\n')
    print(json.dumps({**{key: value for key, value in record.items() if key != 'bases'},
                      'bases': [{key: value for key, value in row.items() if key != 'cover'} for row in bases]}))


if __name__ == '__main__':
    main()
