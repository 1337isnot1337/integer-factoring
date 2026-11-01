#!/usr/bin/env python3
"""Exact algebra and geometric controls for asymmetric crossing pressure.

No supplied profile is claimed to pass the supporting-line threshold.
The interval examples are conditional geometric implications only.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import lcm
from pathlib import Path


def encode(value):
    return [value.numerator, value.denominator]


def profile_sum(orders, base):
    return sum((base ** (-order) for order in orders), Fraction())


def rational_pressure(orders_infinity, orders_zero, gx, gy, *,
                      triple=(3, 4, 5), target=Fraction(7, 5), base=Fraction(3, 2)):
    u, v, c = triple
    assert min(triple) > 0 and u * u + v * v == c * c
    e, f = target.numerator, target.denominator
    assert target > 0 and base > 1
    first = profile_sum(orders_infinity, base ** (v * c * e))
    second = profile_sum(orders_zero, base ** (u * c * e))
    weights_x = [base ** (-u * v * f * g) for g in gx]
    weights_y = [base ** (-u * v * f * g) for g in gy]
    zx, zy = sum(weights_x), sum(weights_y)
    left = first ** (u * f) * second ** (v * f)
    right = (zx * zy) ** (c * e)
    return left, right, first, second, weights_x, weights_y


def gibbs_identity():
    """Exponentiate the full entropy identity, using arbitrary rational types."""
    triple, target, base = (3, 4, 5), Fraction(7, 5), Fraction(3, 2)
    infinity, zero = (-1, 0, 1, 2), (0, 1, 3)
    gx, gy = (0, 1), (-1, 0, 2)
    px = (Fraction(1, 3), Fraction(2, 3))
    py = (Fraction(1, 2), Fraction(1, 3), Fraction(1, 6))
    u, v, c = triple
    e, f = target.numerator, target.denominator
    left_bound, right_bound, first, second, wx, wy = rational_pressure(
        infinity, zero, gx, gy, triple=triple, target=target, base=base)
    optimal_x = [weight / sum(wx) for weight in wx]
    optimal_y = [weight / sum(wy) for weight in wy]
    assert sum(optimal_x) == sum(optimal_y) == 1
    length = lcm(*(p.denominator for p in px + py))
    mean = sum(p * g for p, g in zip(px, gx)) + sum(p * g for p, g in zip(py, gy))
    power = u * v * e * c * f * length * mean
    assert power.denominator == 1
    left = first ** (u * f * length) * second ** (v * f * length) * base ** int(power)
    right = (left_bound / right_bound) ** length
    for probability, optimal in zip(px + py, optimal_x + optimal_y):
        exponent = e * c * length * probability
        assert exponent.denominator == 1
        left *= probability ** int(exponent)
        right *= (probability / optimal) ** int(exponent)
    assert left == right
    return {'profiles_infinity': list(infinity), 'profiles_zero': list(zero),
            'total_heights_x': list(gx), 'total_heights_y': list(gy),
            'weights_triple': list(triple), 'target': encode(target), 'base': encode(base),
            'test_type_x': list(map(encode, px)), 'test_type_y': list(map(encode, py)),
            'gibbs_x': list(map(encode, optimal_x)), 'gibbs_y': list(map(encode, optimal_y)),
            'entropy_identity_equal': True,
            'comparison_integer_bits': max(left.numerator.bit_length(), left.denominator.bit_length()),
            'formal_profile_control_only': True}


def constant_crossing_control():
    """The actual constant 2-by-3 crossing has complete ranks six in both cuts."""
    from check_paired_type_crossing_pressure import flatten, rank
    qx, qy, count = 2, 3, 6
    values = [0] * (qx * qx * qy * qy)
    for x in range(qx):
        for y in range(qy):
            values[x + qx * (y + qy * (y + qy * x))] = 1
    assert [rank(flatten(values, qx, qy, cut)) for cut in (0, 1)] == [count, count]
    records = []
    for target, expected in ((Fraction(4, 3), 1), (Fraction(7, 5), 0), (Fraction(3, 2), -1)):
        original = rational_pressure([0] * count, [0] * count, [0] * qx, [0] * qy,
                                     target=target)
        left, right = original[:2]
        observed = (left > right) - (left < right)
        assert observed == expected
        # K(w)=w^d times the same crossing is still a strict polynomial path.
        for shift in (1, 7, 31):
            shifted = rational_pressure([-shift] * count, [shift] * count,
                                        [0] * qx, [0] * qy, target=target)
            assert shifted[:2] == original[:2]
        records.append({'target': encode(target), 'left_vs_right': observed,
                        'monomial_shifts_checked': [1, 7, 31]})
    return {'alphabets': [qx, qy, qy, qx], 'complete_cut_ranks': [count, count],
            'optimal_support_value': [7, 5], 'weight_sum': [7, 5], 'comparisons': records}


def integer_root(value, degree):
    assert value >= 0 and degree > 0
    if value < 2:
        return value
    low, high = 0, 1 << ((value.bit_length() + degree - 1) // degree)
    while high - low > 1:
        middle = (low + high) // 2
        if middle ** degree <= value:
            low = middle
        else:
            high = middle
    assert low ** degree <= value < (low + 1) ** degree
    return low


def rational_power_interval(base, exponent, denominator=10 ** 12):
    assert base > 0 and exponent > 0
    numerator, degree = exponent.numerator, exponent.denominator
    powered = base ** numerator
    scaled = powered * denominator ** degree
    root = integer_root(scaled.numerator // scaled.denominator, degree)
    low, high = Fraction(root, denominator), Fraction(root + 1, denominator)
    assert low ** degree <= powered < high ** degree
    return low, high


def geometry_controls():
    a, b = Fraction(3, 5), Fraction(4, 5)
    theta = a * a
    assert a * a + b * b == 1 and 1 - theta == b * b
    target = Fraction(9, 10)
    row, column = target * a, target * b
    assert a * row + b * column == target
    assert row / a == column / b == target

    # Exact fifth-power certificate of a visible positive margin.
    lower_row, lower_column = Fraction(27, 50), Fraction(19, 25)
    assert lower_row ** 5 < theta ** 3
    assert lower_column ** 5 < (1 - theta) ** 3
    lower = a * lower_row + b * lower_column
    assert lower == Fraction(233, 250) and target / lower == Fraction(225, 233) < 1

    # The same implication for a one-percent gap and exponent 51/100.
    close_target, exponent = Fraction(99, 100), Fraction(51, 100)
    rx = rational_power_interval(theta, exponent)
    cy = rational_power_interval(1 - theta, exponent)
    dl, du = a * rx[0] + b * cy[0], a * rx[1] + b * cy[1]
    assert close_target < dl <= du < 1
    return {'weights': [encode(a), encode(b)], 'aspect': encode(theta),
            'half_exponent_equalizer': {'target': encode(target), 'row': encode(row),
                                       'column': encode(column), 'exact': True},
            'fifth_power_margin': {'target': encode(target), 'exponent': [3, 5],
                                  'denominator_lower': encode(lower),
                                  'orientation_factor_upper': [225, 233]},
            'close_margin': {'target': encode(close_target), 'exponent': encode(exponent),
                             'row_power_interval': list(map(encode, rx)),
                             'column_power_interval': list(map(encode, cy)),
                             'denominator_interval': [encode(dl), encode(du)],
                             'orientation_factor_upper': encode(close_target / dl),
                             'unit_target_rejected': du < 1},
            'conditional_geometry_only_no_seed_claim': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent
    result = {'status': 'pass', 'scope': 'exact supporting-line algebra and conditional geometry',
              'gibbs_identity': gibbs_identity(), 'constant_crossing': constant_crossing_control(),
              'geometry': geometry_controls(),
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest()
                                for name in ('check_asymmetric_support_line.py',
                                             'check_paired_type_crossing_pressure.py')}}
    serialized = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(json.dumps({'status': result['status'], 'gibbs_identity_equal': True,
                      'constant_crossing_support': [7, 5],
                      'close_gap_denominator': result['geometry']['close_margin']['denominator_interval']}))


if __name__ == '__main__':
    main()
