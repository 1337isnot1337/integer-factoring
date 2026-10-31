#!/usr/bin/env python3
"""Exact global fixed-type dual certificates for unrestricted Gaussian gauges.

No optimizer output proves a sign. The checker recomputes every selected
minor valuation over F2(z), checks nonnegative rational dual stationarity,
and bounds the entropy ratio using outward integer logarithms.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from charged_gaussian_floor import Dyadic
from check_rectangular_defect_support import Z
from check_reflection_unrestricted_gaussian import closed_matrix, rational_smith
from reflection_dual_rail_gaussian import logical_heights


F = Fraction


def logarithm(arithmetic, value):
    """Outward log of a positive rational, reduced to log(1+x), 0<=x<1."""
    value = F(value)
    assert value > 0
    shift = value.numerator.bit_length()-value.denominator.bit_length()
    reduced = value/F(2)**shift
    while reduced < 1:
        reduced *= 2
        shift -= 1
    while reduced >= 2:
        reduced /= 2
        shift += 1
    small = reduced-1
    thermal = arithmetic.log_one_plus(arithmetic.rational(small.numerator, small.denominator))
    return arithmetic.add(thermal, arithmetic.multiply(arithmetic.rational(shift),
                                                       arithmetic.log_one_plus(arithmetic.one)))


def scalar_interval(arithmetic, value):
    value = F(value)
    return arithmetic.rational(value.numerator, value.denominator)


def entropy_interval(arithmetic, probability):
    probability = F(probability)
    assert 0 < probability < 1
    first = arithmetic.multiply(scalar_interval(arithmetic, probability),
                                logarithm(arithmetic, probability))
    second = arithmetic.multiply(scalar_interval(arithmetic, 1-probability),
                                 logarithm(arithmetic, 1-probability))
    return arithmetic.subtract(arithmetic.zero, arithmetic.add(first, second))


def profile_parameters(interval):
    heights = [h for place in ('zero', 'infinity') for h in logical_heights(interval, place)]
    energies = [4*(interval-1-2*i) for _ in range(2) for i in range(interval//2)]
    positive = [4*interval-2-4*i for i in range(interval)]
    return heights, energies, positive


def minor_coefficients(interval, side, row_mask, column_mask):
    modes = interval//2
    legs = (0, 1) if side == 0 else (0, 2)
    rows = [leg*modes+i for leg in legs for i in range(modes)]
    columns = [i for i in range(2*interval) if i not in rows]
    answer = [0]*(4*interval)
    for i in range(2*interval):
        if row_mask >> i & 1:
            answer[side*2*interval+rows[i % interval]] += 1 if i < interval else -1
        if column_mask >> i & 1:
            answer[side*2*interval+columns[i % interval]] += 1 if i < interval else -1
    return answer


def validate_certificate(certificate, arithmetic=None, determinant_checks=True):
    arithmetic = Dyadic(80) if arithmetic is None else arithmetic
    interval = certificate['interval']
    u = F(certificate['exp_negative_scale'])
    probabilities = [F(value) for value in certificate['probabilities']]
    assert interval >= 4 and interval & (interval-1) == 0
    assert 0 < u < 1 and len(probabilities) == interval
    assert all(0 < probability < 1 for probability in probabilities)
    natural, energies, positive = profile_parameters(interval)
    masses = [[F(0)]*interval for _ in range(2)]
    height_sum = [F(0)]*(4*interval)
    scale_sum = F(0)
    coefficients = [(1-u**order)/(2*(1+u**order)) for order in positive]+[F(0)]
    matrices = [closed_matrix(interval, side, infinity=bool(side))[0] for side in (0, 1)]
    if determinant_checks:
        expected = [-value for value in positive]+list(reversed(positive))
        for side in (0, 1):
            matrix, rows, columns = closed_matrix(interval, side, infinity=bool(side))
            heights = natural[side*2*interval:(side+1)*2*interval]
            rw = [heights[i] for i in rows]
            cw = [heights[i] for i in columns]
            rw += [-value for value in rw]
            cw += [-value for value in cw]
            gauged = [[value*Z**(rw[i]+cw[j]) for j, value in enumerate(row)]
                      for i, row in enumerate(matrix)]
            assert rational_smith(gauged) == expected
    seen = set()
    for minor in certificate['dual_minors']:
        side, size = minor['side'], minor['size']
        row_mask, column_mask = int(minor['row_mask'], 16), int(minor['column_mask'], 16)
        multiplier, valuation = F(minor['multiplier']), minor['valuation']
        assert side in (0, 1) and 1 <= size <= interval
        assert 0 < row_mask < 1 << (2*interval) and 0 < column_mask < 1 << (2*interval)
        assert row_mask.bit_count() == column_mask.bit_count() == size
        assert multiplier > 0 and isinstance(valuation, int)
        key = side, size, row_mask, column_mask
        assert key not in seen
        seen.add(key)
        linear = minor_coefficients(interval, side, row_mask, column_mask)
        # Activity is an exact equality at all positive multiples of the
        # natural integer height vector, not a floating near-equality.
        assert valuation+sum(a*b for a, b in zip(linear, natural)) == -sum(positive[:size])
        if determinant_checks:
            rows = [i for i in range(2*interval) if row_mask >> i & 1]
            columns = [i for i in range(2*interval) if column_mask >> i & 1]
            submatrix = [[matrices[side][i][j] for j in columns] for i in rows]
            assert sum(rational_smith(submatrix)) == valuation
        masses[side][size-1] += multiplier
        for i, value in enumerate(linear):
            height_sum[i] += multiplier*value
        scale_sum += multiplier*valuation
    for side in (0, 1):
        assert masses[side] == [coefficients[i]-coefficients[i+1] for i in range(interval)]
    assert scale_sum == -interval*interval
    for side in (0, 1):
        for leg in range(4):
            for mode in range(interval//2):
                index = side*2*interval+leg*(interval//2)+mode
                physical = (0 if leg in (0, 3) else 1)*(interval//2)+mode
                assert height_sum[index] == probabilities[physical]-F(1, 2)
    physical_mean = sum(energy*probability for energy, probability in zip(energies, probabilities))
    assert physical_mean == 2*sum(order*u**order/(1+u**order) for order in positive)
    entropy = arithmetic.zero
    for probability in probabilities:
        entropy = arithmetic.add(entropy, entropy_interval(arithmetic, probability))
    cost = arithmetic.zero
    for order in positive:
        cost = arithmetic.add(cost, arithmetic.multiply(arithmetic.rational(2),
            logarithm(arithmetic, 1+u**order)))
    cost = arithmetic.subtract(cost, arithmetic.multiply(
        scalar_interval(arithmetic, physical_mean), logarithm(arithmetic, u)))
    gap = arithmetic.subtract(cost, arithmetic.multiply(arithmetic.root_two, entropy))
    assert gap[0] > 0
    ratio = arithmetic.divide(cost, arithmetic.multiply(arithmetic.rational(2), entropy))
    return {'interval': interval, 'exact_minor_valuations': len(seen),
            'exact_stationarity_equations': 6*interval+1,
            'natural_bogoliubov_orders_checked': 4*interval,
            'twice_probability_entropy_interval': [str(value) for value in arithmetic.multiply(arithmetic.rational(2), entropy)],
            'global_minimum_cost_interval': [str(value) for value in cost],
            'critical_gap_interval': [str(value) for value in gap],
            'optimal_fixed_type_rho_interval': [str(value) for value in ratio],
            'interval_denominator': str(arithmetic.scale),
            'rho_display_only': float(F(sum(ratio), 2*arithmetic.scale))}


def verify(data):
    arithmetic = Dyadic(80)
    records = [validate_certificate(certificate, arithmetic) for certificate in data['certificates']]
    original = data['certificates'][0]
    altered_multiplier = deepcopy(original)
    altered_multiplier['dual_minors'][0]['multiplier'] = str(F(altered_multiplier['dual_minors'][0]['multiplier'])+1)
    missing_minor = deepcopy(original)
    missing_minor['dual_minors'] = missing_minor['dual_minors'][1:]
    altered_type = deepcopy(original)
    altered_type['probabilities'][0] = str(F(altered_type['probabilities'][0])+F(1, 1000))
    negative_controls = 0
    for altered in (altered_multiplier, missing_minor, altered_type):
        try:
            validate_certificate(altered, arithmetic, determinant_checks=False)
        except AssertionError:
            negative_controls += 1
        else:
            raise AssertionError('A malformed dual certificate was accepted.')
    root = Path(__file__).resolve().parent.parent
    sources = ['research/reflection-gaussian-fixed-type-duals.md',
               'tools/check_reflection_gaussian_fixed_type_duals.py',
               'tools/reflection_unrestricted_gaussian.py',
               'tools/check_reflection_unrestricted_gaussian.py']
    return {'status': 'pass', 'certificates': records,
            'exact_minor_valuations': sum(record['exact_minor_valuations'] for record in records),
            'rejected_negative_controls': negative_controls,
            'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in sources},
            'scope': 'Global minimum over all real independent Gaussian mode gauges and every positive scale for the listed fixed rational physical bit marginals; no all-type or improved-exponent claim.'}


def main():
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificates', type=Path, default=root/'evidence/reflection-gaussian-fixed-type-duals.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify(json.loads(args.certificates.read_text()))
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
