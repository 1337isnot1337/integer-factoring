#!/usr/bin/env python3
"""Exact replay for the canonical alternating Cauchy-packet family.

The main path uses integer root products and exact linear division. The
small complete contractions use the separate rational packet implementation.
No floating-point arithmetic or numerical rank calculation is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from cauchy_packet_polynomials import modular_packet_polynomial


def valuation(value):
    if not value:
        raise ValueError('Zero has no finite valuation.')
    value = abs(value)
    return (value & -value).bit_length() - 1


def multiply_linear(polynomial, constant, linear):
    result = [0] * (len(polynomial) + 1)
    for degree, value in enumerate(polynomial):
        result[degree] += constant * value
        result[degree + 1] += linear * value
    return result


def exact_divide_linear(polynomial, constant, linear):
    remainder = polynomial[:]
    quotient = [0] * (len(polynomial) - 1)
    for degree in range(len(quotient), 0, -1):
        value, rest = divmod(remainder[degree], linear)
        assert not rest
        quotient[degree - 1] = value
        remainder[degree] -= linear * value
        remainder[degree - 1] -= constant * value
    assert not any(remainder)
    return quotient


def row_contents(length, center, state=None):
    if state is None:
        state = sum(1 << x for x in range(0, length, 2))
    dilation = 2 * length
    denominator = 1
    polynomial = [1]
    for root in range(-length, length):
        polynomial = multiply_linear(polynomial, center - root, dilation)
        if root:
            denominator *= -root
    denominator_valuation = valuation(denominator)
    kernels = {
        difference: exact_divide_linear(polynomial, center - difference, dilation)
        for difference in range(-length + 1, length)
    }
    result = []
    for x in range(length):
        numerator = [sum(kernels[y - x][degree]
                         for y in range(length) if state >> y & 1)
                     for degree in range(2 * length)]
        if not (state >> x & 1):
            numerator = [-value for value in numerator]
            numerator[0] += denominator
        if center == 0:
            assert numerator[0] == denominator
        nonzero = [valuation(value) for value in numerator if value]
        assert nonzero, (length, center, state, x)
        result.append(min(nonzero) - denominator_valuation)
    return result


def check():
    records = []
    row_checks = 0
    center_checks = 0
    expanded_checks = 0
    for length in (4, 8, 16, 32, 64):
        a = length.bit_length() - 1
        n = length // 2
        claimed = length * (a + 1) - 1
        profiles = []
        lower_degree_profiles = []
        state = sum(1 << x for x in range(0, length, 2))
        assert sum(row_contents(length, 0)) == 0
        for center in range(-length + 1, length, 2):
            actual = row_contents(length, center)
            assert min(actual) >= 0
            for x in range(0, length, 2):
                value = center + x - n + 1
                clipped = min(valuation(value), a + 1) if value else a + 1
                assert actual[x] == 2 * a + clipped, (length, center, x, actual[x])
                row_checks += 1
            total = sum(actual)
            assert total >= claimed, (length, center, total, claimed)
            profiles.append([center, total])
            lower_degree_profiles.append([center, total - actual[-1]])
            assert total - actual[-1] >= claimed
            center_checks += 1
            if length <= 8:
                column = [0] * (1 << length)
                column[state] = 1
                precision = total + 2
                full = modular_packet_polynomial(length, center, precision, state, column)
                full += [0] * (length * (2 * length - 1) + 1 - len(full))
                observed = min(valuation(value) for value in full if value)
                assert observed == total, (length, center, observed, total)
                expanded_checks += len(full)
        assert profiles[-1] == [length - 1, claimed]
        assert min(value for _, value in profiles) == claimed
        assert lower_degree_profiles[-1] == [length - 1, claimed]
        records.append({'L': length, 'creation_precision': 1, 'original_value': 1,
                        'minimum_content': claimed, 'profiles': profiles,
                        'dropped_last_row_minimum_content': min(value for _, value in lower_degree_profiles)})

    # The theorem deliberately starts at L=4 and controls only odd centers.
    assert sum(row_contents(2, 1)) == 2
    assert sum(row_contents(8, 2)) < 31
    assert sum(row_contents(8, 1, 84)) < 31

    root = Path(__file__).resolve().parent.parent
    sources = ['tools/check_canonical_alternating_packet.py',
               'tools/cauchy_packet_polynomials.py',
               'research/canonical-alternating-packet.md']
    return {'status': 'pass', 'occupied_row_formulas_checked': row_checks,
            'odd_centers_checked': center_checks,
            'independent_expanded_coefficients_checked': expanded_checks,
            'negative_controls': 3, 'records': records,
            'source_sha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                              for name in sources},
            'scope': 'Finite exact replay of an all-size proof for canonical binary packets. '
                     'Only odd centers are controlled; no factoring exponent improvement is claimed.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        result = check()
    except Exception as error:
        if args.output:
            args.output.write_text(json.dumps({'status': 'fail', 'error': repr(error)}, indent=2) + '\n')
        raise
    rendered = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
