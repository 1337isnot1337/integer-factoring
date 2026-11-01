#!/usr/bin/env python3
"""Exact theta-section coefficients, charged states, and common-space witnesses."""

import argparse
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

from lambert_theta_sections import (
    AdditivePhaseAlgebra, coefficient, q3_coefficient, square_section_rank_witness,
    theta_series,
)

ROOT = Path(__file__).resolve().parents[1]


def literal_signed_residues(integer, prime):
    values = [0]*prime
    for divisor in range(1, isqrt(integer)+1):
        if integer % divisor:
            continue
        for d in {divisor, integer//divisor}:
            values[d % prime] ^= 1
            values[-d % prime] ^= 1
    return tuple(values)


def circular_product(left, right, prime):
    """Independent cyclic convolution, followed by Phi_q projection."""
    value = 0
    for i in range(prime):
        for j in range(prime):
            if left >> i & 1 and right >> j & 1:
                value ^= 1 << ((i+j) % prime)
    return value ^ (((1 << prime)-1) if value >> (prime-1) & 1 else 0)


def check_series_identity():
    cases = 0
    for prime in (3, 5, 7, 11, 13):
        algebra = AdditivePhaseAlgebra(prime)
        limit = 64
        denominator = theta_series(algebra, limit)
        numerator = theta_series(algebra, limit, derivative=True)
        series = {0: circular_product(algebra.phase(1), algebra.inverse_theta_constant(1), prime)}
        for integer in range(1, limit+1):
            literal = literal_signed_residues(integer, prime)
            mask = sum(value << k for k, value in enumerate(literal))
            series[integer] = mask ^ (algebra.modulus if mask >> (prime-1) & 1 else 0)
        for integer in range(limit+1):
            value = 0
            for k, a in denominator.items():
                if k <= integer:
                    value ^= circular_product(a, series[integer-k], prime)
            assert value == numerator.get(integer, 0)
            cases += 1
    return cases


def check_coefficients():
    cases = products = 0
    for prime in (3, 5, 7, 11, 13, 17, 19, 23):
        for integer in range(1, 258):
            result = coefficient(integer, prime)
            assert result['signed_residues'] == literal_signed_residues(integer, prime)
            if prime == 3:
                assert result['value'] == q3_coefficient(integer)
            cases += 1
            products += result['ring_products']
    return {"cases": cases, "prime_levels": [3, 5, 7, 11, 13, 17, 19, 23],
            "integer_range": [1, 257], "charged_ring_products": products}


def check_rank_witnesses():
    records = []
    cross_checks = literal_checks = engine_checks = 0
    for k in range(4, 15):
        witness = square_section_rank_witness(k)
        modulus, roots = witness['modulus'], witness['roots']
        assert witness['dimension'] >= witness['uniform_lower_bound']
        assert len(set(witness['rows'])) == len(roots)
        assert len(set(witness['columns'])) == len(roots)
        for x, row in zip(roots, witness['rows']):
            for y, column in zip(roots, witness['columns']):
                integer = modulus*row+column
                assert integer % 3 == 1
                value = int(isqrt(integer)**2 == integer)
                assert value == (x == y)
                cross_checks += 1
                if k <= 7:
                    assert literal_signed_residues(integer, 3) == (0, value, value)
                    literal_checks += 1
                if k <= 6:
                    assert coefficient(integer, 3)['value'] == value
                    engine_checks += 1
        records.append(witness)
    return {"witnesses": records, "all_identity_entries_checked": cross_checks,
            "literal_divisor_entries": literal_checks, "theta_engine_entries": engine_checks}


def check_large_closed_controls():
    cases = []
    for two_exponent, five_exponent, three_exponent in ((64, 0, 0), (91, 23, 17), (128, 79, 47)):
        square = (1 << two_exponent)*5**five_exponent
        for extra in (1, 7):
            integer = square*square * 3**three_exponent * extra
            expected = int(extra == 1)
            assert q3_coefficient(integer) == expected
            cases.append({"integer": str(integer), "bits": integer.bit_length(),
                          "value": expected,
                          "known_prime_exponents": {"2": 2*two_exponent, "3": three_exponent,
                                                    "5": 2*five_exponent, "7": int(extra == 7)}})
    return cases


def profiles(bits):
    results = []
    for b in bits:
        integer = (1 << b)-1
        for prime in (3, 5, 7, 11):
            result = coefficient(integer, prime)
            assert result['signed_residues'] == literal_signed_residues(integer, prime)
            results.append({"bits": b, "integer": integer, "prime": prime, **result})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile-bits', type=int, nargs='*', default=[8, 10, 12])
    args = parser.parse_args()
    result = {
        "scope": "Exact theta-section engine and common-linear-space lower bound; no improved factoring rate",
        "series_identity_coefficients": check_series_identity(),
        "coefficient_controls": check_coefficients(),
        "rank_controls": check_rank_witnesses(),
        "large_level_three_controls": check_large_closed_controls(),
        "charged_profiles": profiles(args.profile_bits),
        "source_sha256": {name: sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('tools/lambert_theta_sections.py',
                                       'tools/check_lambert_theta_sections.py',
                                       'tools/cyclotomic_character_decoder.py')},
    }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
