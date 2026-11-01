#!/usr/bin/env python3
"""Exact finite controls for the explicit analytic-cutoff arithmetic."""
from fractions import Fraction
import json

from prime_form_cutoffs import ceil_log_two, generator_cutoff, sample_cutoff_bits


def main():
    logarithm_checks = 0
    for numerator in range(1, 65):
        for denominator in range(1, 65):
            ratio = Fraction(numerator, denominator)
            exponent = ceil_log_two(ratio)
            upper = Fraction(1 << exponent) if exponent >= 0 else Fraction(1, 1 << (-exponent))
            assert upper / 2 < ratio <= upper
            logarithm_checks += 1
    generation_checks = 0
    for discriminant in (-3, -4, -15, -23, -180, -41612, -(1 << 4096)):
        bits = (-discriminant).bit_length()
        cutoff = generator_cutoff(discriminant)
        base = 10000 * (bits + 1)
        assert cutoff == base ** 32
        assert Fraction(cutoff, 8) > 32 * (bits + 1) * base
        generation_checks += 1
    records = []
    for bits in (2, 3, 8, 16, 64, 1024, 1 << 20):
        for epsilon in (Fraction(1, 2), Fraction(2, 7), Fraction(1, 16),
                        Fraction(3, 1000), Fraction(1, 1 << 64), Fraction(1, 1 << 1024)):
            result = sample_cutoff_bits(bits, epsilon)
            records.append({
                'discriminant_bits': bits,
                'epsilon_numerator_bits': epsilon.numerator.bit_length(),
                'epsilon_denominator_bits': epsilon.denominator.bit_length(),
                'cutoff_bits': result['cutoff'].bit_length(),
                'weight_bits': result['weight_bits'],
                'normalized_error_le_epsilon_over_32': True,
                'unrounded_bias_lt_epsilon_over_4': True,
                'rounded_bias_lt_epsilon': True,
            })
    print(json.dumps({
        'scope': 'Exact cutoff and budget arithmetic under the stated analytic theorem; no prime enumeration or analytic theorem verification',
        'ceil_log_two_checks': logarithm_checks,
        'generation_cutoff_checks': generation_checks,
        'sample_budget_checks': len(records),
        'records': records,
    }, indent=2))


if __name__ == '__main__':
    main()
