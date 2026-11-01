#!/usr/bin/env python3
"""Finite combinatorial and exact parameter checks for powered samples."""
import cmath
from fractions import Fraction
import json
from math import comb, pi
from pathlib import Path

from powered_prime_norm_parameters import powered_prime_parameters
from powered_prime_norm_parameters import powered_parameter_summary
from prime_form_cutoffs import sample_cutoff_bits


def binomial_checks():
    cases = 0
    for population in range(1, 129):
        for selected in range(1, population + 1):
            assert comb(population, selected) * selected ** selected >= population ** selected
            cases += 1
    return cases


def convolution_checks():
    """Exact convolution followed by a diagnostic complex Fourier check.

    The distributions and convolution coefficients are rational. The
    Fourier identity is also checked by multiplication in the exact
    quotient ring Q[z]/(z^m-1), before evaluating any complex numbers.
    """
    examples = []
    for modulus in range(3, 13):
        weights = [1 + (index * index + 3 * index) % 7 for index in range(modulus)]
        total = sum(weights)
        distribution = [Fraction(weight, total) for weight in weights]
        current = [Fraction(int(index == 0)) for index in range(modulus)]
        fourier = [sum(distribution[index] * cmath.exp(2j * pi * frequency * index / modulus)
                       for index in range(modulus)) for frequency in range(modulus)]
        initial_bias = max(abs(value) for value in fourier[1:])
        for steps in range(1, 7):
            convolution = [sum(current[index] * distribution[(target - index) % modulus]
                               for index in range(modulus)) for target in range(modulus)]
            # Independent polynomial product and reduction modulo z^m-1.
            polynomial = [Fraction(0) for _ in range(2 * modulus - 1)]
            for left, left_weight in enumerate(current):
                for right, right_weight in enumerate(distribution):
                    polynomial[left + right] += left_weight * right_weight
            quotient = [Fraction(0) for _ in range(modulus)]
            for degree, coefficient in enumerate(polynomial):
                quotient[degree % modulus] += coefficient
            assert convolution == quotient
            assert sum(convolution) == 1 and min(convolution) >= 0
            actual = [sum(convolution[index] * cmath.exp(2j * pi * frequency * index / modulus)
                          for index in range(modulus)) for frequency in range(modulus)]
            assert all(abs(actual[frequency] - fourier[frequency] ** steps) < 1e-12
                       for frequency in range(modulus))
            assert max(abs(value) for value in actual[1:]) <= initial_bias ** steps + 1e-12
            current = convolution
        examples.append({'cyclic_order': modulus, 'powers_checked': 6})
    return {'exact_rational_convolution_examples': examples,
            'complex_checks': 'Diagnostic only; tolerance 1e-12, not used by parameter proofs.'}


def parameter_checks():
    examples = []
    eligible = 0
    for bits in [2, 1024, 1000000, 10000000, 100000000, 10000000000]:
        values = powered_prime_parameters(bits)
        v = values['smoothness_exponent']
        w = values['root_base_exponent']
        assert v == 32 * w
        assert 2048 * (w - 1) ** 2 < bits <= 2048 * w * w
        if values['eligible']:
            eligible += 1
            k = values['smooth_form_multiplicity']
            r = values['radius_exponent']
            h = values['sample_bias_bits']
            density = values['density_denominator']
            steps = values['walk_steps']
            assert k * v <= r < (k + 1) * v
            assert density == (bits + 2) * 2 ** (((bits + 1) // 2) - k * v) * (32 * v * k) ** k
            assert 4 * density <= 1 << (h * steps)
            assert 1 << (h * (steps - 1)) < 4 * density
            assert values['sample_cutoff'] <= 1 << v
            independent = sample_cutoff_bits(bits, Fraction(1, 1 << h))
            assert values['sample_cutoff'] == independent['cutoff']
            assert values['norm_multiplicity'] == 2 * k + steps
            assert values['relation_diagonal'] == values['norm_multiplicity'] + 1
            assert values['accepted_norm_log_two_bound'] <= values['guaranteed_norm_log_two_bound']
            assert values['guaranteed_norm_log_two_bound'] <= bits + steps * v
            assert 2 * r + 2 < bits - 1
        examples.append(powered_parameter_summary(values))
    assert eligible == 3

    custom = []
    for width in [4096, 8192, 16384, 32768]:
        values = powered_prime_parameters(100000000, width)
        assert values['eligible']
        assert values['sample_cutoff'] <= 1 << width
        assert 4 * values['density_denominator'] <= 1 << (
            values['sample_bias_bits'] * values['walk_steps'])
        custom.append(powered_parameter_summary(values))
    return {'eligible_balanced_cases': eligible, 'balanced_cases': examples,
            'custom_width_cases': custom}


def main():
    evidence = {
        'claim_scope': 'Existence parameters and decoding costs; no target-norm selector.',
        'binomial_inequality_cases': binomial_checks(),
        'convolution': convolution_checks(),
        'parameters': parameter_checks(),
    }
    destination = Path(__file__).resolve().parents[1] / 'evidence' / 'powered-prime-norm-parameters.json'
    destination.write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps({'status': 'pass',
                      'binomial_cases': evidence['binomial_inequality_cases'],
                      'cyclic_convolution_cases': 60,
                      'eligible_balanced_cases': evidence['parameters']['eligible_balanced_cases'],
                      'evidence': str(destination)}, indent=2))


if __name__ == '__main__':
    main()
