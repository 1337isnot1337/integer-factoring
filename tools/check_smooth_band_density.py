#!/usr/bin/env python3
"""Exact checks for the multiband smooth-form density parameter theorem."""
from fractions import Fraction
import json
from math import factorial
from pathlib import Path

from prime_form_cutoffs import sample_cutoff_bits
from smooth_band_density import balanced_band_parameters
from smooth_band_density import largest_band_coefficient, parameter_summary


def kernel_checks():
    q = Fraction(4, 3)
    first = Fraction(21, 16)
    second = Fraction(27, 64)
    assert 1 - first * q + second * q * q == 0
    assert -1 + first * q * q - second * q ** 4 == 0
    assert 1 - first * q == -1 / q
    assert -1 + first * q * q == q
    integral = (1 - first / q + second / (q * q)) / 6
    assert integral == Fraction(259, 6144)
    assert integral - Fraction(1, 1000) > Fraction(1, 32)
    assert 1 + first + second == Fraction(175, 64) < 3

    def weight(value):
        return value * (1 - value) if 0 <= value <= 1 else 0

    for numerator in range(2049):
        value = Fraction(numerator, 2048)
        psi = (weight(value) - first * weight(q * value)
               + second * weight(q * q * value))
        expected = (0 if value <= Fraction(9, 16) else
                    q * value * (value - Fraction(9, 16))
                    if value <= Fraction(3, 4) else weight(value))
        assert psi == expected
        assert 0 <= psi <= Fraction(1, 4)
    tested_bits = [2, 3, 7, 32, 1024, 262144, 1000000, 100000000]
    for b in tested_bits:
        z = 1000000 * (b + 1)
        error = (Fraction(300 * (b + 1), z) + Fraction(48, z ** 15)
                 + Fraction(48 * (b + 1), z ** 31))
        assert error < Fraction(1, 1000)
    return {'exact_kernel_grid_points': 2049,
            'shell_error_discriminant_bits': tested_bits}


def coefficient_checks():
    summaries = []
    for s in range(1, 65):
        index, maximum, coefficients = largest_band_coefficient(s)
        assert len(coefficients) == s * (s - 1) + 1
        assert sum(coefficients) == s ** s
        assert maximum * s * s >= s ** s
        assert coefficients == list(reversed(coefficients))
        if s in (1, 2, 8, 32, 64):
            summaries.append({'bands': s, 'largest_index': index,
                              'largest_coefficient_bits': maximum.bit_length()})
    return {'tested_band_counts': 64, 'examples': summaries}


def weighted_selection_checks():
    """Count distinct weighted choices independently by a subset DP.

    Each band has 8s separate labelled objects, with weights 1, 2, or 3.
    The DP counts unordered subsets by size and sum of band indices.
    Multiplication by s! recovers the ordered distinct-choice mass.
    """
    examples = []
    for s in range(1, 9):
        index, pattern_count, _ = largest_band_coefficient(s)
        bands = [[1 + (item + band) % 3 for item in range(8 * s)]
                 for band in range(s)]
        least_mass = min(map(sum, bands))
        largest_weight = max(max(weights) for weights in bands)
        assert least_mass >= 2 * s * largest_weight
        # Independent bivariate generating function, one factor 1+w*x*z^j
        # for each labelled object. Decreasing subset sizes avoid reuse.
        counts = {(0, 0): 1}
        for band, weights in enumerate(bands):
            for weight in weights:
                additions = {}
                for (size, index_sum), count in counts.items():
                    if size < s and index_sum + band <= index:
                        key = (size + 1, index_sum + band)
                        additions[key] = additions.get(key, 0) + weight * count
                for key, count in additions.items():
                    counts[key] = counts.get(key, 0) + count
        ordered_mass = factorial(s) * counts.get((s, index), 0)
        direct_lower = pattern_count * (least_mass - (s - 1) * largest_weight) ** s
        half_mass_lower = pattern_count * least_mass ** s
        assert ordered_mass >= direct_lower
        assert ordered_mass * 2 ** s >= half_mass_lower
        examples.append({'bands': s, 'objects': 8 * s * s,
                         'ordered_mass_bits': ordered_mass.bit_length(),
                         'direct_lower_bound_bits': direct_lower.bit_length()})
    return {'distinct_choice_examples': examples}


def parameter_checks():
    cases = [2, 64, 1024, 65536, 131072, 262144, 1000000, 10000000, 100000000]
    summaries = []
    eligible_count = 0
    for b in cases:
        values = balanced_band_parameters(b)
        v = values['smoothness_exponent']
        r = values['radius_exponent']
        t = values['ceil_log_two_bits']
        assert 2 ** (t - 1) < b <= 2 ** t
        assert (v - 1) ** 2 < 8 * b * t <= v * v
        assert r == (b + 1) // 2 - 3
        if values['eligible']:
            eligible_count += 1
            q = values['eta_denominator']
            s = values['band_count']
            assert q >= 2
            assert v * 2 ** (q - 1) < r <= v * 2 ** q
            assert (s - 1) * q * v < (q + 1) * r <= s * q * v
            assert v >= 4 * s * q
            assert 1 << v >= (1000000 * (b + 1)) ** 128
            assert 1 << v >= (512 * s * s * q) ** 4
            assert 1 << (5 * v) >= (64 * s * v) ** 8
            denominator = (256 * (b + 2) * (s + 1) * s * s * v
                           * (256 * q) ** s * factorial(s))
            assert values['density_denominator'] == denominator
            epsilon = Fraction(1, 4 * denominator)
            independent = sample_cutoff_bits(b, epsilon)
            assert values['sample_base'] == independent['base']
            assert values['sample_cutoff'] == independent['cutoff']
            assert values['largest_prime_cutoff'] == max(1 << v, independent['cutoff'])
            assert values['norm_multiplicity'] == 2 * s + 3
            assert values['relation_diagonal'] == 2 * s + 4
            # All b-bit discriminants exceed 4R^2; compare exponents only.
            assert 2 * r + 2 < b - 1
        summaries.append(parameter_summary(values))
    assert eligible_count == 4
    return {'eligible_parameter_cases': eligible_count,
            'parameter_cases': summaries}


def main():
    evidence = {
        'claim_scope': 'Parameter and finite combinatorial checks; no norm finder.',
        'kernel': kernel_checks(),
        'band_coefficients': coefficient_checks(),
        'weighted_selection': weighted_selection_checks(),
        'parameters': parameter_checks(),
    }
    destination = Path(__file__).resolve().parents[1] / 'evidence' / 'smooth-band-density.json'
    destination.write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps({
        'status': 'pass',
        'kernel_grid_points': evidence['kernel']['exact_kernel_grid_points'],
        'band_counts': evidence['band_coefficients']['tested_band_counts'],
        'weighted_selection_examples': len(evidence['weighted_selection']['distinct_choice_examples']),
        'eligible_parameter_cases': evidence['parameters']['eligible_parameter_cases'],
        'evidence': str(destination),
    }, indent=2))


if __name__ == '__main__':
    main()
