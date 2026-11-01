#!/usr/bin/env python3
"""Finite constructor checks; this does not test the Hecke zero-free theorem."""
from collections import deque
from fractions import Fraction
from functools import lru_cache
import json

from binary_quadratic_forms import compose_forms, inverse_form, reduced_forms
from prime_form_bias import (
    log_lower_dyadic, log_rational_bounds, prime_form_sample,
    uniform_distance_certificate,
)


def fraction_record(value):
    return {'numerator': value.numerator, 'denominator': value.denominator}


def generated_subgroup(identity, generators):
    multiply = lru_cache(None)(compose_forms)
    reached = {identity}
    pending = deque([identity])
    while pending:
        current = pending.popleft()
        for generator in generators:
            candidate = multiply(current, generator)
            if candidate not in reached:
                reached.add(candidate)
                pending.append(candidate)
    return reached


def check_discriminant(discriminant):
    forms = reduced_forms(discriminant)
    identity = next(form for form in forms if form[0] == 1)
    attempts = []
    for cutoff in (17, 37, 97, 257, 1021, 4093):
        masses, records = prime_form_sample(discriminant, cutoff)
        if not records:
            attempts.append({'cutoff': cutoff, 'empty': True})
            continue
        assert all(form in forms for form in masses)
        assert all(masses[inverse_form(form)] == mass for form, mass in masses.items())
        assert all(discriminant % prime for prime, _, _, _ in records)
        # Every source record is an actual prime form; no class-group
        # coordinates, factorization of Delta, or class-number input was used.
        for prime, middle, _, mass in records:
            assert (middle * middle - discriminant) % (4 * prime) == 0
            assert mass == log_lower_dyadic(prime, 12) * prime * (cutoff - prime)
        distance = uniform_distance_certificate(forms, masses)
        attempts.append({'cutoff': cutoff, 'l1_bound': fraction_record(distance)})
        if distance <= Fraction(3, 4):
            closure = generated_subgroup(identity, tuple(masses))
            assert closure == set(forms)
            # The l1 bound is a direct rational certificate of character bias;
            # it does not rely on floating Fourier transforms.
            return {
                'discriminant': discriminant,
                'class_number': len(forms),
                'cutoff': cutoff,
                'oriented_prime_records': len(records),
                'distinct_support_classes': len(masses),
                'total_mass_bits': sum(masses.values()).bit_length(),
                'certified_nonprincipal_bias_at_most': fraction_record(distance),
                'generated_group_order': len(closure),
                'attempts': attempts,
            }
    raise AssertionError(('finite cutoff list did not certify bias', discriminant, attempts))


def main():
    logarithm_checks = 0
    for integer in (1, 2, 3, 5, 7, 17, 31, 97, 257, 1021, 4093, (1 << 61) - 1):
        for bits in (1, 3, 12, 40):
            low, high = log_rational_bounds(integer, bits)
            lower = Fraction(log_lower_dyadic(integer, bits), 1 << bits)
            fine_low, fine_high = log_rational_bounds(integer, bits + 12)
            assert low <= fine_low <= fine_high <= high
            assert lower <= fine_low and fine_high - lower <= Fraction(2, 1 << bits)
            logarithm_checks += 1
    discriminants = (
        -3, -4, -7, -8, -11, -15, -16, -19, -20, -23, -24, -27,
        -28, -31, -36, -39, -40, -47, -48, -56, -63, -64, -68,
        -71, -75, -84, -100, -120, -144, -163, -180, -196, -225 * 4,
        -300, -400, -500, -675, -1155 * 4, -10403 * 4,
    )
    records = [check_discriminant(discriminant) for discriminant in discriminants]
    print(json.dumps({
        'scope': 'Exact prime-form and rational-weight controls; small all-class enumeration certifies displayed biases only; analytic theorem and uniform cutoff are not verified by this program',
        'logarithm_interval_checks': logarithm_checks,
        'discriminants_checked': len(records),
        'maximum_class_number': max(row['class_number'] for row in records),
        'maximum_cutoff': max(row['cutoff'] for row in records),
        'oriented_prime_records_checked': sum(row['oriented_prime_records'] for row in records),
        'records': records,
    }, indent=2))


if __name__ == '__main__':
    main()
