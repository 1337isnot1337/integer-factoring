#!/usr/bin/env python3
"""Actual-form counterexamples to replacing pair extraction by single hits."""
from collections import Counter
from math import gcd
import json

from binary_quadratic_forms import compose_forms, reduce_form
from prime_form_bias import prime_form_sample, primes_below


def smooth_remainder(integer, smooth_bound):
    remaining = integer
    for prime in (2,) + primes_below(smooth_bound + 1):
        while remaining % prime == 0:
            remaining //= prime
    return remaining


def recognition_failure(form, discriminant, radius, smooth_bound):
    leading = form[0]
    if leading > radius:
        return 'leading_coefficient_exceeds_radius'
    if gcd(leading, discriminant) != 1:
        return 'leading_coefficient_not_coprime'
    if smooth_remainder(leading, smooth_bound) != 1:
        return 'leading_coefficient_not_smooth'
    return None


def check_target(prime):
    assert prime in primes_below(prime + 1)
    discriminant = 1 - 4 * prime ** 4
    radius = prime ** 2 - 1
    target = (prime, 1, prime ** 3)
    assert reduce_form(target) == target
    assert 4 * radius ** 2 < -discriminant
    assert (prime - 1) * radius < prime ** 3
    masses, source_records = prime_form_sample(discriminant, prime, 8)
    assert masses  # Avoid vacuous examples with no split sampling prime.
    source_primes = {form: p for p, _, form, _ in source_records}
    reasons = Counter()
    prime_controls = []
    for form in masses:
        translated = compose_forms(target, form)
        assert translated[0] == prime * source_primes[form]
        reason = recognition_failure(translated, discriminant, radius, prime - 1)
        assert reason == 'leading_coefficient_not_smooth'
        reasons[reason] += 1
        prime_controls.append({'sample': form, 'translate': translated, 'failure': reason})
    word_count = 0
    word_values = set()
    for p, _, left, _ in source_records:
        for q, _, right, _ in source_records:
            if p * q * radius >= prime ** 3:
                continue
            word = compose_forms(left, right)
            translated = compose_forms(target, word)
            assert translated[0] == prime * word[0]
            reason = recognition_failure(translated, discriminant, radius, prime - 1)
            assert reason == 'leading_coefficient_not_smooth'
            word_count += 1
            word_values.add(word)
    coordinate_checks = 0
    for y in range(-4, 5):
        for x in range(-2 * prime, 2 * prime + 1):
            if gcd(x, y) != 1:
                continue
            value = prime * x * x + x * y + prime ** 3 * y * y
            if 0 < value < prime ** 3:
                assert y == 0 and abs(x) == 1 and value == prime
            coordinate_checks += 1
    return {
        'prime_parameter': prime,
        'discriminant': discriminant,
        'radius': radius,
        'smooth_bound': prime - 1,
        'target': target,
        'prime_sample_classes': len(masses),
        'oriented_prime_records': len(source_records),
        'single_translate_failures': dict(reasons),
        'two_prime_words_checked': word_count,
        'distinct_two_prime_word_values': len(word_values),
        'primitive_coordinate_checks': coordinate_checks,
        'prime_controls': prime_controls,
    }


def main():
    records = [check_target(prime) for prime in (13, 17, 19, 23, 29, 31, 43, 101, 257, 1009)]
    print(json.dumps({
        'scope': 'Exact finite arithmetic counterexamples to the stated single-smooth-form prime-translate rule; not a factoring lower bound',
        'targets_checked': len(records),
        'prime_translates_checked': sum(row['prime_sample_classes'] for row in records),
        'two_prime_words_checked': sum(row['two_prime_words_checked'] for row in records),
        'primitive_coordinate_checks': sum(row['primitive_coordinate_checks'] for row in records),
        'records': records,
    }, indent=2))


if __name__ == '__main__':
    main()
