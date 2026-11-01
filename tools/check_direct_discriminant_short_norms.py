#!/usr/bin/env python3
"""Exact HNF, prime absorption and fixed-discriminant torsion controls."""

from functools import lru_cache
from math import gcd, isqrt, prod
import json
import random

from arithmetic_two_group import power, project_to_two_part, two_group_basis
from binary_quadratic_forms import (
    ambiguous_factor, compose_forms, inverse_form, reduce_form, reduced_forms,
)
from bounded_smooth_norm import smooth_factors
from quadratic_ideal_hnf import (
    compose_forms_hnf, form_to_ideal, ideal_to_form, multiply_ideals,
    primitive_ideal, reduce_with_transform, represented_norm_from_word,
)
from smooth_norm_words import norm_to_word, word_value


def prime_powers(integer):
    result = []
    candidate = 2
    while candidate * candidate <= integer:
        if integer % candidate == 0:
            value = 1
            while integer % candidate == 0:
                integer //= candidate
                value *= candidate
            result.append(value)
        candidate += 1
    if integer > 1:
        result.append(integer)
    return result


def unitary_pairs(integer):
    powers = prime_powers(integer)
    divisors = [1]
    for value in powers:
        divisors += [divisor * value for divisor in divisors]
    return [(value, integer // value) for value in sorted(divisors)
            if value <= integer // value]


def prime_forms(discriminant, cutoff):
    result = []
    for prime in range(2, cutoff + 1):
        if any(prime % divisor == 0 for divisor in range(2, isqrt(prime) + 1)):
            continue
        if discriminant % prime == 0:
            continue
        for root in range(2 * prime):
            if (root * root - discriminant) % (4 * prime) == 0:
                form = (prime, root, (root * root - discriminant) // (4 * prime))
                result.append((prime, root, reduce_form(form)))
    return result


def raw_prime_terms(form, bound, sign=1):
    factors = smooth_factors(form[0], bound)
    assert factors is not None
    return [(prime, sign * factors.count(prime), form[1] % (2 * prime))
            for prime in sorted(set(factors))]


def hnf_controls():
    discriminants = (
        -3, -4, -7, -8, -11, -15, -19, -20, -23, -24, -31, -39,
        -40, -43, -47, -56, -67, -68, -71, -84, -120, -163,
        -144, -180, -348, -420, -999, -1007, -1679, -4620, -10007,
    )
    products, transforms, word_controls, cancellations = 0, 0, 0, 0
    generator = random.Random(20261031)
    for discriminant in discriminants:
        forms = reduced_forms(discriminant)
        identity = next(form for form in forms if form[0] == 1)
        for left in forms:
            assert reduce_with_transform(ideal_to_form(form_to_ideal(left), discriminant))[0] == left
            for right in forms:
                assert compose_forms_hnf(left, right) == compose_forms(left, right)
                products += 1
            a, b, c = left
            for shift in (-19, -2, 3, 27):
                changed = (a + b * shift + c * shift * shift, b + 2 * c * shift, c)
                reduced, matrix = reduce_with_transform(changed)
                assert reduced == left
                for x, y in ((1, 0), (0, 1), (2, -3)):
                    xx = matrix[0][0] * x + matrix[0][1] * y
                    yy = matrix[1][0] * x + matrix[1][1] * y
                    assert a * x * x + b * x * y + c * y * y == (
                        changed[0] * xx * xx + changed[1] * xx * yy + changed[2] * yy * yy
                    )
                transforms += 1
        choices = prime_forms(discriminant, 19)
        if not choices:
            continue
        for _ in range(40):
            terms, expected = [], identity
            for _ in range(generator.randrange(1, 6)):
                prime, root, value = generator.choice(choices)
                exponent = generator.randrange(-6, 7)
                terms.append((prime, exponent, root))
                if exponent < 0:
                    value = inverse_form(value)
                expected = compose_forms_hnf(
                    expected, power(value, abs(exponent), identity, compose_forms_hnf)
                )
            representation = represented_norm_from_word(discriminant, terms)
            assert representation.form == expected
            assert representation.norm <= representation.product_norm
            converted = norm_to_word(expected, *representation.vector, 19)
            assert converted is not None
            assert word_value(converted, discriminant) == expected
            word_controls += 1
            cancellations += representation.removed_integer_content > 1
    large = represented_norm_from_word(-999, ((5, 1000, 1), (5, -999, 1)))
    assert large.norm == 5 and large.removed_integer_content == 5 ** 999
    return {
        "discriminants_checked": len(discriminants),
        "pair_products_compared_with_independent_implementation": products,
        "retained_reduction_matrices_checked": transforms,
        "supplied_prime_words_checked": word_controls,
        "words_with_nontrivial_integer_content": cancellations,
        "large_content_control": {
            "discriminant": -999,
            "terms": [[5, 1000, 1], [5, -999, 1]],
            "product_norm_bits": large.product_norm.bit_length(),
            "primitive_norm": large.norm,
            "removed_content_bits": large.removed_integer_content.bit_length(),
        },
    }


def absorption_controls():
    records = []
    count, cancellations, above_old_radius = 0, 0, 0
    for discriminant in (-84, -144, -999, -1007, -1679, -10007):
        radius = isqrt(-discriminant - 1) // 2
        smooth_bound, prime_cutoff = 7, 19
        budget = (radius * radius).bit_length()
        smooth_forms = [form for form in reduced_forms(discriminant)
                        if form[0] <= radius and gcd(form[0], discriminant) == 1
                        and smooth_factors(form[0], smooth_bound) is not None]
        # This constructor replay deliberately caps the number of supplied
        # pairs; it is not an exhaustive witness-search or density check.
        selected = smooth_forms[:6]
        samples = prime_forms(discriminant, prime_cutoff)
        local_count = local_cancellations = 0
        for left in selected:
            for right in selected:
                for prime, root, sample in samples:
                    terms = (raw_prime_terms(left, smooth_bound)
                             + raw_prime_terms(right, smooth_bound, -1)
                             + [(prime, -1, root)])
                    representation = represented_norm_from_word(discriminant, terms)
                    expected = compose_forms_hnf(
                        compose_forms_hnf(left, inverse_form(right)), inverse_form(sample)
                    )
                    assert representation.form == expected
                    assert representation.product_norm == prime * left[0] * right[0]
                    assert representation.norm <= prime_cutoff * radius * radius
                    assert gcd(representation.norm, discriminant) == 1
                    factors = smooth_factors(representation.norm, prime_cutoff)
                    assert factors is not None and len(factors) <= budget
                    direct = norm_to_word(expected, *representation.vector, prime_cutoff, budget)
                    assert direct is not None and word_value(direct, discriminant) == expected
                    local_count += 1
                    local_cancellations += representation.removed_integer_content > 1
                    above_old_radius += representation.norm > radius * radius
        count += local_count
        cancellations += local_cancellations
        records.append({
            "discriminant": discriminant,
            "radius": radius,
            "smooth_bound": smooth_bound,
            "prime_cutoff": prime_cutoff,
            "full_small_smooth_set_size": len(smooth_forms),
            "selected_small_smooth_forms": len(selected),
            "oriented_sample_primes": len(samples),
            "absorbed_norms_checked": local_count,
            "integer_content_removals": local_cancellations,
        })
    assert cancellations > 0 and above_old_radius > 0
    return {
        "absorbed_norms_checked": count,
        "integer_content_removals": cancellations,
        "absorbed_norms_above_old_R_squared_bound": above_old_radius,
        "records": records,
    }


def fixed_discriminant_controls():
    targets = list(range(3, 2002, 2)) + [10403, 10807, 11025, 16637, 25025]
    ambiguous_count, proper_count, maximum_rank = 0, 0, 0
    selected_records = []
    basis_targets = {15, 21, 45, 63, 75, 105, 225, 385, 1155, 1365,
                     10403, 10807, 11025, 16637, 25025}
    basis_records = []
    for integer in targets:
        discriminant = -4 * integer
        forms = reduced_forms(discriminant)
        identity = (1, 0, integer)
        ambiguous = {form for form in forms if inverse_form(form) == form}
        boundary = {form for form in forms if form[1] == 0 or form[1] == form[0] or form[0] == form[2]}
        assert ambiguous == boundary
        pairs = unitary_pairs(integer)
        zero_middle = {(left, 0, right) for left, right in pairs}
        assert zero_middle == {form for form in ambiguous if form[1] == 0}
        expected = set(zero_middle)
        if integer % 4 == 1:
            for left, right in pairs:
                assert right != 3 * left
                if right > 3 * left:
                    form = (2 * left, 2 * left, (left + right) // 2)
                else:
                    form = ((left + right) // 2, right - left, (left + right) // 2)
                assert reduce_form(form) == form
                expected.add(form)
        assert expected == ambiguous
        omega = len(prime_powers(integer))
        rank = omega - (integer % 4 == 3)
        assert len(ambiguous) == 1 << rank
        nonfactoring = {form for form in ambiguous if ambiguous_factor(form, integer) is None}
        expected_kernel = {identity}
        if integer % 4 == 1:
            expected_kernel.add((2, 2, (integer + 1) // 2))
        assert nonfactoring == expected_kernel
        assert all(compose_forms_hnf(left, right) in nonfactoring
                   for left in nonfactoring for right in nonfactoring)
        if omega >= 2:
            assert ambiguous - nonfactoring
        ambiguous_count += len(ambiguous)
        proper_count += len(ambiguous - nonfactoring)
        maximum_rank = max(maximum_rank, rank)
        if integer in basis_targets:
            multiply = lru_cache(None)(compose_forms_hnf)
            # Enumeration supplies this annihilator only for the finite
            # control; no efficient class-number algorithm is claimed.
            projected = project_to_two_part(forms, len(forms), identity, multiply)
            basis = two_group_basis(projected, identity, multiply, inverse_form, len(forms).bit_length())
            tops = [power(form, 1 << (order - 1), identity, multiply) for form, order in basis]
            assert len(tops) == rank
            factors = [ambiguous_factor(form, integer) for form in tops]
            assert any(factor is not None for factor in factors)
            basis_records.append({
                "integer": integer,
                "class_number": len(forms),
                "two_rank": rank,
                "basis_tops": tops,
                "proper_factors_from_tops": [factor for factor in factors if factor is not None],
            })
        if integer in (3, 5, 9, 15, 21, 45, 63, 225, 1155, 11025, 25025):
            selected_records.append({
                "integer": integer,
                "class_number": len(forms),
                "distinct_prime_factors_for_independent_check": omega,
                "two_rank": rank,
                "ambiguous_classes": len(ambiguous),
                "nonfactoring_subgroup": sorted(nonfactoring),
            })
    return {
        "integers_checked": len(targets),
        "complete_odd_interval": [3, 2001],
        "additional_targets": targets[1000:],
        "ambiguous_classes_checked": ambiguous_count,
        "proper_factor_classes_checked": proper_count,
        "maximum_two_rank": maximum_rank,
        "selected_classification_records": selected_records,
        "complete_two_basis_factor_controls": basis_records,
    }


def main():
    print(json.dumps({
        "scope": "Exact polynomial ideal-HNF constructors, absorbed supplied prime words, and finite complete two-torsion classification at Delta=-4N",
        "hnf_controls": hnf_controls(),
        "prime_absorption_controls": absorption_controls(),
        "fixed_discriminant_controls": fixed_discriminant_controls(),
        "fast_relation_or_smooth_norm_search_supplied": False,
        "new_factoring_exponent_established": False,
    }, indent=2))


if __name__ == "__main__":
    main()
