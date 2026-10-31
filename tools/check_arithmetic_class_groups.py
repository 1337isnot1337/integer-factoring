#!/usr/bin/env python3
"""Actual quadratic-form controls; all-class enumeration is explicitly charged."""
from functools import lru_cache
from math import prod
import json

from arithmetic_two_group import power, project_to_two_part, two_group_basis
from binary_quadratic_forms import (
    ambiguous_factor, compose_forms, inverse_form, reduced_forms, reduce_form, smooth_prime_word,
)


def check_discriminant(discriminant, integer=None):
    forms = reduced_forms(discriminant)
    elements = set(forms)
    identity = next(form for form in forms if form[0] == 1)
    multiply = lru_cache(None)(compose_forms)
    inverse = inverse_form
    smooth_count = 0
    smooth_word_length = 0
    for form in forms:
        assert reduce_form(form) == form
        assert inverse(form) in elements and multiply(form, inverse(form)) == identity
        a, b, c = form
        for shift in (-3, -1, 1, 4):
            assert reduce_form((a+b*shift+c*shift*shift, b+2*c*shift, c)) == form
        assert reduce_form((c, -b, a)) == form
        for other in forms:
            assert multiply(form, other) in elements
            assert multiply(form, other) == multiply(other, form)
            for third in forms:
                assert multiply(multiply(form, other), third) == multiply(form, multiply(other, third))
        word = smooth_prime_word(form, 7)
        if word is not None:
            value = identity
            for prime_form, exponent in word:
                value = multiply(value, power(prime_form if exponent > 0 else inverse(prime_form),
                                              abs(exponent), identity, multiply))
            assert value == form
            smooth_count += 1
            smooth_word_length = max(smooth_word_length, sum(abs(exponent) for _, exponent in word))
    # Supplying every class makes the finite generator claim verifiable. This
    # is exponential enumeration in the input bit length, not the open fast
    # generator primitive. Orders are found by finite cyclic enumeration.
    generators = [form for form in forms if form != identity]
    orders = []
    for form in generators:
        current = form
        order = 1
        while current != identity:
            current = multiply(current, form)
            order += 1
            assert order <= len(forms)
        orders.append(order)
    # The diagonal matrix of individual orders is a full-rank relation
    # matrix. Its index is generally much larger than the actual class number.
    annihilator = prod(orders)
    assert all(power(form, annihilator, identity, multiply) == identity for form in forms)
    projected = project_to_two_part(generators, annihilator, identity, multiply)
    basis = two_group_basis(projected, identity, multiply, inverse, len(forms).bit_length())
    tops = [power(form, 1 << (order-1), identity, multiply) for form, order in basis]
    ambiguous = {form for form in forms if multiply(form, form) == identity}
    span = {identity}
    for form in tops:
        span |= {multiply(value, form) for value in span}
    assert span == ambiguous
    assert all(b == 0 or b == a or a == c for a, b, c in ambiguous)
    factor = None
    if integer is not None:
        factor = next((factor for form in tops for factor in [ambiguous_factor(form, integer)] if factor), None)
        assert factor is not None and integer % factor == 0 and 1 < factor < integer
    return {'discriminant': discriminant, 'class_number': len(forms),
            'two_rank': len(tops), 'two_part_order': 1 << sum(order for _, order in basis),
            'relation_determinant_bits': annihilator.bit_length(),
            'relation_index_exceeds_class_number': annihilator > len(forms),
            'factor_target': integer, 'recovered_factor': factor,
            'smooth_bound': 7, 'smooth_forms_checked': smooth_count,
            'maximum_smooth_word_length': smooth_word_length,
            'pair_products_checked': len(forms)**2, 'associative_triples_checked': len(forms)**3}


def main():
    records = [check_discriminant(value) for value in (
        -3, -4, -7, -8, -11, -15, -19, -20, -23, -24, -31, -39,
        -40, -43, -47, -56, -67, -68, -71, -84, -120, -163,
    )]
    targets = (15, 21, 33, 35, 45, 55, 65, 75, 77, 105, 165, 231, 385, 1155, 1365, 10403, 10807, 16637)
    records += [check_discriminant(-4*integer, integer) for integer in targets]
    assert all(record['class_number'] == 1 for record in records[:5])
    assert next(row for row in records if row['discriminant'] == -23)['class_number'] == 3
    print(json.dumps({
        'scope': 'Exact small quadratic-form class groups; generator/order search enumerates all classes and is not a fast factoring algorithm',
        'discriminants_checked': len(records),
        'pair_products_checked': sum(row['pair_products_checked'] for row in records),
        'associative_triples_checked': sum(row['associative_triples_checked'] for row in records),
        'factoring_targets_checked': len(targets),
        'maximum_class_number': max(row['class_number'] for row in records),
        'records': records,
    }, indent=2))


if __name__ == '__main__':
    main()
