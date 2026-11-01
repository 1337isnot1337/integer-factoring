#!/usr/bin/env python3
"""Check all-representation allocation matching against complete ellipses."""

import argparse
from dataclasses import asdict
from hashlib import sha256
from itertools import product
import json
from math import gcd, isqrt, prod
from pathlib import Path

from arbitrary_norm_representation import (
    balanced_radix_partition, represent_arbitrary_factored_norm,
)
from binary_quadratic_forms import reduced_forms
from quadratic_ideal_hnf import represented_norm_from_word


def primes_through(bound):
    return [p for p in range(2, bound+1)
            if all(p % divisor for divisor in range(2, isqrt(p)+1))]


def base_for(discriminant, bound):
    split, inert = [], []
    for prime in primes_through(bound):
        if gcd(prime, discriminant) != 1:
            continue
        roots = [root for root in range(2*prime)
                 if (root*root-discriminant) % (4*prime) == 0]
        if roots:
            split.append((prime, roots[0]))
        elif prime != 2:
            assert pow(discriminant % prime, (prime-1)//2, prime) == prime-1
            inert.append((prime, None))
    return split, inert


def ellipse_status(form, norm):
    """Return existence and primitive existence by exhaustive integer rows."""
    a, middle, c = form
    absolute_discriminant = 4*a*c-middle*middle
    limit = isqrt((4*a*norm)//absolute_discriminant)
    any_found = False
    for y in range(-limit, limit+1):
        square = 4*a*norm-absolute_discriminant*y*y
        root = isqrt(square)
        if root*root != square:
            continue
        for signed in (root, -root):
            numerator = signed-middle*y
            if numerator % (2*a):
                continue
            x = numerator//(2*a)
            assert a*x*x+middle*x*y+c*y*y == norm
            any_found = True
            if gcd(x, y) == 1:
                return True, True
    return any_found, False


def smallest_primitive_divisor(form, norm):
    """Minimize n/gcd(x,y)^2 over every representation of the fixed n."""
    a, middle, c = form
    absolute_discriminant = 4*a*c-middle*middle
    limit = isqrt((4*a*norm)//absolute_discriminant)
    smallest = None
    for y in range(-limit, limit+1):
        square = 4*a*norm-absolute_discriminant*y*y
        root = isqrt(square)
        if root*root != square:
            continue
        for signed in (root, -root):
            numerator = signed-middle*y
            if numerator % (2*a):
                continue
            x = numerator//(2*a)
            content = gcd(x, y)
            value = norm//(content*content)
            smallest = value if smallest is None else min(smallest, value)
    return smallest


def check_case(form, factors, complete=True):
    norm = prod(prime**exponent for prime, exponent, _ in factors)
    result = represent_arbitrary_factored_norm(form, factors)
    full, primitive = result['representation'], result['primitive_representation']
    expected = ellipse_status(form, norm) if complete else None
    if expected is not None:
        assert (full is not None) == expected[0], (form, factors, result, expected)
    if full is not None:
        content = result['integer_content']
        assert full.norm == norm and full.norm == content*content*primitive.norm
        assert gcd(*full.vector) == content and gcd(*primitive.vector) == 1
        assert tuple(content*entry for entry in primitive.vector) == full.vector
        assert prod(prime**abs(exponent) for prime, exponent, _ in result['prime_terms']) == primitive.norm
        a, middle, c = form
        for represented in (full, primitive):
            x, y = represented.vector
            assert a*x*x+middle*x*y+c*y*y == represented.norm
    else:
        assert primitive is None and result['prime_terms'] is None
    stats = result['stats']
    assert stats['left_capacity']*stats['right_capacity'] == stats['choice_product']
    assert stats['right_capacity']**2 <= stats['choice_product']*stats['maximum_radix']
    assert stats['right_products_examined'] <= stats['right_capacity']
    if not stats['odd_inert_exponent']:
        assert stats['left_products'] == stats['left_capacity']
    else:
        assert full is None and stats['group_multiplications'] == 0
    return result, expected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/arbitrary-norm-representation.json')
    args = parser.parse_args()
    partition_controls = 0
    for length in range(7):
        for radices in product(range(2, 6), repeat=length):
            indices, sizes = balanced_radix_partition(radices)
            assert sorted(indices[0]+indices[1]) == list(range(length))
            assert sizes[0]*sizes[1] == prod(radices)
            assert sizes[1] <= max(radices, default=1)*sizes[0]
            partition_controls += 1
    discriminants = (-3, -4, -7, -8, -15, -20, -23, -24, -31, -39,
                     -40, -47, -56, -84, -120, -163, -180, -220,
                     -324, -400, -41612)
    cases = represented = primitive_possible = imprimitive_outputs = 0
    inert_odd = noncanonical = bounded_controls = 0
    records = []
    for discriminant in discriminants:
        split, inert = base_for(discriminant, 43)
        base = sorted(split[:2]+inert[:2])
        factors_list = []
        for exponents in product(range(4), repeat=len(base)):
            if sum(exponents) > 6:
                continue
            factors = tuple((prime, exponent, root)
                            for (prime, root), exponent in zip(base, exponents) if exponent)
            if prod(prime**exponent for prime, exponent, _ in factors) <= 12000:
                factors_list.append(factors)
        forms = reduced_forms(discriminant)
        local_cases = local_successes = local_imprimitive_only = 0
        for form in forms:
            for factors in factors_list:
                result, expected = check_case(form, factors)
                local_cases += 1
                local_successes += expected[0]
                local_imprimitive_only += expected[0] and not expected[1]
                primitive_possible += expected[1]
                imprimitive_outputs += result['integer_content'] is not None and result['integer_content'] > 1
                inert_odd += result['stats']['odd_inert_exponent']
                if local_cases % 17 == 0:
                    norm = prod(prime**exponent for prime, exponent, _ in factors)
                    smallest = smallest_primitive_divisor(form, norm)
                    bounds = [norm] if smallest is None else [smallest]
                    if smallest is not None and smallest > 1:
                        bounds.append(smallest-1)
                    for bound in bounds:
                        bounded = represent_arbitrary_factored_norm(
                            form, factors, maximum_primitive_norm=bound)
                        assert (bounded['representation'] is not None) == (smallest is not None and smallest <= bound)
                        if bounded['representation'] is not None:
                            assert bounded['primitive_representation'].norm <= bound
                        bounded_controls += 1
                if local_cases % 37 == 0:
                    a, middle, c = form
                    shift = -3
                    changed = (a, middle+2*a*shift, c+middle*shift+a*shift*shift)
                    changed_result, changed_expected = check_case(changed, factors)
                    assert changed_expected == expected
                    assert (changed_result['representation'] is not None) == expected[0]
                    noncanonical += 1
        cases += local_cases
        represented += local_successes
        records.append({'discriminant': discriminant, 'class_number': len(forms),
                        'base': base, 'cases': local_cases,
                        'represented': local_successes,
                        'norms_requiring_imprimitive_representation': local_imprimitive_only})
    split_only, _ = check_case((1, 1, 6), ((2, 2, 1),))
    assert split_only['integer_content'] == 2
    assert split_only['primitive_representation'].norm == 1
    forced_inert, _ = check_case((1, 1, 6), ((5, 2, None),))
    assert forced_inert['integer_content'] == 5
    assert forced_inert['primitive_representation'].norm == 1

    competing_factors = ((7, 16, 8), (3, 2, None))
    first = represent_arbitrary_factored_norm((7, -6, 16), competing_factors)
    gated = represent_arbitrary_factored_norm(
        (7, -6, 16), competing_factors, maximum_primitive_norm=7**4)
    impossible_gate = represent_arbitrary_factored_norm(
        (7, -6, 16), competing_factors, maximum_primitive_norm=7**4-1)
    assert first['primitive_representation'].norm == 7**14
    assert gated['primitive_representation'].norm == 7**4
    assert gated['stats']['over_bound_matches'] > 0
    assert impossible_gate['representation'] is None
    left_collision = represent_arbitrary_factored_norm(
        (1, 0, 103), ((7, 16, 8), (13, 12, 2)), maximum_primitive_norm=1)
    assert left_collision['primitive_representation'].norm == 1
    assert left_collision['stats']['left_products'] == 13
    assert left_collision['stats']['over_bound_matches'] > 0

    discriminant = 1-4*1009**4
    split, inert = base_for(discriminant, 401)
    split = split[:12]
    factors = tuple((prime, 1+index % 3, root)
                    for index, (prime, root) in enumerate(split))
    signed = tuple((prime, exponent if index % 2 else -exponent, root)
                   for index, (prime, exponent, root) in enumerate(factors))
    known = represented_norm_from_word(discriminant, signed)
    factors += ((inert[0][0], 6, None),)
    large, _ = check_case(known.form, factors, complete=False)
    assert large['representation'] is not None
    assert large['integer_content'] % inert[0][0]**3 == 0
    malformed = [((1, 1, 6), ((2, 1, None),)),
                 ((1, 1, 6), ((3, 1, None),)),
                 ((1, 1, 6), ((5, 1, 1),)),
                 ((1, 1, 6), ((2, 1, 1), (2, 1, 1))),
                 ((2, 2, 12), ()), ((1, 3, 1), ())]
    for form, factors in malformed:
        try:
            represent_arbitrary_factored_norm(form, factors)
        except ValueError:
            pass
        else:
            raise AssertionError('Malformed representation input was accepted.')
    root = Path(__file__).resolve().parents[1]
    report = {
        'status': 'pass',
        'scope': 'Complete representation of a supplied factored norm; no smooth-norm selection bound.',
        'partition_controls': partition_controls,
        'complete_ellipse_comparisons': cases,
        'representations': represented,
        'primitive_representations_possible': primitive_possible,
        'norms_requiring_imprimitive_representation': represented-primitive_possible,
        'imprimitive_outputs': imprimitive_outputs,
        'odd_inert_exponent_decisions': inert_odd,
        'noncanonical_form_controls': noncanonical,
        'bounded_primitive_norm_ellipse_controls': bounded_controls,
        'malformed_input_controls': len(malformed),
        'split_square_only_control': {'full': asdict(split_only['representation']),
                                      'primitive': asdict(split_only['primitive_representation'])},
        'forced_inert_square_control': {'full': asdict(forced_inert['representation']),
                                        'primitive': asdict(forced_inert['primitive_representation'])},
        'competing_content_control': {
            'first_primitive_norm': first['primitive_representation'].norm,
            'bounded_primitive_norm': gated['primitive_representation'].norm,
            'bound': 7**4, 'one_smaller_bound_rejected': True,
            'stats': gated['stats'],
        },
        'classwise_left_minimum_control': {
            'primitive_representation': asdict(left_collision['primitive_representation']),
            'stats': left_collision['stats'],
        },
        'large_known_word_control': {'norm_bits': large['representation'].norm.bit_length(),
                                     'integer_content': large['integer_content'],
                                     'representation': asdict(large['representation']),
                                     'primitive_representation': asdict(large['primitive_representation']),
                                     'stats': large['stats']},
        'records': records,
        'engine_sha256': sha256((root/'tools/arbitrary_norm_representation.py').read_bytes()).hexdigest(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(args.output)
    if not output.is_absolute():
        output = root/output
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in (
        'status', 'complete_ellipse_comparisons', 'representations',
        'norms_requiring_imprimitive_representation', 'noncanonical_form_controls',
    )}))


if __name__ == '__main__':
    main()
