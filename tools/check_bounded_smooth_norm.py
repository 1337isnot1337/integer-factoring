#!/usr/bin/env python3
"""Exact supplied-norm to pair conversion; witness search here is exhaustive."""
from math import gcd, isqrt
import json

from binary_quadratic_forms import compose_forms, inverse_form, reduced_forms, smooth_prime_word
from bounded_smooth_norm import norm_to_pair, smooth_factors


def main():
    records, representations, scaled_controls, outside_smooth_set, unbalanced = [], 0, 0, 0, 0
    for discriminant in (-348, -420, -4620, -5460, -10007, -41612, -43228, -66548):
        bound, radius = 7, isqrt(-discriminant-1)//2
        maximum_norm = radius*radius
        forms = reduced_forms(discriminant)
        targets, non_smooth_targets, reached = 0, 0, set()
        for form in forms:
            a, b, c = form
            y_bound = isqrt((4*a*maximum_norm)//(-discriminant))
            found = False
            for y in range(-y_bound, y_bound+1):
                root = isqrt(4*a*maximum_norm+discriminant*y*y)
                for x in range(-((root+b*y)//(2*a)), (root-b*y)//(2*a)+1):
                    norm = a*x*x+b*x*y+c*y*y
                    if not norm or gcd(norm, discriminant) != 1 or smooth_factors(norm, bound) is None:
                        continue
                    pair = norm_to_pair(form, x, y, bound, radius)
                    if pair is None:
                        unbalanced += 1
                        continue
                    assert all(value[0] <= radius and smooth_prime_word(value, bound) is not None for value in pair)
                    representations += 1
                    if gcd(2, discriminant) == 1 and 4*norm <= maximum_norm:
                        assert norm_to_pair(form, 2*x, 2*y, bound, radius) == pair
                        scaled_controls += 1
                    found = True
            if found:
                targets += 1
                reached.add(form)
                if smooth_prime_word(form, bound) is None:
                    non_smooth_targets += 1
        smooth_set = [form for form in forms if form[0] <= radius and smooth_prime_word(form, bound) is not None]
        quotients = {compose_forms(left, inverse_form(right)) for left in smooth_set for right in smooth_set}
        assert reached == quotients
        outside_smooth_set += non_smooth_targets
        records.append({'discriminant': discriminant, 'class_number': len(forms),
                        'radius': radius, 'maximum_supplied_norm': maximum_norm,
                        'classes_with_a_bounded_smooth_norm': targets,
                        'exact_smooth_quotient_set_size': len(quotients),
                        'such_classes_outside_original_smooth_set': non_smooth_targets})
    assert outside_smooth_set > 0 and scaled_controls > 0
    # All factors are small, but greedy balancing misses the 14*15 split.
    pair = norm_to_pair((2, -1, 210), 0, 1, 7, 15)
    assert sorted(form[0] for form in pair) == [14, 15]
    # A smooth represented norm alone need not admit a balanced split.
    assert norm_to_pair((2, -1, 125), 0, 1, 7, 12) is None
    print(json.dumps({
        'scope': 'Exact conversion from a supplied bounded smooth norm to a smooth-form quotient; finite witness search is exhaustive',
        'representations_checked': representations,
        'nonprimitive_scaled_controls': scaled_controls,
        'classes_reached_outside_original_smooth_set': outside_smooth_set,
        'unbalanced_smooth_representations_rejected': unbalanced,
        'exact_divisor_search_control': {'discriminant': -1679, 'norm': 210, 'pair_norms': [form[0] for form in pair]},
        'unbalanced_control': {'discriminant': -999, 'norm': 125, 'radius': 12},
        'smooth_bound': 7, 'records': records,
        'fast_witness_search_supplied': False,
    }, indent=2))


if __name__ == '__main__':
    main()
