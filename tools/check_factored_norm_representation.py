#!/usr/bin/env python3
"""Compare orientation matching with exhaustive integer ellipse rows."""

from dataclasses import asdict
from itertools import product
import json
from math import gcd, isqrt, prod

from binary_quadratic_forms import reduced_forms
from factored_norm_representation import represent_factored_norm
from quadratic_ideal_hnf import represented_norm_from_word


def primes_at_most(bound):
    return [p for p in range(2, bound+1)
            if all(p % divisor for divisor in range(2, isqrt(p)+1))]


def split_factor_base(discriminant, bound):
    result = []
    for prime in primes_at_most(bound):
        if gcd(prime, discriminant) != 1:
            continue
        roots = [root for root in range(prime+1)
                 if (root*root-discriminant) % (4*prime) == 0]
        if roots:
            result.append((prime, roots[0]))
    return result


def complete_primitive_representation(form, norm):
    """Independent complete integer check via 4an=(2ax+by)^2+D y^2."""
    a, middle, c = form
    absolute_discriminant = 4*a*c-middle*middle
    limit = isqrt((4*a*norm)//absolute_discriminant)
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
            if gcd(x, y) == 1:
                assert a*x*x+middle*x*y+c*y*y == norm
                return x, y
    return None


def check_case(form, factors, use_complete_control=True):
    norm = prod(prime**exponent for prime, exponent, _ in factors)
    result = represent_factored_norm(form, factors)
    witness = result["representation"]
    if use_complete_control:
        independent = complete_primitive_representation(form, norm)
        assert (witness is not None) == (independent is not None), (form, factors, result)
    if witness is not None:
        assert witness.norm == norm and witness.removed_integer_content == 1
        x, y = witness.vector
        assert gcd(x, y) == 1
        assert form[0]*x*x+form[1]*x*y+form[2]*y*y == norm
        assert sum(abs(exponent) for _, exponent, _ in result["prime_terms"]) == sum(exponent for _, exponent, _ in factors)
    stats = result["stats"]
    count = len(factors)
    assert stats["left_products"] == 2**(count//2)
    assert stats["right_products_examined"] <= 2**((count+1)//2)
    # Binary prime powering, two product trees, and one comparison product
    # per right orientation. Sorting comparisons are not group operations.
    bound = (2*sum(exponent.bit_length() for _, exponent, _ in factors)
             +5*2**((count+1)//2))
    assert stats["group_multiplications"] <= bound
    return result


def main():
    discriminants = (-3, -4, -7, -8, -15, -20, -23, -24, -31, -39,
                     -40, -47, -56, -84, -120, -163, -180, -220,
                     -324, -400, -41612)
    records = []
    total_cases = total_successes = repeated_power_cases = noncanonical_cases = 0
    for discriminant in discriminants:
        base = split_factor_base(discriminant, 43)[:4]
        candidate_factors = []
        for exponents in product(range(3), repeat=len(base)):
            if sum(exponents) > 5:
                continue
            factors = tuple((prime, exponent, root)
                            for (prime, root), exponent in zip(base, exponents) if exponent)
            if prod(prime**exponent for prime, exponent, _ in factors) <= 10000:
                candidate_factors.append(factors)
        forms = reduced_forms(discriminant)
        cases = successes = 0
        for form in forms:
            for factors in candidate_factors:
                result = check_case(form, factors)
                cases += 1
                successes += result["representation"] is not None
                repeated_power_cases += any(exponent > 1 for _, exponent, _ in factors)
                if cases % 17 == 0:
                    a, middle, c = form
                    shift = -3
                    changed = (a, middle+2*a*shift, c+middle*shift+a*shift*shift)
                    changed_result = check_case(changed, factors)
                    assert (changed_result["representation"] is not None) == (result["representation"] is not None)
                    noncanonical_cases += 1
        total_cases += cases
        total_successes += successes
        records.append({"discriminant": discriminant, "class_number": len(forms),
                        "factor_base": base, "fixed_norm_cases": len(candidate_factors),
                        "form_norm_cases": cases, "primitive_representations": successes})

    repeated_control = check_case((2, -1, 3), ((2, 2, 1),))
    assert repeated_control["representation"] is not None
    assert check_case((2, -1, 3), ())["representation"] is None

    large_discriminant = 1-4*1009**4
    base = split_factor_base(large_discriminant, 401)[:18]
    assert len(base) == 18
    factors = tuple((prime, 1+index % 3, root)
                    for index, (prime, root) in enumerate(base))
    signed = tuple((prime, exponent if index % 2 else -exponent, root)
                   for index, (prime, exponent, root) in enumerate(factors))
    generated = represented_norm_from_word(large_discriminant, signed)
    large = check_case(generated.form, factors, use_complete_control=False)
    assert large["representation"] is not None
    print(json.dumps({
        "scope": "Exact representation of a supplied factored norm; the algorithm does not select the norm or enumerate classes",
        "discriminants": len(discriminants),
        "complete_ellipse_comparisons": total_cases,
        "primitive_representations": total_successes,
        "negative_fixed_norm_decisions": total_cases-total_successes,
        "cases_with_repeated_prime_powers": repeated_power_cases,
        "additional_noncanonical_form_controls": noncanonical_cases,
        "prime_square_control": {"form": (2, -1, 3), "norm": 4,
                                 "representation": asdict(repeated_control["representation"]),
                                 "squarefree_part_one_not_represented": True},
        "large_known_word_control": {"discriminant": large_discriminant,
                                     "norm_bits": large["representation"].norm.bit_length(),
                                     "representation": asdict(large["representation"]),
                                     "stats": large["stats"]},
        "records": records,
    }, indent=2))


if __name__ == "__main__":
    main()
