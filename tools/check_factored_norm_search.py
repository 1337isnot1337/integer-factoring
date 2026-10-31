#!/usr/bin/env python3
"""Independent exact controls for ideal composition and candidate-norm search."""

from itertools import combinations
from math import gcd, isqrt, prod
import json

from binary_quadratic_forms import compose_forms, reduced_forms, reduce_form
from bounded_smooth_norm import smooth_factors
from factored_norm_search import (
    combine_roots, lift_prime_root, represent_smooth_norm, represent_with_prime_powers,
)
from quadratic_form_arithmetic import (
    compose_primitive, composition_data, inverse_primitive, reduce_with_transform,
)


def transformed_form(form, matrix):
    a, b, c = form
    (u, v), (w, z) = matrix
    return (a*u*u+b*u*w+c*w*w,
            2*a*u*v+b*(u*z+v*w)+2*c*w*z,
            a*v*v+b*v*z+c*z*z)


def check_arithmetic():
    pair_cases = transform_cases = common_content_cases = 0
    for absolute in range(3, 1000):
        discriminant = -absolute
        if discriminant % 4 not in (0, 1):
            continue
        forms = reduced_forms(discriminant)
        for form in forms:
            for shift in (-31, -1, 0, 2, 29):
                a, b, c = form
                changed = (a+b*shift+c*shift*shift, b+2*c*shift, c)
                reduced, matrix = reduce_with_transform(changed)
                assert reduced == form == reduce_form(changed)
                assert transformed_form(changed, matrix) == reduced
                assert matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0] == 1
                transform_cases += 1
            for other in forms:
                data = composition_data(form, other)
                assert compose_primitive(form, other) == compose_forms(form, other)
                vectors = data["ideal_generators"]
                minor_gcd = 0
                for left, right in combinations(vectors, 2):
                    minor_gcd = gcd(minor_gcd, left[0]*right[1]-left[1]*right[0])
                assert minor_gcd == form[0]*other[0]
                content = gcd(*(entry for vector in vectors for entry in vector))
                assert content == data["content"]
                index = data["primitive_ideal_index"]
                horizontal, one = data["vertical_one_after_content"]
                assert one == 1 and index*content*content == minor_gcd
                for x, y in vectors:
                    assert x % content == y % content == 0
                    assert (x//content-horizontal*(y//content)) % index == 0
                common_content_cases += content > 1
                pair_cases += 1
    return {"discriminant_absolute_value_below": 1000,
            "independent_composition_pairs": pair_cases,
            "composition_pairs_with_nontrivial_content": common_content_cases,
            "exact_reduction_transform_controls": transform_cases}


def brute_primitive_norms(form, bound):
    a, b, c = form
    discriminant = b*b-4*a*c
    output = {}
    y_bound = isqrt((4*a*bound)//(-discriminant))
    for y in range(-y_bound, y_bound+1):
        root = isqrt(4*a*bound+discriminant*y*y)
        for x in range(-((root+b*y)//(2*a)), (root-b*y)//(2*a)+1):
            if gcd(x, y) == 1:
                norm = a*x*x+b*x*y+c*y*y
                if 1 <= norm <= bound:
                    output.setdefault(norm, (x, y))
    return output


def check_norm_decisions():
    cases = positive = negative = inert = transformed = 0
    records = []
    discriminants = [-value for value in range(3, 257) if (-value) % 4 in (0, 1)]
    discriminants += [-999, -1679, -4799, -10007]
    for discriminant in discriminants:
        forms = reduced_forms(discriminant)
        local_cases = local_positive = 0
        for form in forms:
            direct = brute_primitive_norms(form, 128)
            for integer in range(1, 129):
                if gcd(integer, discriminant) != 1 or smooth_factors(integer, 7) is None:
                    continue
                result = represent_smooth_norm(form, integer, 7)
                exists = integer in direct
                assert (result["status"] == "represented") == exists, (form, integer, result)
                assert result["status"] in ("represented", "no_primitive_representation")
                if exists:
                    x, y = result["primitive_vector"]
                    a, b, c = form
                    assert gcd(x, y) == 1 and a*x*x+b*x*y+c*y*y == integer
                    positive += 1
                    local_positive += 1
                    if transformed < 500:
                        shift = 17
                        changed = (a+b*shift+c*shift*shift, b+2*c*shift, c)
                        replay = represent_smooth_norm(changed, integer, 7)
                        assert replay["status"] == "represented"
                        xx, yy = replay["primitive_vector"]
                        assert changed[0]*xx*xx+changed[1]*xx*yy+changed[2]*yy*yy == integer
                        transformed += 1
                else:
                    negative += 1
                    inert += result["reason"] == "inert_prime"
                if "peak_half_states" in result:
                    assert max(result["peak_half_states"]) <= 1 << ((result["omega"]+1)//2)
                cases += 1
                local_cases += 1
        if discriminant in (-4, -23, -999, -1679, -4799, -10007):
            records.append({"discriminant": discriminant, "classes": len(forms),
                            "candidate_decisions": local_cases, "positive": local_positive})
    # An inert even power has nonprimitive representations but no primitive one.
    inert_square = represent_smooth_norm((1, 0, 1), 9, 3)
    assert inert_square["status"] == "no_primitive_representation"
    assert inert_square["reason"] == "inert_prime"
    assert represent_smooth_norm((1, 0, 1), 13, 7)["status"] == "not_smooth"
    return {"discriminants_checked": len(discriminants), "norm_cutoff": 128,
            "smooth_bound": 7, "exact_candidate_decisions": cases,
            "positive_decisions": positive, "negative_decisions": negative,
            "inert_prime_rejections": inert, "nonreduced_target_controls": transformed,
            "selected_records": records,
            "inert_square_control": {"form": [1, 0, 1], "norm": 9,
                                     "nonprimitive_vector": [3, 0],
                                     "primitive_representation_exists": False}}


def primes_through_count(count):
    primes = []
    candidate = 2
    while len(primes) < count:
        if all(candidate % prime for prime in primes if prime*prime <= candidate):
            primes.append(candidate)
        candidate += 1
    return primes


def check_large_controls():
    records = []
    specs = [([(prime, 1) for prime in primes_through_count(24)], 512),
             (list(zip(primes_through_count(8), (48, 32, 24, 16, 12, 8, 6, 4))), 2048)]
    for powers, bits in specs:
        prime_product = prod(prime for prime, _ in powers)
        discriminant = 1-4*prime_product*(1 << bits)
        roots = [lift_prime_root(discriminant, prime, exponent, 1)
                 for prime, exponent in powers]
        mask = sum(1 << index for index in range(len(powers)) if index % 3 == 1)
        norm = prod(prime**exponent for prime, exponent in powers)
        root = combine_roots(discriminant, powers, roots, mask)
        target, _ = reduce_with_transform((norm, root, (root*root-discriminant)//(4*norm)))
        result = represent_with_prime_powers(target, powers, {prime: 1 for prime, _ in powers})
        assert result["status"] == "represented"
        x, y = result["primitive_vector"]
        assert gcd(x, y) == 1 and target[0]*x*x+target[1]*x*y+target[2]*y*y == norm
        assert compose_primitive(target, inverse_primitive(target))[0] == 1
        records.append({"discriminant_bits": (-discriminant).bit_length(),
                        "norm_bits": norm.bit_length(), "omega": len(powers),
                        "total_prime_occurrences": sum(exponent for _, exponent in powers),
                        "all_orientation_count": 1 << len(powers),
                        "half_state_counts": result["half_state_counts"],
                        "composition_calls": result["composition_calls"],
                        "primitive_vector": [str(x), str(y)],
                        "direct_norm_identity_checked": True})
    return records


def main():
    print(json.dumps({
        "scope": "Polynomial ideal arithmetic and exact supplied-candidate norm representation; candidate selection is not provided",
        "arithmetic_controls": check_arithmetic(),
        "candidate_norm_controls": check_norm_decisions(),
        "large_integer_controls": check_large_controls(),
        "candidate_norm_selection_supplied": False,
        "improved_complete_factoring_bound_claimed": False,
    }, indent=2))


if __name__ == "__main__":
    main()
