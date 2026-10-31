#!/usr/bin/env python3
"""Independent coordinate enumeration for the exact square-box solver."""

from itertools import product
from math import isqrt, prod
import json

from binary_quadratic_forms import reduced_forms
from check_smooth_norm_relation import factor_base
from factored_square_box_search import SquareBoxSearch, balanced_rectangles
from quadratic_form_arithmetic import compose_primitive, power_primitive


def represented_norms(form, bound):
    a, b, c = form
    discriminant = b*b-4*a*c
    values = set()
    for y in range(-isqrt(4*a*bound//(-discriminant)), isqrt(4*a*bound//(-discriminant))+1):
        root = isqrt(4*a*bound+discriminant*y*y)
        for x in range(-((root+b*y)//(2*a)), (root-b*y)//(2*a)+1):
            norm = a*x*x+b*x*y+c*y*y
            if 1 <= norm <= bound:
                values.add(norm)
    return values


def small_controls():
    cases = positive = negative = nonprimitive = remainder_hits = transformed = 0
    records = []
    for absolute in range(3, 193):
        discriminant = -absolute
        if discriminant % 4 not in (0, 1):
            continue
        forms = reduced_forms(discriminant)
        base = factor_base(discriminant, 7)
        primes = [prime for prime, _ in base]
        roots = {prime: next(value for value in range(prime+1)
                             if (value*value-discriminant) % (4*prime) == 0)
                 for prime in primes}
        direct = {form: represented_norms(form, 4096) for form in forms}
        local_cases = local_positive = 0
        for caps in product(range(4), repeat=len(primes)):
            root_norm = prod(prime**cap for prime, cap in zip(primes, caps))
            if root_norm > 64:
                continue
            search = SquareBoxSearch(discriminant, primes, roots, caps)
            for form in forms:
                answer = search.solve(form)
                expected = root_norm*root_norm in direct[form]
                assert answer["represented"] == expected, (form, primes, caps, answer)
                cases += 1
                local_cases += 1
                if not expected:
                    negative += 1
                    continue
                positive += 1
                local_positive += 1
                witness = answer["representation"]
                assert witness["root_norm"] == root_norm
                x, y = witness["vector"]
                assert form[0]*x*x+form[1]*x*y+form[2]*y*y == root_norm*root_norm
                nonprimitive += witness["content"] > 1
                remainder_hits += answer["rectangle_index"] > 0
                if transformed < 100:
                    a, b, c = form
                    changed = (a+17*b+17*17*c, b+34*c, c)
                    replay = search.solve(changed)
                    assert replay["represented"]
                    xx, yy = replay["representation"]["vector"]
                    assert changed[0]*xx*xx+changed[1]*xx*yy+changed[2]*yy*yy == root_norm*root_norm
                    transformed += 1
        records.append({"discriminant": discriminant, "classes": len(forms),
                        "factor_base_primes": primes, "decisions": local_cases,
                        "positive_decisions": local_positive})
    return {"exact_decisions": cases, "positive_decisions": positive,
            "negative_decisions": negative, "nonprimitive_returned_witnesses": nonprimitive,
            "remainder_rectangle_hits": remainder_hits,
            "nonreduced_target_controls": transformed, "records": records}


def large_controls():
    # One coordinate with a large cap: splitting the interval is essential
    # to avoid putting its entire numeric length on one MITM side.
    discriminant = 1-8*((1 << 256)+(1 << 128)+12345)
    cap = 10**6
    search = SquareBoxSearch(discriminant, [2], {2: 1}, [cap])
    supplied_exponent = 123456
    target = power_primitive(search.generators[0], 2*supplied_exponent)
    result = search.solve(target, with_representation=False)
    assert result["represented"] and abs(result["half_word"][0]) <= cap
    assert power_primitive(search.generators[0], 2*result["half_word"][0]) == target
    first = {"discriminant_bits": (-discriminant).bit_length(), "caps": [cap],
             "statistics": search.statistics(), "recovered_half_word": result["half_word"],
             "query_group_multiplications": result["query_group_multiplications"],
             "coordinates_requested": False}

    primes = []
    candidate = 2
    while len(primes) < 16:
        if all(candidate % prime for prime in primes if prime*prime <= candidate):
            primes.append(candidate)
        candidate += 1
    discriminant = 1-4*prod(primes)*(1 << 64)
    search = SquareBoxSearch(discriminant, primes, {prime: 1 for prime in primes}, [1]*16)
    supplied_word = [(index % 3)-1 for index in range(16)]
    target = search.identity
    for generator, exponent in zip(search.generators, supplied_word):
        target = compose_primitive(target, power_primitive(generator, 2*exponent))
    result = search.solve(target)
    assert result["represented"]
    witness = result["representation"]
    x, y = witness["vector"]
    assert target[0]*x*x+target[1]*x*y+target[2]*y*y == witness["fixed_norm"]
    second = {"discriminant_bits": (-discriminant).bit_length(), "prime_support_size": 16,
              "caps": [1]*16, "statistics": search.statistics(),
              "recovered_half_word": result["half_word"],
              "query_group_multiplications": result["query_group_multiplications"],
              "fixed_norm_bits": witness["fixed_norm"].bit_length(),
              "nonprimitive_content": witness["content"], "vector": witness["vector"]}
    return [first, second]


def main():
    # These highly unbalanced boxes exercise both rectangles and caps zero.
    shapes = [(1,), (2000001,), (3, 101), (1, 3, 101, 1), (9,)*9, (7, 11, 13)]
    rectangle_controls = [{"radices": radices, "rectangles": balanced_rectangles(radices)}
                          for radices in shapes]
    print(json.dumps({
        "scope": "Exact deterministic variable-cap smooth-square MITM search, with square-root exponent-box cost fully charged",
        "small_coordinate_controls": small_controls(),
        "large_controls": large_controls(), "balanced_rectangle_controls": rectangle_controls,
        "polynomial_representation_oracle_claimed": False,
        "improved_complete_factoring_bound_claimed": False,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
