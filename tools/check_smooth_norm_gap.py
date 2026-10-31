#!/usr/bin/env python3
"""Exact controls for the radius-gap converter; no fast norm search is tested."""

import json
from math import gcd, isqrt, prod

from binary_quadratic_forms import compose_forms, inverse_form, reduced_forms
from bounded_smooth_norm import smooth_factors
from smooth_norm_gap import (
    balanced_divisor_with_gap,
    norm_to_pair_with_gap,
    pair_from_primitive_divisor,
    trimmed_lower_divisor,
)


def trial_factors(integer):
    result = []
    candidate = 2
    while candidate*candidate <= integer:
        while integer % candidate == 0:
            result.append(candidate)
            integer //= candidate
        candidate += 1
    if integer > 1:
        result.append(integer)
    return result


def exact_divisors(factors):
    divisors = {1}
    for prime in factors:
        divisors |= {value*prime for value in divisors}
    return sorted(divisors)


def check_divisors():
    cases = reordered_cases = 0
    worst = (1, 1)
    for integer in range(1, 4097):
        factors = trial_factors(integer)
        divisors = exact_divisors(factors)
        optimum = max(value for value in divisors if value*value <= integer)
        inner = integer // optimum
        for z in (1, 2, 4, 16):
            trace = trimmed_lower_divisor(factors, z, with_trace=True)
            value = trace["divisor"]
            assert integer % value == 0 and value <= optimum
            assert z*optimum <= (z+1)*value
            assert max(trace["retained_sizes"]) <= 1 + 2*trace["trim_denominator"]*integer.bit_length()
            outer = ((z+1)*inner + z-1)//z
            supplied = balanced_divisor_with_gap(factors, outer, z)
            assert supplied is not None and max(supplied, integer//supplied) <= outer
            if optimum*worst[1] > worst[0]*value:
                worst = (optimum, value)
            cases += 1
            if integer <= 256:
                reverse = trimmed_lower_divisor(list(reversed(factors)), z)
                assert integer % reverse == 0 and reverse <= optimum
                assert z*optimum <= (z+1)*reverse
                reordered_cases += 1
    # The approximate negative result must not be treated as exact-radius failure.
    factors = trial_factors(220)
    assert 220 == 11*20
    assert trimmed_lower_divisor(factors, 1) == 10
    assert balanced_divisor_with_gap(factors, 20, 1) is None
    assert balanced_divisor_with_gap(factors, 40, 1) == 10
    return {"integers_through": 4096, "approximation_cases": cases,
            "reversed_item_controls": reordered_cases,
            "worst_observed_optimum_to_output_ratio": list(worst),
            "exact_radius_negative_control": {
                "integer": 220, "prime_factors": factors,
                "valid_exact_pair": [11, 20], "trimmed_pair": [10, 22],
                "exact_radius": 20, "gap_radius": 40}}


def check_forms():
    records = []
    total = scaled = shared_prime_controls = 0
    discriminants = (-348, -420, -999, -1679, -2639, -9679,
                    -10007, -41612, -43228, -66548)
    for discriminant in discriminants:
        outer = isqrt(-discriminant-1)//2
        forms = reduced_forms(discriminant)
        for z in (1, 2, 4):
            inner = z*outer//(z+1)
            accepted = target_classes = 0
            assert (z+1)*inner <= z*outer
            for form in forms:
                a, b, c = form
                y_bound = isqrt((4*a*inner*inner)//(-discriminant))
                found = False
                for y in range(-y_bound, y_bound+1):
                    root = isqrt(4*a*inner*inner + discriminant*y*y)
                    for x in range(-((root+b*y)//(2*a)), (root-b*y)//(2*a)+1):
                        if gcd(x, y) != 1:
                            continue
                        norm = a*x*x+b*x*y+c*y*y
                        if norm <= 0 or gcd(norm, discriminant) != 1:
                            continue
                        factors = smooth_factors(norm, 7)
                        if factors is None:
                            continue
                        possible = [value for value in exact_divisors(factors)
                                    if value <= inner and norm//value <= inner]
                        if not possible:
                            continue
                        pair = norm_to_pair_with_gap(form, x, y, 7, outer, z)
                        assert pair is not None
                        assert all(part[0] <= outer for part in pair)
                        assert pair[0][0]*pair[1][0] == norm
                        assert compose_forms(pair[0], inverse_form(pair[1])) == form
                        assert norm_to_pair_with_gap(form, 3*x, 3*y, 7, outer, z) == pair
                        scaled += 1
                        # Independently test the direct construction at every
                        # exact inner split, including repeated prime factors.
                        for divisor in possible:
                            direct = pair_from_primitive_divisor(form, x, y, divisor, outer)
                            assert compose_forms(direct[0], inverse_form(direct[1])) == form
                            if gcd(divisor, norm//divisor) > 1:
                                shared_prime_controls += 1
                        accepted += 1
                        total += 1
                        found = True
                target_classes += found
            records.append({"discriminant": discriminant, "inverse_slack": z,
                            "inner_radius": inner, "outer_radius": outer,
                            "promised_representations": accepted,
                            "target_classes": target_classes})
    assert total > 0 and shared_prime_controls > 0
    # The exact-radius rejection above also occurs for an actual primitive form.
    form = (11, 1, 220)
    assert norm_to_pair_with_gap(form, 0, 1, 11, 20, 1) is None
    pair = norm_to_pair_with_gap(form, 0, 1, 11, 40, 1)
    assert [part[0] for part in pair] == [10, 22]
    assert compose_forms(pair[0], inverse_form(pair[1])) == form
    return {"total_promised_representations": total,
            "scaled_nonprimitive_controls": scaled,
            "shared_prime_direct_controls": shared_prime_controls,
            "records": records,
            "actual_form_gap_control": {"form": list(form), "norm": 220,
                                        "exact_radius": 20, "gap_radius": 40,
                                        "pair": [list(part) for part in pair]}}


def check_large_square():
    primes = []
    candidate = 2
    while len(primes) < 32:
        if all(candidate % prime for prime in primes if prime*prime <= candidate):
            primes.append(candidate)
        candidate += 1
    # The complete divisor set has 3^32 members.  Its optimum is known exactly,
    # since the integer is a square; it is never enumerated by this control.
    factors = [prime for prime in primes for _ in range(2)]
    square_root = prod(primes)
    trace = trimmed_lower_divisor(factors, 8, with_trace=True)
    output = trace["divisor"]
    assert trace["integer"] == square_root*square_root
    assert output <= square_root and 8*square_root <= 9*output
    return {"distinct_primes": len(primes), "prime_occurrences": len(factors),
            "input_bits": trace["integer"].bit_length(),
            "complete_divisor_count": 3**len(primes),
            "inverse_slack": 8, "maximum_retained_list": max(trace["retained_sizes"]),
            "total_candidates_processed": trace["candidate_count"],
            "exact_optimum": str(square_root), "returned_divisor": str(output),
            "complete_divisor_set_enumerated": False}


def main():
    print(json.dumps({
        "scope": "Polynomial supplied-norm conversion with a certified inner/outer radius gap; no fast witness search or improved factoring score",
        "divisor_controls": check_divisors(),
        "form_controls": check_forms(),
        "large_square_control": check_large_square(),
        "fast_witness_search_supplied": False,
        "exact_radius_decision_claimed": False,
    }, indent=2))


if __name__ == "__main__":
    main()
