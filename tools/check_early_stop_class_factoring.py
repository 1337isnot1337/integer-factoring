#!/usr/bin/env python3
"""Finite exact controls for early-stop class-group factor extraction."""

from collections import Counter
from math import gcd
import json

from binary_quadratic_forms import reduced_forms
from early_stop_class_factoring import (
    ambiguous_factorization, factor_from_annihilated_generators,
)
from quadratic_form_arithmetic import compose_primitive, power_primitive


def factorization(value):
    result, prime = {}, 2
    while prime*prime <= value:
        while value % prime == 0:
            result[prime] = result.get(prime, 0)+1
            value //= prime
        prime += 1
    if value > 1:
        result[value] = result.get(value, 0)+1
    return result


def check(odd_n):
    forms = reduced_forms(-4*odd_n)
    identity = (1, 0, odd_n)
    group_order = len(forms)
    known_factors = factorization(odd_n)
    expected_factor = len(known_factors) > 1
    ambiguous = [form for form in forms if compose_primitive(form, form) == identity]
    unhelpful = []
    shapes = Counter()
    for form in ambiguous:
        entry = ambiguous_factorization(form, odd_n)
        shapes[entry["shape"]] += 1
        if entry["proper_factor"] is None:
            unhelpful.append(form)
    assert identity in unhelpful
    assert len(unhelpful) == (2 if odd_n % 4 == 1 else 1)
    assert set(unhelpful) <= {identity, (2, 2, (odd_n+1)//2)}
    if expected_factor:
        assert len(ambiguous) > len(unhelpful)
    runs = []
    for label, generators in (("canonical", forms), ("reverse", list(reversed(forms)))):
        # An intentionally oversaturated annihilator: the algorithm must not
        # require an exact order, its odd factorization, or a minimal 2-part.
        annihilator = group_order*15*(1 << (odd_n % 7))
        answer = factor_from_annihilated_generators(odd_n, generators, annihilator, with_trace=True)
        assert (answer["factor"] is not None) == expected_factor
        stats = answer["statistics"]
        assert stats["maximum_retained_basis_rank"] <= 1
        if answer["factor"]:
            factor = answer["factor"]
            assert 1 < factor < odd_n and odd_n % factor == 0
            assert gcd(factor, odd_n//factor) == 1
        for entry in answer["trace"]:
            top = tuple(entry["top"])
            assert top in ambiguous and top != identity
            assert power_primitive(tuple(entry["candidate"]),
                                   1 << (entry["candidate_order_log"]-1)) == top
        runs.append({"ordering": label, "factor": answer["factor"], "statistics": stats})
    return {"odd_n": odd_n, "class_number": group_order,
            "distinct_prime_factors_in_checker": len(known_factors),
            "non_squarefree": any(exponent > 1 for exponent in known_factors.values()),
            "ambiguous_classes": len(ambiguous), "unhelpful_classes": unhelpful,
            "ambiguous_shapes": dict(shapes), "runs": runs}


def main():
    records = [check(odd_n) for odd_n in range(3, 1024, 2)]
    stats = [run["statistics"] for record in records for run in record["runs"]]
    mixed_generators = [(2, 2, 81), (3, 2, 54), (5, 4, 33)]
    states = {(1, 0, 161)}
    while True:
        enlarged = states | {compose_primitive(value, generator)
                             for value in states for generator in mixed_generators}
        if enlarged == states:
            break
        states = enlarged
    assert states == set(reduced_forms(-644))
    mixed_control = factor_from_annihilated_generators(161, mixed_generators, 240,
                                                      with_trace=True)
    assert mixed_control["factor"] == 7
    assert mixed_control["statistics"]["order_swaps"] == 1
    assert mixed_control["statistics"]["dependent_top_cancellations"] == 3
    # A proper subgroup may have only the exceptional ambiguous direction.
    # This must not be mistaken for a primality or prime-power certificate.
    subgroup_control = factor_from_annihilated_generators(21, [(2, 2, 11)], 2)
    assert subgroup_control["factor"] is None
    assert subgroup_control["status"] == "no_factor_in_supplied_subgroup"
    invalid_rejected = False
    try:
        factor_from_annihilated_generators(21, [(2, 2, 11)], 1)
    except ValueError as error:
        assert "does not annihilate" in str(error)
        invalid_rejected = True
    assert invalid_rejected
    assert sum(item["dependent_top_cancellations"] for item in stats) > 0
    assert sum(item["order_swaps"] for item in stats) > 0
    print(json.dumps({
        "scope": "Polynomial early-stop factor extraction from supplied complete class-group generators and an annihilator; finite checker obtains those inputs by full enumeration",
        "odd_integers_checked": len(records), "complete_generator_runs": len(stats),
        "non_squarefree_inputs": sum(record["non_squarefree"] for record in records),
        "composite_non_prime_power_inputs": sum(record["distinct_prime_factors_in_checker"] > 1 for record in records),
        "ambiguous_forms_checked": sum(record["ambiguous_classes"] for record in records),
        "dependent_top_cancellations": sum(item["dependent_top_cancellations"] for item in stats),
        "order_swaps": sum(item["order_swaps"] for item in stats),
        "maximum_retained_basis_rank": max(item["maximum_retained_basis_rank"] for item in stats),
        "mixed_top_cancellation_control": mixed_control,
        "proper_subgroup_negative_control": subgroup_control,
        "invalid_annihilator_rejected": invalid_rejected,
        "records": records,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
