#!/usr/bin/env python3
"""Independent finite controls for the norm sieve and arithmetic pipeline.

The brute-force fixture uses the two coordinate bounds separately and
factors each small norm by trial division. It does not reuse the sieve's
row bounds, root progressions, or smoothness decisions.
"""

from dataclasses import asdict
from math import gcd, isqrt
import json

from arithmetic_short_norm_factor import factor_with_short_norms
from arithmetic_two_group import power
from binary_quadratic_forms import reduced_forms
from quadratic_ideal_hnf import compose_forms_hnf, reduce_with_transform
from smooth_norm_sieve import SmoothNormSieve


def independent_multiplicity(integer, smooth_bound):
    remaining, multiplicity = integer, 0
    divisor = 2
    while divisor * divisor <= remaining:
        while remaining % divisor == 0:
            if divisor > smooth_bound:
                return None
            remaining //= divisor
            multiplicity += 1
        divisor += 1
    if remaining > 1:
        if remaining > smooth_bound:
            return None
        multiplicity += 1
    return multiplicity


def brute_witnesses(form, maximum_norm, smooth_bound):
    a, b, c = form
    discriminant = b * b - 4 * a * c
    x_bound = isqrt((4 * c * maximum_norm) // -discriminant)
    y_bound = isqrt((4 * a * maximum_norm) // -discriminant)
    witnesses = []
    for x in range(-x_bound, x_bound + 1):
        for y in range(-y_bound, y_bound + 1):
            if gcd(x, y) != 1:
                continue
            norm = a * x * x + b * x * y + c * y * y
            if not 0 < norm <= maximum_norm or gcd(norm, discriminant) != 1:
                continue
            multiplicity = independent_multiplicity(norm, smooth_bound)
            if multiplicity is not None:
                witnesses.append((x, y, norm, multiplicity))
    return witnesses


def main():
    discriminants = (-23, -31, -47, -71, -84, -144, -999, -1007)
    smooth_bounds = (2, 3, 5, 7, 11)
    caps = (20, 100, 500)
    budgets = (0, 1, 2, 4, 8)
    controls = successes = exhausted = root_controls = all_root_rows = 0
    for discriminant in discriminants:
        forms = reduced_forms(discriminant)
        for smooth_bound in smooth_bounds:
            sieve = SmoothNormSieve(discriminant, smooth_bound)
            for form in forms:
                a, b, c = form
                for prime, root in sieve.prime_roots:
                    for y in range(2 * prime):
                        residues = sieve._row_roots(form, y, prime, root)
                        actual = set(range(prime)) if residues is None else set(residues)
                        expected = {x for x in range(prime)
                                    if (a*x*x + b*x*y + c*y*y) % prime == 0}
                        assert actual == expected
                        root_controls += 1
                        all_root_rows += residues is None
                for maximum_norm in caps:
                    witnesses = brute_witnesses(form, maximum_norm, smooth_bound)
                    for budget in budgets:
                        found = sieve.search(form, maximum_norm, budget)
                        expected_exists = any(row[3] <= budget for row in witnesses)
                        assert (found.vector is not None) == expected_exists
                        assert found.exhausted != expected_exists
                        if found.vector is not None:
                            x, y = found.vector
                            assert (x, y, found.norm, found.multiplicity) in witnesses
                            assert y >= 0 and (y != 0 or x == 1)
                            assert found.multiplicity <= budget
                            successes += 1
                        else:
                            exhausted += 1
                        controls += 1

    # This is a concrete infeasibility control, not just a search failure.
    # y != 0 forces Q(x,y) >= N/81 > 31^3; at y=0 primitiveness
    # forces x=+/-1, whose norm 81 has four prime factors with multiplicity.
    obstructed_integer = 10000001
    obstructed_delta = -4 * obstructed_integer
    prime_form = reduce_with_transform((3, 2, (4-obstructed_delta)//12))[0]
    obstructed_form = power(prime_form, 4, (1, 0, obstructed_integer), compose_forms_hnf)
    assert obstructed_form == (81, 8, 123457)
    assert obstructed_integer > 81 * 31**3
    obstruction = SmoothNormSieve(obstructed_delta, 31).search(obstructed_form, 31**3, 3)
    assert obstruction.exhausted and obstruction.coordinate_points_visited == 1

    # The factors exceed the chosen prime cutoff. Each successful output
    # therefore comes from searched relations and ambiguous extraction.
    integers = (1022117, 100160063, 10002200057)
    pipeline = []
    for integer in integers:
        record = factor_with_short_norms(integer, smooth_bound=127,
                                        maximum_expansions=8, keep_relations=True)
        assert record["status"] == "factor"
        assert record["method"].startswith("searched short norms")
        assert integer % record["factor"] == 0
        assert not record["complete_class_group_enumerated"]
        assert not record["generator_orders_searched"]
        assert record["ambiguous_extraction"]["stats"]["peak_retained_cyclic_directions"] <= 1
        pipeline.append(record)
    prime_control = factor_with_short_norms(1009, smooth_bound=31, maximum_expansions=4)
    assert prime_control["status"] == "unresolved"

    print(json.dumps({
        "scope": "Complete finite norm-search controls and bounded exact factoring runs; no asymptotic search or generation claim",
        "discriminants": discriminants,
        "smooth_bounds": smooth_bounds,
        "norm_caps": caps,
        "multiplicity_budgets": budgets,
        "brute_force_existence_controls": controls,
        "successful_search_controls": successes,
        "exhausted_search_controls": exhausted,
        "modular_root_controls": root_controls,
        "all_root_progression_controls": all_root_rows,
        "multiplicity_obstruction": {
            "integer": obstructed_integer,
            "form": obstructed_form,
            "search": asdict(obstruction),
            "nonzero_second_coordinate_norm_lower_bound": [obstructed_integer, 81],
        },
        "factoring_pipeline_records": pipeline,
        "prime_input_control": prime_control,
    }, indent=2))


if __name__ == "__main__":
    main()
