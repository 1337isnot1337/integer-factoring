#!/usr/bin/env python3
"""Finite exact checks of conic cancellation and a descent obstruction."""

from collections import Counter
from math import gcd
import json

from binary_quadratic_forms import reduced_forms
from conic_denominator_descent import parameter_point, recover_square_root, square_target_from_root
from quadratic_form_arithmetic import compose_primitive, inverse_primitive


def main():
    total = Counter()
    content_cases = Counter()
    examples = {}
    records = []
    parameters = [(r, s) for r in range(-5, 6) for s in range(0, 6) if gcd(r, s) == 1]
    for odd_n in range(3, 256, 2):
        counts = Counter()
        for root in reduced_forms(-4*odd_n):
            prepared = square_target_from_root(root, odd_n)
            counts["root_forms"] += 1
            if prepared["status"] == "factored":
                counts["initial_gcd_factors"] += 1
                divisor = prepared["factor"]
                assert 1 < divisor < odd_n and odd_n % divisor == 0
                continue
            counts["prepared_roots"] += 1
            for r, s in parameters:
                point = parameter_point(prepared["target"], r, s, prepared["root"])
                counts["parameter_points"] += 1
                if point["proper_factor"]:
                    counts["linear_gcd_factors"] += 1
                    content_cases["proper_factor"] += 1
                    continue
                case = ("N" if point["linear_gcd_with_N"] == odd_n else "1")
                if point["parity_factor"] == 2:
                    case = "2*"+case
                content_cases[case] += 1
                examples.setdefault(case, {"odd_n": odd_n, "root": root,
                                            "target": prepared["target"], "point": point})
                if "recovered_root" in point:
                    counts["unramified_root_recoveries"] += 1
                if not point["residual_norm_form_primitive"]:
                    counts["even_imprimitive_residuals"] += 1
        total.update(counts)
        records.append({"odd_n": odd_n, **dict(counts)})

    # S={3} generates this entire order-five class group.  The known root's
    # minimum represented denominator is 7, whereas its first S-smooth one is 9.
    root = (7, 6, 8)
    prepared = square_target_from_root(root, 47)
    assert prepared["target"] == (49, 20, 3)
    assert prepared["reduced_target"] == (3, -2, 16)
    point = parameter_point(prepared["target"], 2, 1, root)
    assert point["vector"] == (-1, 8) and point["denominator"] == 9
    generator = (3, 2, 16)
    powers, value = set(), (1, 0, 47)
    while value not in powers:
        powers.add(value)
        value = compose_primitive(value, generator)
    assert powers == set(reduced_forms(-188)) and len(powers) == 5
    transform = prepared["target_reduction_transform"]
    inv = ((transform[1][1], -transform[0][1]),
           (-transform[1][0], transform[0][0]))
    reduced_vector = (inv[0][0]*point["vector"][0]+inv[0][1]*point["vector"][1],
                      inv[1][0]*point["vector"][0]+inv[1][1]*point["vector"][1])
    assert reduced_vector == (5, 1)
    obstruction = {"odd_n": 47, "discriminant": -188, "root": root,
                   "split_support": [3], "class_number": 5,
                   "support_generates_full_group": True,
                   "initial_denominator": 7, "root_minimum_represented_value": 7,
                   "first_supported_denominator": 9, "prepared_target": prepared["target"],
                   "reduced_target": prepared["reduced_target"],
                   "supported_parameter": point, "reduced_supported_vector": reduced_vector}

    root = (7, -6, 44)
    prepared = square_target_from_root(root, 299)
    supported_point = parameter_point(prepared["target"], -3, 1, root)
    factoring_point = parameter_point(prepared["target"], -3, 2, root)
    assert supported_point["raw_point"] == (3, -8, 45) and supported_point["content"] == 1
    assert factoring_point["linear_gcd_with_N"] == 13 and factoring_point["denominator"] == 15
    other_root = recover_square_root(prepared["target"], factoring_point["vector"],
                                    factoring_point["denominator"])
    difference = compose_primitive(other_root, inverse_primitive(root))
    assert other_root == (15, -2, 20) and difference == (13, 0, 23)
    generators = [(3, 2, 100), (5, 2, 60)]
    classes = {(1, 0, 299)}
    while True:
        enlarged = classes | {compose_primitive(value, generator)
                              for value in classes for generator in generators}
        if enlarged == classes:
            break
        classes = enlarged
    assert classes == set(reduced_forms(-1196)) and len(classes) == 24
    composite_obstruction = {
        "odd_n": 299, "discriminant": -1196, "root": root,
        "split_support": [3, 5], "class_number": 24,
        "support_generates_full_group": True, "initial_denominator": 7,
        "root_minimum_represented_value": 7,
        "first_supported_denominator_in_same_root_class": 45,
        "prepared_target": prepared["target"], "reduced_target": prepared["reduced_target"],
        "supported_parameter": supported_point, "factoring_parameter": factoring_point,
        "recovered_other_root": other_root, "root_difference": difference,
    }
    assert all(content_cases[case] > 0 for case in ("1", "2*1", "N", "2*N", "proper_factor"))
    print(json.dumps({
        "scope": "Exact rational-conic denominator cancellation and its reduction to root-form smooth values; no smooth-value search algorithm is asserted",
        "odd_discriminants_checked": len(records), "parameter_stencil_size": len(parameters),
        "totals": dict(total), "content_cases": dict(content_cases),
        "content_case_examples": examples, "monotone_descent_obstruction": obstruction,
        "composite_no_factor_branch_obstruction": composite_obstruction,
        "records": records,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
