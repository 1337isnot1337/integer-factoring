#!/usr/bin/env python3
"""Exact norm identities, independent conic points, and content split controls."""

from math import gcd, isqrt
from pathlib import Path
import json

from binary_quadratic_forms import ambiguous_factor, compose_forms, reduced_forms
from bounded_smooth_norm import smooth_factors
from class_square_conic import (
    content_factor, cross, evaluate, modular_unit_point, parameters_from_point, parametrized_point,
    representation_from_denominator, square_parametrization,
)
from quadratic_form_arithmetic import power_primitive


def polynomial_product(left, right):
    result = [0]*5
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return result


def check_small_forms():
    forms_checked = identities = forward = inverse = normalized = nontrivial = 0
    discriminants = 0
    all_contents = set()
    for absolute in range(3, 1001):
        discriminant = -absolute
        if discriminant % 4 not in (0, 1):
            continue
        discriminants += 1
        for root in reduced_forms(discriminant):
            data = square_parametrization(root)
            target = data["target"]
            assert target == compose_forms(root, root)
            first, second, third = data["polynomials"]
            a, b, c = target
            assert sum(x*y for x, y in zip(first, cross(second, third))) == -discriminant
            # The stronger identity uses three independent variables, before
            # imposing the rank-one equation r0*r2=r1^2.
            for i in range(3):
                for j in range(i, 3):
                    coefficient = (2*a*first[i]*first[j]+b*(first[i]*second[j]+first[j]*second[i])
                                   +2*c*second[i]*second[j]-2*third[i]*third[j])
                    expected = (-2*discriminant if i == j == 1 else
                                discriminant if (i, j) == (0, 2) else 0)
                    assert coefficient == expected
            terms = [polynomial_product(first, first), polynomial_product(first, second),
                     polynomial_product(second, second), polynomial_product(third, third)]
            assert all(a*aa+b*bb+c*cc == zz for aa, bb, cc, zz in zip(*terms))
            identities += 1
            for u in range(-4, 5):
                for v in range(-4, 5):
                    if gcd(u, v) != 1:
                        continue
                    result = parametrized_point(data, u, v)
                    all_contents.add(result["content"])
                    nontrivial += result["content"] > 1
                    forward += 1
                    if forward % 19 == 0:
                        x, y, denominator = result["primitive_point"]
                        witness = representation_from_denominator(data, u, v, 6*denominator)
                        assert witness["representation"] == (6*x, 6*y)
                        normalized += 1
            # Independently enumerate the conic, without using the forward map.
            for x in range(-10, 11):
                for y in range(-10, 11):
                    value = evaluate(target, x, y)
                    z = isqrt(value)
                    if not z or z*z != value:
                        continue
                    result = parameters_from_point(data, (x, y, z))
                    u, v = result["parameters"]
                    assert gcd(u, v) == 1
                    assert evaluate(root, u, v) == result["content"]*z//result["input_content"]
                    inverse += 1
            forms_checked += 1
    return {"absolute_discriminant_maximum": 1000, "discriminants": discriminants,
            "forms_and_coefficient_identities": forms_checked,
            "independent_composition_comparisons": identities,
            "primitive_parameter_evaluations": forward,
            "nontrivial_content_evaluations": nontrivial,
            "distinct_contents_observed": len(all_contents),
            "independently_enumerated_conic_points": inverse,
            "divisible_denominator_scalings": normalized}


def check_splitting():
    roots = split = kernels = 0
    examples = []
    for integer in range(3, 1000, 2):
        discriminant = -4*integer
        for root in reduced_forms(discriminant):
            data = square_parametrization(root)
            if data["target"] != (1, 0, integer):
                continue
            points = [content_factor(data, (sign, 0, 1), integer) for sign in (1, -1)]
            factors = [point["factor"] for point in points if point["factor"]]
            expected = ambiguous_factor(root, integer)
            assert bool(factors) == bool(expected)
            assert all(integer % factor == 0 for factor in factors)
            roots += 1
            split += bool(factors)
            kernels += not factors
            if integer in (15, 21, 65) and root[0] != 1:
                examples.append({"integer": integer, "root": root,
                                 "unit_point_contents": [p["content"] for p in points],
                                 "proper_factors": factors})
    return {"odd_integer_maximum": 999, "ambiguous_roots": roots,
            "roots_giving_a_proper_factor": split, "kernel_roots": kernels,
            "independent_boundary_formula_comparisons": roots, "examples": examples}


def check_large_and_unreduced():
    state = Path(__file__).resolve().parents[1]/"evidence"
    records = []
    for bits in (60, 80):
        certificate = json.loads((state/f"arithmetic-norm-word-sieve-{bits}bit.json").read_text())
        integer = int(certificate["integer"])
        for encoded, order_log in certificate["extracted_two_group_basis"]:
            root = power_primitive(tuple(map(int, encoded)), 1 << (order_log-1))
            data = square_parametrization(root)
            assert data["target"] == (1, 0, integer)
            points = [content_factor(data, (sign, 0, 1), integer) for sign in (1, -1)]
            factors = sorted({p["factor"] for p in points if p["factor"]})
            records.append({"integer_bits": bits, "root": root,
                            "unit_point_contents": [p["content"] for p in points],
                            "proper_factors": factors})
        assert any(record["proper_factors"] for record in records if record["integer_bits"] == bits)
    transformed = 0
    for root in ((2, 1, 2), (4, 0, 3), (3, 2, 3), (7, 1, 8)):
        a, b, c = root
        for shift in (17, -1009, 1 << 80):
            changed = (a+b*shift+c*shift*shift, b+2*c*shift, c)
            data = square_parametrization(changed)
            for u, v in ((1, 0), (0, 1), (1, -1), (19, 7)):
                point = parametrized_point(data, u, v)
                inverse = parameters_from_point(data, point["primitive_point"])
                assert inverse["parameters"] in ((u, v), (-u, -v))
                transformed += 1
    # A known root gives an immediate rational point with a rough denominator.
    rough_root = (703135247427, -646306737026, 1867853905160)
    data = square_parametrization(rough_root)
    point = parametrized_point(data, 1, 0)
    denominator = point["primitive_point"][2]
    assert smooth_factors(denominator, 1000) is None
    return {"large_certificate_unit_points": records,
            "unreduced_root_evaluations": transformed,
            "known_root_rough_denominator": {"root": rough_root,
                "denominator": denominator, "smooth_bound_rejected": 1000}}


def check_modular_points():
    cases = 0
    partition_sizes = set()
    for absolute in range(3, 257):
        if (-absolute) % 4 not in (0, 1):
            continue
        for root in reduced_forms(-absolute):
            data = square_parametrization(root)
            for modulus in (2, 8, 9, 16, 25, 72, 1001, (1 << 127)*3**45*5**19):
                result = modular_unit_point(data, modulus)
                x, y = result["vector"]
                assert 0 <= x < modulus and 0 <= y < modulus
                assert (evaluate(data["target"], x, y)-1) % modulus == 0
                partition_sizes.add(sum(part > 1 for part in result["coprime_parts"]))
                cases += 1
    # Q is a nonprincipal square of a class of order three. It represents
    # one over every Z/mZ, while its minimum nonzero integer norm is two.
    data = square_parametrization((2, -1, 3))
    assert data["target"] == (2, 1, 3)
    minimum = min(evaluate(data["target"], x, y) for x in range(-10, 11)
                  for y in range(-10, 11) if x or y)
    assert minimum == 2
    examples = [modular_unit_point(data, modulus) for modulus in
                (2, 9, 18, 2**81*3**67*23**11)]
    return {"uniform_construction_cases": cases,
            "nontrivial_partition_sizes_observed": sorted(partition_sizes),
            "every_modulus_but_no_integer_unit_point": {
                "root": data["root"], "target": data["target"],
                "integer_minimum": minimum, "examples": examples}}


def main():
    print(json.dumps({"scope": "Exact conic transport and supplied-point extraction; no search bound",
                      "small_forms": check_small_forms(), "splitting": check_splitting(),
                      "large_and_unreduced": check_large_and_unreduced(),
                      "modular_unit_points": check_modular_points()}, indent=2))


if __name__ == "__main__":
    main()
