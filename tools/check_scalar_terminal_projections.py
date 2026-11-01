#!/usr/bin/env python3
"""Exact controls for the maximal scalar terminal-line projection engine."""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
import json
from pathlib import Path

import check_cancellation_module_automata as arithmetic
import check_minimal_boundary_modules as linear
import check_triangular_auxiliary_row_fiber as reflection
from check_protected_subtensor_pressure import BASE_PROFILE, LETTER_PROFILES
from scalar_terminal_projection import (
    terminal_scalar_projection, scalar_pressure,
    boundary_scalar_projections, boundary_rate_lower_bound,
)


ROOT = Path(__file__).resolve().parent.parent
ZERO, ONE = arithmetic.ZERO, arithmetic.ONE


def tensor_product(left, right):
    return tuple(tuple(arithmetic.rmul(a, b) for a in row_left for b in row_right)
                 for row_left in left for row_right in right)


def rational_kernel_control():
    local, degree = arithmetic.elementary(2)
    left = (3, 2, 4, 8)  # 1+z,z,z^2,z^3
    right = (1, 3, 0, 2)
    for (x, y, u, v), matrix in local.items():
        a, b = 2 * x + y, 2 * u + v
        matrix[1][0] = arithmetic.pmul(3, arithmetic.pmul(left[a], right[b]))
        matrix[0][1] = 3 * int(x == y == u == v)
    degree = max(value.bit_length() - 1 for matrix in local.values()
                 for row in matrix for value in row)
    certificate = terminal_scalar_projection(local, degree)
    projection, scalar = certificate["projection"], certificate["scalar_matrix"]
    assert certificate["projection_rows"] == 3
    assert certificate["orders"] == (0, 0, 0)
    assert any(value[1] != 1 for row in projection for value in row), \
        "this control must exercise non-polynomial regular kernel coordinates"

    width = 4
    word = [[ZERO] * (width * width) for _ in range(width * width)]
    labels = tuple(product(range(2), repeat=4))
    for first, second in product(labels, repeat=2):
        a1, b1 = 2 * first[0] + first[1], 2 * first[2] + first[3]
        a2, b2 = 2 * second[0] + second[1], 2 * second[2] + second[3]
        matrices = [tuple(tuple(arithmetic.rat(value) for value in row)
                          for row in local[key]) for key in (first, second)]
        word[width * a1 + a2][width * b1 + b2] = \
            linear.matrix_product(*matrices)[0][0]
    projected_word = linear.matrix_product(tensor_product(projection, projection), word)
    assert projected_word == tensor_product(scalar, scalar)
    return {
        "auxiliary_dimension": 2,
        "physical_width": 4,
        "maximal_projection_rows": 3,
        "orders": certificate["orders"],
        "regular_nonpolynomial_coordinates_used": True,
        "local_terminal_identity_entries": certificate["regular_terminal_identity_entries"],
        "literal_two_site_projected_entries": 9 * 16,
    }


def transform_endpoint(local, heights, infinity):
    result = {}
    degree = 0
    for labels, matrix in local.items():
        height = reflection.interpolation.height(labels, heights)
        changed = []
        for row in matrix:
            changed_row = []
            for value in row:
                polynomial = 0
                for k in range(value.bit_length()):
                    if value >> k & 1:
                        exponent = height - k if infinity else height + k
                        assert exponent >= 0
                        polynomial ^= 1 << exponent
                degree = max(degree, polynomial.bit_length() - 1)
                changed_row.append(polynomial)
            changed.append(changed_row)
        result[labels] = changed
    return result, degree


def reflection_control():
    values = reflection.reflection_seed()
    local = {}
    endpoint_checks = 0
    for labels, value in zip(reflection.LABELS, values):
        x, y, u, v = labels
        upper = 3 * int(x == y == u == v == 0)
        lower = 3 * int(x == 3 and y == u)
        matrix = [[value, upper], [lower, value]]
        for i, j in product(range(2), repeat=2):
            assert matrix[i][j].bit_count() % 2 == int(i == j and x == v and y == u)
            endpoint_checks += 1
        local[labels] = matrix
    matrices = [tuple(tuple(arithmetic.rat(value) for value in row)
                      for row in local[key]) for key in ((0, 0, 0, 0), (3, 0, 0, 0))]
    assert linear.matrix_product(*matrices) != linear.matrix_product(*reversed(matrices))
    records, rates = [], []
    point = Fraction(5)
    for heights, cut, infinity in (
            (reflection.ZERO_HEIGHTS, 0, False),
            (reflection.INFINITY_HEIGHTS, 1, True)):
        transformed, degree = transform_endpoint(local, heights, infinity)
        certificate = terminal_scalar_projection(transformed, degree, alphabet=4, cut=cut)
        assert certificate["projection_rows"] == certificate["scalar_rank"] == 12
        expected = BASE_PROFILE if infinity else tuple(value + 4 for value in BASE_PROFILE)
        assert certificate["orders"] == expected
        rates.append(scalar_pressure(certificate, point))
        records.append({
            "endpoint": "infinity" if infinity else "zero",
            "projection_rows": certificate["projection_rows"],
            "orders": certificate["orders"],
            "local_terminal_identity_entries": certificate["regular_terminal_identity_entries"],
        })
    g = (1 + point ** -1) * (1 + point ** -3) * (1 + point ** -5 + point ** -7)
    assert point ** 4 * rates[0] * rates[1] == g * g
    return {
        "auxiliary_dimension": 2,
        "strict_endpoint_entries": endpoint_checks,
        "auxiliary_symbols_do_not_commute": True,
        "upper_offdiagonal_alters_vacuum_letter": True,
        "lower_offdiagonal_supported_on_last_letter": True,
        "full_matrix_scalar_rectangle_required": False,
        "endpoints": records,
        "point": str(point),
        "raw_product_lower_bound": str(g * g),
    }


def boundary_extraction_control():
    local, degree = arithmetic.elementary(2)
    scalar, _ = arithmetic.elementary(1)
    for (x, y, u, v), matrix in local.items():
        defect = 3 * int(x == u and y == v)
        for i, j in product(range(2), repeat=2):
            matrix[i][j] ^= defect
    # N is the all-one 2-by-2 matrix, N^2=0.  The actual terminal e0
    # has no nonzero left scalar projection in this physical cut, while
    # v=(1,1) is a reachable and visible scalar line for every symbol.
    direct = terminal_scalar_projection(local, degree, cut=0, infinity=True)
    assert direct["projection_rows"] == 0
    screen = boundary_scalar_projections(local, degree, cut=0, infinity=True)
    assert screen["reachable_dimension"] == screen["observable_dimension"] == 2
    assert len(screen["certificates"]) == 3 and not screen["excluded_lines"]
    recovered = next(item["certificate"] for item in screen["certificates"] if item["line"] == 3)
    assert recovered["projection_rows"] == 4 and recovered["orders"] == (0, 1, 1, 2)
    point = Fraction(5)
    assert boundary_rate_lower_bound(screen, point) == point + 2 + 1 / point

    suffix = (0, 1, 0, 1)  # C=1,E=0, so K_suffix=(1+z)X.
    suffix_matrix = tuple(tuple(arithmetic.rat(value) for value in row) for row in local[suffix])
    factor = arithmetic.rat(3)
    labels = tuple(product(range(2), repeat=4))
    identities = 0
    for first, second in product(labels, repeat=2):
        matrices = [tuple(tuple(arithmetic.rat(value) for value in row)
                          for row in local[key]) for key in (first, second)]
        word = linear.matrix_product(*matrices)
        prolonged = linear.matrix_product(word, suffix_matrix)
        actual = arithmetic.radd(arithmetic.rmul(factor, word[0][0]), prolonged[0][0])
        expected = arithmetic.rmul(factor, arithmetic.rat(arithmetic.pmul(
            scalar[first][0][0], scalar[second][0][0])))
        assert actual == expected
        identities += 1

    # A second scalar component cannot be used when the actual boundaries
    # never reach or observe it, even though it has a full local projection.
    hidden, _ = arithmetic.elementary(2)
    for matrix in hidden.values():
        matrix[1][1] = 1
    hidden_screen = boundary_scalar_projections(hidden, 1)
    assert hidden_screen["reachable_dimension"] == hidden_screen["observable_dimension"] == 1
    assert [item["line"] for item in hidden_screen["certificates"]] == [1]
    assert len(hidden_screen["excluded_lines"]) == 2
    return {
        "direct_terminal_projection_rows": 0,
        "extracted_scalar_projection_rows": 4,
        "reachable_dimension": 2,
        "observable_dimension": 2,
        "extracted_line": 3,
        "uniform_rate_lower_bound_at_five": str(boundary_rate_lower_bound(screen, point)),
        "different_length_identity": "(1+z) Gamma_t + Gamma_(t+1,suffix) = (1+z) T^tensor(t)",
        "literal_different_length_identity_entries": identities,
        "hidden_unreachable_component_rejected": True,
    }


def vacuum_extension_controls():
    """A visible scalar line survives arbitrary one-letter perturbations.

    The polynomial implementation is z times the Laurent vacuum extension.
    A separate perturbation on x=a makes the auxiliary symbols noncommute.
    """
    scalar = dict(zip(reflection.LABELS, reflection.reflection_seed()))
    records = []
    point = Fraction(5)
    for letter in (1, 2, 3):
        local = {}
        for key, value in scalar.items():
            x, y, u, v = key
            vacuum = arithmetic.pmul(3, value) if x == 0 else 0
            perturbation = 3 * int(x == letter and y == u)
            local[key] = [[(value << 1) ^ vacuum, vacuum],
                          [vacuum ^ perturbation, (value << 1) ^ vacuum]]
            for i, j in product(range(2), repeat=2):
                assert local[key][i][j].bit_count() % 2 == \
                    int(i == j and x == v and y == u)
        # The unchanged vacuum symbol extracts w=(1,1), independently
        # of the perturbation on the nonzero physical letter.
        suffix = local[(0, 0, 0, 0)]
        assert suffix == [[16, 48], [48, 16]]
        first = tuple(tuple(arithmetic.rat(value) for value in row) for row in suffix)
        second = tuple(tuple(arithmetic.rat(value) for value in row)
                       for row in local[(letter, 0, 0, 0)])
        assert linear.matrix_product(first, second) != linear.matrix_product(second, first)

        endpoint_records, rates = [], []
        infinity_heights = tuple(
            tuple(value + int(leg == 0) for value in heights)
            for leg, heights in enumerate(reflection.INFINITY_HEIGHTS))
        for heights, cut, infinity in (
                (reflection.ZERO_HEIGHTS, 0, False),
                (infinity_heights, 1, True)):
            transformed, degree = transform_endpoint(local, heights, infinity)
            screen = boundary_scalar_projections(transformed, degree, alphabet=4, cut=cut)
            assert screen["reachable_dimension"] == screen["observable_dimension"] == 2
            certificate = next(item["certificate"] for item in screen["certificates"]
                               if item["line"] == 3)
            expected = LETTER_PROFILES[letter] if infinity else \
                tuple(value + 5 for value in LETTER_PROFILES[letter])
            assert certificate["projection_rows"] == certificate["scalar_rank"] == 12
            assert certificate["orders"] == expected
            rates.append(scalar_pressure(certificate, point))
            endpoint_records.append({
                "endpoint": "infinity" if infinity else "zero",
                "projection_rows": certificate["projection_rows"],
                "orders": certificate["orders"],
                "reachable_dimension": screen["reachable_dimension"],
                "observable_dimension": screen["observable_dimension"],
            })
        pressure = sum(point ** -order for order in LETTER_PROFILES[letter])
        assert point ** 5 * rates[0] * rates[1] == pressure ** 2

        protected = [key for key in reflection.LABELS if key[0] != letter]
        identities = 0
        for left, right in product(protected, repeat=2):
            a, b = local[left], local[right]
            word00 = arithmetic.pmul(a[0][0], b[0][0]) ^ arithmetic.pmul(a[0][1], b[1][0])
            word01 = arithmetic.pmul(a[0][0], b[0][1]) ^ arithmetic.pmul(a[0][1], b[1][1])
            prolonged00 = arithmetic.pmul(word00, suffix[0][0]) ^ \
                arithmetic.pmul(word01, suffix[1][0])
            actual = prolonged00 ^ (word00 << 5)
            expected = arithmetic.pmul(48, arithmetic.pmul(scalar[left], scalar[right]) << 2)
            assert actual == expected
            identities += 1
        records.append({
            "perturbed_nonzero_letter": letter,
            "auxiliary_symbols_do_not_commute": True,
            "strict_endpoint_entries": 1024,
            "extracted_line": 3,
            "endpoints": endpoint_records,
            "different_length_identity":
                "Gamma_(t+1,suffix) + z^5 Gamma_t = (z^4+z^5) (z T)^tensor(t) on x != a",
            "literal_two_site_different_length_entries": identities,
            "raw_pressure_product_at_five": str(pressure ** 2),
            "all_scale_rho_floor_from_protected_letter_certificate":
                "5/7" if letter == 2 else "8/11",
        })
    return records


def run():
    scalar, degree = arithmetic.elementary(2)
    control = terminal_scalar_projection(scalar, degree, cut=1)
    assert control["projection_rows"] == 4 and control["orders"] == (0, 1, 1, 2)
    return {
        "scope": "maximal regular left projection for a specified auxiliary terminal line",
        "scalar_control": {"projection_rows": 4, "orders": control["orders"]},
        "rational_kernel_control": rational_kernel_control(),
        "reflection_terminal_line_control": reflection_control(),
        "bounded_boundary_extraction": boundary_extraction_control(),
        "vacuum_extension_single_letter_controls": vacuum_extension_controls(),
        "floating_arithmetic_used": False,
        "factoring_exponent_improved": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "evidence" / "scalar-terminal-projections.json")
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
