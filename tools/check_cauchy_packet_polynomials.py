#!/usr/bin/env python3
"""Replay the affine Cauchy realization and exact finite packet witnesses."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import random

from cauchy_packet_polynomials import (
    add, affine_realization, cauchy_polynomials, content,
    direct_boolean_affine_polynomial, evaluate, modular_packet_polynomial,
    multiply, odd_root_polynomial, reconstructed_affine_polynomial,
    scale, shift, trim,
)


def binomial_polynomial(length, order):
    result = [1]
    for offset in range(order):
        result = multiply(result, [length - offset, 1])
    return scale(result, Fraction(1, factorial(order)))


def check():
    rng = random.Random(824)
    realization_cases = 0
    identity_coefficients = 0
    literal_boolean_cases = 0
    odd_center_cases = 0
    negative_controls = 0
    for length in (2, 4, 8, 16):
        basis = [binomial_polynomial(length, order) for order in range(length + 1)]
        cases = [basis[-1], odd_root_polynomial(length), [Fraction(1)]]
        for _ in range(4):
            polynomial = [0]
            for term in basis:
                polynomial = add(polynomial, scale(term, rng.randrange(-5, 6)))
            cases.append(polynomial)
        kernels = cauchy_polynomials(length)
        assert all(evaluate(polynomial, 0) == (difference == 0)
                   for difference, polynomial in kernels.items())
        for polynomial in cases:
            realization = affine_realization(polynomial, length)
            observed = reconstructed_affine_polynomial(realization)
            assert observed == trim(polynomial)
            identity_coefficients += max(len(observed), 2 * length)
            realization_cases += 1
        odd_polynomial = cases[1]
        assert evaluate(odd_polynomial, 0) == 1
        for center in range(-length, length):
            expected = 2 * length - 1 if center & 1 else 0
            assert content(shift(odd_polynomial, center, 2 * length)) == expected
            odd_center_cases += 1
        if length <= 4:
            realization = affine_realization(odd_polynomial, length)
            assert direct_boolean_affine_polynomial(realization) == odd_polynomial
            literal_boolean_cases += 1
            bad = dict(realization)
            bad["column"] = list(reversed(realization["column"]))
            assert reconstructed_affine_polynomial(bad) != odd_polynomial
            negative_controls += 1

    for invalid_polynomial, length in (([Fraction(1, 2)], 2), ([0, 0, 0, 1], 2)):
        try:
            affine_realization(invalid_polynomial, length)
        except ValueError:
            negative_controls += 1
        else:
            raise AssertionError("An invalid affine realization was accepted.")

    # A binary packet needs no noncanonical large coefficients: both tensors
    # are indicators of the alternating state 0101, and the input precision is 1.
    binary_column = [0] * 16
    binary_column[5] = 1
    binary_contents = []
    binary_coefficients = 0
    for center in (-3, -1, 1, 3):
        polynomial = modular_packet_polynomial(4, center, 20, 5, binary_column)
        nonzero = [value for value in polynomial if value]
        order = min((value & -value).bit_length() - 1 for value in nonzero)
        assert order >= 11
        binary_contents.append(order)
        binary_coefficients += 4 * (2 * 4 - 1) + 1
    assert min(binary_contents) == 11

    evidence_path = Path(__file__).resolve().parents[1] / "evidence/cauchy-packet-witness.json"
    witness = json.loads(evidence_path.read_text())
    length = witness["L"]
    row_state = witness["row_state"]
    column = witness["column_tensor"]
    required = witness["required_content"]
    assert column[row_state] & 1
    assert all(0 <= value < 1 << witness["coefficient_modulus_bits"] for value in column)
    centers = [center for center in range(-length, length)
               if center and (abs(center) & -abs(center)).bit_length() - 1 < witness["u"]]
    assert centers == witness["centers"]
    lifted_coefficients = 0
    for center in centers:
        polynomial = modular_packet_polynomial(length, center, required + 1, row_state, column)
        assert all(value % (1 << required) == 0 for value in polynomial)
        lifted_coefficients += length * (2 * length - 1) + 1

    # A literal canonical-lift negative control. For L=2 the exact affine
    # realization of 1-T^2 has A=(-3,-4,-2,-3), B=(3,-3,1,0).
    # Canonical representatives modulo2 instead have A=(1,0,0,1), B=(1,1,1,0).
    # At T=1, W_01=-3 is the only nonzero matrix entry. The value is then 6,
    # which is not divisible by8, although every original shifted coefficient is.
    realization = affine_realization([1, 0, -1], 2)
    row_values = [realization["constant"] + sum(
        realization["row"][x] for x in range(2) if state >> x & 1
    ) for state in range(4)]
    column_values = [realization["column_zero"], *realization["column"], 0]
    assert row_values == [-3, -4, -2, -3]
    assert column_values == [3, -3, 1, 0]
    canonical_row = [int(value) % 2 for value in row_values]
    canonical_column = [int(value) % 2 for value in column_values]
    matrix = [[evaluate(cauchy_polynomials(2)[y - x], 1) for y in range(2)] for x in range(2)]
    canonical_value = 0
    for x_state, column_weight in enumerate(canonical_column):
        images = [sum(matrix[x][y] for y in range(2) if x_state >> y & 1) for x in range(2)]
        for y_state, row_weight in enumerate(canonical_row):
            term = row_weight * column_weight
            for x, image in enumerate(images):
                term *= image if y_state >> x & 1 else 1 - image
            canonical_value += term
    assert canonical_value == 6
    assert canonical_value % 8 != 0
    negative_controls += 1

    sources = (Path(__file__), Path(__file__).with_name("cauchy_packet_polynomials.py"), evidence_path)
    return {
        "status": "pass",
        "affine_realizations": realization_cases,
        "global_identity_coefficients": identity_coefficients,
        "literal_boolean_sums": literal_boolean_cases,
        "odd_root_local_contents": odd_center_cases,
        "binary_packet_centers": [-3, -1, 1, 3],
        "binary_packet_contents": binary_contents,
        "binary_packet_coefficients_checked": binary_coefficients,
        "eight_wire_centers": centers,
        "eight_wire_coefficients_checked": lifted_coefficients,
        "eight_wire_required_content": required,
        "eight_wire_creation_precision": witness["coefficient_modulus_bits"],
        "eight_wire_origin_value": column[row_state],
        "canonical_lift_negative_control_value": int(canonical_value),
        "negative_controls": negative_controls,
        "source_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        "scope": (
            "The affine realization has an all-size proof. The two packet witnesses "
            "are finite exact certificates. Neither yields an asymptotic precision "
            "lower bound or an improved factoring score. The 33-bit witness uses "
            "33-bit input representatives and does not prove a 33-bit added loss."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = check()
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
