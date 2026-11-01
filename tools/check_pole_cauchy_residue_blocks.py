#!/usr/bin/env python3
"""Exact Hermite selectors and arbitrary-precision inverse-block controls."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import random

from check_pole_cauchy_conditioning import expand, inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from pole_cauchy_quadratic_selectors import inverse_residue_block_parity
from pole_cauchy_residue_blocks import inverse_residue_block, selector_rounding_coefficients

ROOT = Path(__file__).resolve().parents[1]


def controls():
    coefficient_identities = scalar_selector_values = quadratic_selector_values = 0
    for precision in range(1, 17):
        coefficients = selector_rounding_coefficients(precision)
        direct = [0] * (2 * precision)
        for j in range(precision):
            for ell in range(j + 1):
                direct[precision + ell] += comb(precision + j - 1, j) * comb(j, ell) * (-1) ** ell
        assert coefficients == {j: value for j, value in enumerate(direct) if value}
        around_one = [sum(value * comb(j, ell) * (-1) ** ell
                          for j, value in coefficients.items() if j >= ell)
                      for ell in range(precision)]
        assert around_one == [1] + [0] * (precision - 1)
        coefficient_identities += 1
        for value in range(-16, 17):
            result = sum(coefficient * pow(value, degree, 1 << precision)
                         for degree, coefficient in coefficients.items()) % (1 << precision)
            assert result == value % 2
            scalar_selector_values += 1
        for real in range(-3, 4):
            for imaginary in range(-3, 4):
                value = Quad(Fraction(real), Fraction(2 * imaginary))
                power = Quad(Fraction(1))
                result = Quad()
                for degree in range(2 * precision):
                    if degree in coefficients:
                        result += coefficients[degree] * power
                    power *= value
                assert mod_pair(result, precision) == (real % 2, 0)
                quadratic_selector_values += 1
    rng = random.Random(2026110112)
    dense_blocks = high_bit_unit_controls = parity_agreements = 0
    examples = []
    for index in range(18):
        size = rng.randrange(1, 4)
        left = ((rng.randrange(1, 8, 2), size),)
        right = tuple((b, 1) for b in sorted(rng.sample(list(range(1, 16, 2)), size)))
        aa, bb = expand(left), expand(right)
        n = rng.randrange(21, 150, 2)
        k = max(abs(n - a * b) for a in aa for b in bb).bit_length() + 1
        precision = 1 + index % 3
        row_bits, column_bits = index % 4, (index // 4) % 4
        row_residue, column_residue = rng.randrange(8), rng.randrange(8)
        dense = inverse(matrix_at(n, aa, bb, k))
        actual = inverse_residue_block(n, left, right, k, row_residue, row_bits,
                                       column_residue, column_bits, precision)
        expected = sum((dense[j][i] for i, a in enumerate(aa) for j, b in enumerate(bb)
                        if (a - row_residue) % (1 << row_bits) == 0
                        and (b - column_residue) % (1 << column_bits) == 0), Quad())
        expected *= 1 << actual["inverse_normalization_exponent"]
        assert tuple(actual["residue"]) == mod_pair(expected, precision)
        assert actual["repeated_pole_contractions"] == precision ** 2
        dense_blocks += 1
        if precision == 1:
            earlier = inverse_residue_block_parity(n, left, right, k, row_residue, row_bits,
                                                   column_residue, column_bits)
            assert earlier["parity"] == actual["residue"]
            parity_agreements += 1
        examples.append({"n": n, "size": size, "precision": precision,
                         "residue": actual["residue"],
                         "inverse_normalization_exponent": actual["inverse_normalization_exponent"],
                         "maximum_pole_order": actual["maximum_pole_order"],
                         "maximum_working_precision": max(x["working_precision"]
                                                           for x in actual["contraction_budgets"])})
    # At precision L+1 the original unit directions need not vanish.
    # These are access controls, not an algorithm to count those directions.
    domain = ((1, 3),)
    dense = inverse(matrix_at(15, (1, 3, 5), (1, 3, 5), 5))
    high_bit_examples = []
    for a, b in ((1, 1), (3, 5), (5, 3)):
        actual = inverse_residue_block(15, domain, domain, 5, a, 3, b, 3, 5)
        expected = dense[(b - 1) // 2][(a - 1) // 2] * 16
        assert tuple(actual["residue"]) == mod_pair(expected, 5)
        if a * b == 15:
            assert actual["residue"] == [16, 0]
        high_bit_examples.append({"a": a, "b": b, "precision": 5,
                                  "normalized_inverse_residue": actual["residue"]})
        high_bit_unit_controls += 1
    return {"exact_selector_polynomial_identities": coefficient_identities,
            "integer_selector_values": scalar_selector_values,
            "quadratic_selector_values": quadratic_selector_values,
            "dense_residue_block_sums": dense_blocks,
            "earlier_parity_interface_agreements": parity_agreements,
            "higher_precision_unit_direction_controls": high_bit_unit_controls,
            "examples": examples, "higher_precision_examples": high_bit_examples}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"claim": "Arbitrary-precision residue-block sums from polynomial-degree rational selectors",
              "controls": controls(), "factoring_exponent_improved": False,
              "source_sha256": {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "pole_cauchy_residue_blocks.py",
                                             "pole_cauchy_quadratic_selectors.py")},
              "scope": "An implicit block-sum oracle does not itself count Smith unit directions."}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["controls"], indent=2))


if __name__ == "__main__":
    main()
