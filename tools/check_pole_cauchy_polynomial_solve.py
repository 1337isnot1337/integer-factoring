#!/usr/bin/env python3
"""Independent exact checks for implicit polynomial Cauchy solves."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import random

from check_pole_cauchy_conditioning import expand, inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from pole_cauchy_conditioning import balanced_padding
from pole_cauchy_polynomial_solve import (inverse_constant_vector_entry,
                                         inverse_polynomial_vector_entry,
                                         inverse_polynomial_bilinear_form)


def scaled(value: Quad, exponent: int) -> Quad:
    return value / (1 << exponent) if exponent >= 0 else value * (1 << -exponent)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(8311)
    constant_checks = 0
    polynomial_checks = 0
    bilinear_checks = 0
    for index in range(150):
        n = rng.randrange(31, 200, 2)
        if index < 100:
            aa, bb, _ = balanced_padding(
                n, (rng.randrange(1, 8, 2), rng.randrange(1, 6)),
                (rng.randrange(1, 8, 2), rng.randrange(1, 6)))
        else:
            size = rng.randrange(1, 6)
            aa = tuple((a, 1) for a in sorted(rng.sample(range(1, 32, 2), size)))
            bb = tuple((b, 1) for b in sorted(rng.sample(range(1, 32, 2), size)))
        rows, cols = expand(aa), expand(bb)
        h = max(1, max(abs(n - a * b) for a in rows for b in cols).bit_length())
        k = h + 3
        inv = inverse(matrix_at(n, rows, cols, k))
        for column, b in enumerate(cols):
            actual = sum(inv[column], Quad())
            oracle = inverse_constant_vector_entry(n, aa, bb, b, k, 8)
            assert actual.valuation() == oracle["valuation"]
            assert mod_pair(scaled(actual, oracle["valuation"]), 8) == tuple(oracle["unit"])
            constant_checks += 1
            for degree in (0, 1, 3, 6):
                coefficients = tuple(rng.randrange(-5, 6) for _ in range(degree + 1))
                values = [sum(c * a ** d for d, c in enumerate(coefficients)) for a in rows]
                actual = sum((entry * value for entry, value in zip(inv[column], values)), Quad())
                oracle = inverse_polynomial_vector_entry(n, aa, bb, b, k, coefficients, 8)
                assert mod_pair(scaled(actual, oracle["scale_exponent"]), 8) == tuple(oracle["residue"])
                polynomial_checks += 1
        for row_degree, column_degree in ((0, 0), (0, 4), (4, 0), (3, 5)):
            row_coefficients = tuple(rng.randrange(-5, 6) for _ in range(row_degree + 1))
            column_coefficients = tuple(rng.randrange(-5, 6) for _ in range(column_degree + 1))
            row_values = [sum(c * a ** d for d, c in enumerate(row_coefficients)) for a in rows]
            column_values = [sum(c * b ** e for e, c in enumerate(column_coefficients)) for b in cols]
            actual = sum((inv[j][i] * row_values[i] * column_values[j]
                          for i in range(len(rows)) for j in range(len(cols))), Quad())
            oracle = inverse_polynomial_bilinear_form(n, aa, bb, k, row_coefficients, column_coefficients, 8)
            assert mod_pair(actual * (1 << k), 8) == tuple(oracle["residue"])
            bilinear_checks += 1
    n = (1 << 401) + 45
    aa, bb, modulus = balanced_padding(n, (1, (1 << 180) + 7), (3, (1 << 200) + 13))
    max_a = aa[-1][0] + 2 * (aa[-1][1] - 1)
    max_b = bb[-1][0] + 2 * (bb[-1][1] - 1)
    k = max(abs(n - aa[0][0] * bb[0][0]), abs(n - max_a * max_b)).bit_length() + 12
    coefficients = (3, -2, 0, 5, 1, 0, -1, 0, 2)
    long = []
    for b in (3, 5, 3 + (1 << 151)):
        low = inverse_polynomial_vector_entry(n, aa, bb, b, k, coefficients, 12)
        high = inverse_polynomial_vector_entry(n, aa, bb, b, k, coefficients, 16)
        assert low["scale_exponent"] == high["scale_exponent"]
        assert low["residue"] == [v % (1 << 12) for v in high["residue"]]
        long.append({"column": b, **low})
    bilinear = inverse_polynomial_bilinear_form(n, aa, bb, k, coefficients, coefficients[::-1], 12)
    bilinear_high = inverse_polynomial_bilinear_form(n, aa, bb, k, coefficients, coefficients[::-1], 16)
    assert bilinear["residue"] == [v % (1 << 12) for v in bilinear_high["residue"]]
    absolute_bilinear = inverse_polynomial_bilinear_form(n, aa, bb, k, coefficients, coefficients[::-1], k + 12)
    assert bilinear["residue"] == [v % (1 << 12) for v in absolute_bilinear["residue"]]
    absolute_bilinear["absolute_precision"] = 12
    result = {"claim": "Polynomial-cost entries of inverse pole Cauchy matrices applied to polynomial vectors",
              "exact_constant_vector_entries": constant_checks,
              "exact_polynomial_vector_entries": polynomial_checks,
              "exact_polynomial_bilinear_forms": bilinear_checks,
              "matrix_controls": 150,
              "implicit_control": {"input_bits": n.bit_length(), "domain_size": (1 << 200) + 13,
                                    "left": aa, "right": bb, "k": k,
                                    "coefficients": coefficients, "queries": long,
                                    "bilinear_form": bilinear,
                                    "absolute_precision_bilinear_form": absolute_bilinear},
              "scope": "Polynomial-weight inverse moments sum both axes; the Hadamard-square trace and factor count remain open."}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
