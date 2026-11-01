#!/usr/bin/env python3
"""Sharp dyadic Bezout denominator controls for short arithmetic domains."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, prod
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[1]


def v2(value):
    if not value:
        return None
    return (value & -value).bit_length() - 1


def multiply(a, b, modulus=None):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return [x % modulus for x in result] if modulus else result


def factors_product(factors):
    result = [1]
    for factor in factors:
        result = multiply(result, factor)
    return result


def remainder(a, polynomial, modulus):
    assert polynomial[-1] == 1
    degree = len(polynomial) - 1
    result = [x % modulus for x in a]
    for j in range(len(result) - 1, degree - 1, -1):
        value = result[j]
        for i in range(degree + 1):
            result[j - degree + i] = (result[j - degree + i] - value * polynomial[i]) % modulus
    return (result + [0] * degree)[:degree]


def scaled_lagrange_inverse(n, aa, bb):
    """Compute 2^L Q^-1 mod P directly from all exact rational root values."""
    m = len(aa)
    terms, resultant_valuation = [], 0
    for a in aa:
        value = prod(n - a * b for b in bb)
        resultant_valuation += v2(value)
        derivative = prod(a - other for other in aa if other != a)
        denominator = value * derivative
        terms.append((a, denominator, v2(denominator)))
    scale = max(item[2] for item in terms)
    modulus = 1 << (scale + 1)
    polynomial = factors_product([[-a, 1] for a in aa])
    inverse = [0] * m
    for a, denominator, order in terms:
        quotient = [0] * m
        quotient[-1] = polynomial[-1] % modulus
        for i in range(m - 2, -1, -1):
            quotient[i] = (polynomial[i + 1] + a * quotient[i + 1]) % modulus
        coefficient = ((1 << (scale - order)) * pow(denominator >> order, -1, modulus)) % modulus
        inverse = [(x + coefficient * y) % modulus for x, y in zip(inverse, quotient)]
    q = factors_product([[n, -b] for b in bb])
    assert remainder(multiply(q, inverse, modulus), polynomial, modulus) == [1 << scale] + [0] * (m - 1)
    minimum = min(v2(x) for x in inverse if x)
    assert resultant_valuation == m * m
    return {"dyadic_denominator": scale - minimum, "interpolation_scale": scale,
            "ordinary_resultant_valuation": resultant_valuation,
            "scaled_inverse_modulus_bits": scale + 1, "scaled_inverse_coefficients": inverse}


def normalized_hensel_inverse(n, aa, bb):
    """Independent quotient-algebra inverse after X=1+2Z."""
    m = len(aa)
    polynomial = factors_product([[-((a - 1) // 2), 1] for a in aa])
    q = factors_product([[(n - b) // 2, -b] for b in bb])
    assert [x % 2 for x in polynomial] == [0] * m + [1]
    inverse = [comb(m + i - 1, i) % 2 for i in range(m)]
    precision, target = 1, m + 1
    while precision < target:
        precision = min(target, 2 * precision)
        modulus = 1 << precision
        error = remainder(multiply(q, inverse, modulus), polynomial, modulus)
        correction = [(-x) % modulus for x in error]
        correction[0] = (correction[0] + 2) % modulus
        inverse = remainder(multiply(inverse, correction, modulus), polynomial, modulus)
        assert remainder(multiply(q, inverse, modulus), polynomial, modulus) == [1] + [0] * (m - 1)
    result = max(m + i - v2(coefficient) for i, coefficient in enumerate(inverse) if coefficient)
    return {"dyadic_denominator": result, "normalized_inverse_precision": target,
            "normalized_inverse_coefficients": inverse,
            "normalized_inverse_parities": [x % 2 for x in inverse]}


def exact_fraction_inverse(n, aa, bb):
    m = len(aa)
    result = [Fraction(0)] * m
    for a in aa:
        quotient = factors_product([[-other, 1] for other in aa if other != a])
        denominator = prod(n - a * b for b in bb) * prod(a - other for other in aa if other != a)
        result = [x + Fraction(y, denominator) for x, y in zip(result, quotient)]
    return max(v2(x.denominator) for x in result if x)


def controls():
    rng = random.Random(2026110129)
    checked = sharp = fractions = 0
    examples = []
    for m in list(range(2, 19)) + [33]:
        schedules = [(3, [1 + 4 * i for i in range(m)], [1 + 4 * i for i in range(m)])]
        for _ in range(2):
            schedules.append((3 + 4 * rng.randrange(1, 100),
                              sorted(1 + 4 * i for i in rng.sample(range(3 * m), m)),
                              sorted(1 + 4 * i for i in rng.sample(range(4 * m), m))))
        for n, aa, bb in schedules:
            first = scaled_lagrange_inverse(n, aa, bb)
            second = normalized_hensel_inverse(n, aa, bb)
            assert first["dyadic_denominator"] == second["dyadic_denominator"]
            assert m <= first["dyadic_denominator"] <= 2 * m - 2
            assert second["normalized_inverse_parities"] == [comb(m + i - 1, i) % 2 for i in range(m)]
            if (m - 1) & (m - 2) == 0:
                assert first["dyadic_denominator"] == 2 * m - 2
                assert second["normalized_inverse_parities"][m - 2] == 1
                assert second["normalized_inverse_parities"][m - 1] == 0
                sharp += 1
            if m <= 5:
                assert exact_fraction_inverse(n, aa, bb) == first["dyadic_denominator"]
                fractions += 1
            checked += 1
            if n == 3:
                examples.append({"n": n, "dimension": m, "a_progression": [1, 4, m],
                                 "b_progression": [1, 4, m],
                                 "dyadic_denominator": first["dyadic_denominator"],
                                 "ordinary_resultant_valuation": first["ordinary_resultant_valuation"],
                                 "parity_profile": second["normalized_inverse_parities"],
                                 "interpolation_scale": first["interpolation_scale"]})
    irregular = {"n": 3, "rows": [1, 5, 9, 13], "columns": [1, 5, 9, 17]}
    irregular["lagrange"] = scaled_lagrange_inverse(irregular["n"], irregular["rows"], irregular["columns"])
    irregular["hensel"] = normalized_hensel_inverse(irregular["n"], irregular["rows"], irregular["columns"])
    assert irregular["lagrange"]["dyadic_denominator"] == irregular["hensel"]["dyadic_denominator"] == 6
    assert irregular["hensel"]["normalized_inverse_coefficients"] == [3, 6, 0, 6]
    s = 200
    m = (1 << s) + 1
    implicit = {"n": 3, "dimension": m, "a_progression": [1, 4, m],
                "b_progression": [1, 4, m], "endpoint_encoding_bits": (1 + 4 * (m - 1)).bit_length(),
                "exact_dyadic_denominator": 2 * m - 2,
                "ordinary_resultant_valuation": m * m,
                "denominator_exponent_encoding_bits": (2 * m - 2).bit_length(),
                "enumerated_nodes": 0,
                "scope": "Closed-form application of the proved family, not a dense inverse computation."}
    return {"independent_lagrange_and_hensel_cases": checked,
            "sharp_power_plus_one_cases": sharp,
            "exact_rational_coefficient_cases": fractions,
            "examples": examples, "irregular_node_counterexample": irregular, "implicit_family": implicit}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"claim": "Sharp dyadic Bezout denominator 2m-2 for m=2^s+1 in the declared residue-ball family",
              "controls": controls(), "factoring_exponent_improved": False,
              "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result["controls"].items() if isinstance(value, int)}, indent=2))


if __name__ == "__main__":
    main()
