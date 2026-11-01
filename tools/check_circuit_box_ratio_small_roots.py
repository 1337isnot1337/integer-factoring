#!/usr/bin/env python3
"""Independent exact controls for the classical quarter-bit endpoint."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, isqrt
from pathlib import Path
import random

from circuit_box_ratio_small_roots import (
    bounded_integer_roots,
    dyadic_residue_factors,
    residue_short_polynomial,
    smallest_factor_by_small_roots,
    small_root_parameters,
)


def _multiply(first, second):
    result = [0] * (len(first) + len(second) - 1)
    for i, x in enumerate(first):
        for j, y in enumerate(second):
            result[i + j] += x * y
    return result


def _evaluate(polynomial, argument):
    return sum(coefficient * argument**i for i, coefficient in enumerate(polynomial))


def _smallest_factor(number):
    for divisor in range(2, isqrt(number) + 1):
        if number % divisor == 0:
            return divisor
    return number


def _bound_budget(bits, factor_exponent, precision, dimension, multiplicity):
    length_exponent = factor_exponent + 1 - precision
    return (
        4 * factor_exponent * multiplicity * dimension
        - dimension * (dimension - 1)
        - 2 * bits * multiplicity * (multiplicity + 1)
        - 2 * length_exponent * dimension * (dimension - 1)
    )


def run_checks(complete=True):
    rng = random.Random(2026110164283)
    counts = {"universal_parameter_inequalities": 0, "selected_parameter_inequalities": 0, "literal_integer_root_lists": 0, "repeated_and_rational_root_controls": 0, "literal_dyadic_residue_factor_lists": 0, "smallest_prime_factor_controls": 0, "independent_short_polynomial_certificates": 0}
    for bits in list(range(2, 513)) + [1024, 2048, 4096, 16384]:
        precision = (bits + 3) // 4 + 2
        dimension = 2 * bits
        overhead = (dimension ** (2 * dimension)).bit_length()
        exponents = range(precision, bits) if bits <= 512 else sorted({precision, bits // 2, 3 * bits // 4, bits - 1})
        for exponent in exponents:
            multiplicity = 2 * exponent
            assert 1 <= multiplicity < dimension
            budget = _bound_budget(bits, exponent, precision, dimension, multiplicity)
            assert budget >= overhead
            # A separate exact square-completion check proves the input
            # bound X < N**((k/n)**2), without floating powers.
            difference = Fraction((bits - 1) * exponent * exponent, bits * bits) - (exponent + 1 - precision)
            assert difference > 0
            counts["universal_parameter_inequalities"] += 1
    for _ in range(160):
        bits = rng.randrange(8, 260)
        precision = (bits + 3) // 4 + 2
        if precision >= bits:
            continue
        exponent = rng.randrange(precision, bits)
        dimension, multiplicity = small_root_parameters(bits, exponent, precision)
        budget = _bound_budget(bits, exponent, precision, dimension, multiplicity)
        assert dimension ** (2 * dimension) < 1 << budget
        assert dimension <= 2 * bits
        counts["selected_parameter_inequalities"] += 1

    for _ in range(800):
        polynomial = [rng.randrange(-15, 16) for _ in range(rng.randrange(1, 10))]
        if not any(polynomial):
            polynomial[0] = 1
        lower, upper = sorted(rng.randrange(-30, 31) for _ in range(2))
        actual = bounded_integer_roots(polynomial, lower, upper)
        expected = [x for x in range(lower, upper) if _evaluate(polynomial, x) == 0]
        assert actual["roots"] == expected
        counts["literal_integer_root_lists"] += 1
    for _ in range(120):
        polynomial, expected = [1], set()
        for _ in range(rng.randrange(1, 7)):
            root = rng.randrange(-20, 21)
            expected.add(root)
            for _ in range(rng.randrange(1, 4)):
                polynomial = _multiply(polynomial, [-root, 1])
        polynomial = _multiply(polynomial, [1, 2])
        polynomial = _multiply(polynomial, [1, 2])
        polynomial = _multiply(polynomial, [-2, 0, 1])
        actual = bounded_integer_roots(polynomial, -30, 31)
        assert actual["roots"] == sorted(expected)
        counts["repeated_and_rational_root_controls"] += 1
    huge_polynomial = [1]
    for factor in ([1, 2], [1, 2], [-1, 1], [-1, 1], [-1, 1], [-2, 0, 1], [-(1 << 256), 1]):
        huge_polynomial = _multiply(huge_polynomial, factor)
    huge_isolation = bounded_integer_roots(huge_polynomial, -(1 << 300), 1 << 300)
    assert huge_isolation["roots"] == [1, 1 << 256]
    assert huge_isolation["visited_intervals"] <= 2 * (len(huge_polynomial) - 1) * 602 + 1

    for _ in range(48):
        number = rng.randrange(3, 5000) | 1
        exponent = rng.randrange(number.bit_length())
        precision = (number.bit_length() + 3) // 4 + 2
        modulus = 1 << precision
        residue = rng.randrange(1, modulus, 2)
        actual = dyadic_residue_factors(number, exponent, residue)
        expected = [a for a in range(1 << exponent, min(1 << (exponent + 1), number + 1)) if a % modulus == residue and number % a == 0]
        assert actual["factors"] == expected
        counts["literal_dyadic_residue_factor_lists"] += 1

    for number in range(2, 2000):
        assert smallest_factor_by_small_roots(number)["factor"] == _smallest_factor(number)
        counts["smallest_prime_factor_controls"] += 1

    large_controls = []
    for exponent_three, exponent_five, guard in ((8, 6, 2), (24, 16, 8), (60, 40, 16), (80, 55, 16)):
        first, second = 3**exponent_three, 5**exponent_five
        number = first * second
        exponent = first.bit_length() - 1
        precision = (number.bit_length() + 3) // 4 + guard
        modulus = 1 << precision
        residue = first % modulus
        witness = residue_short_polynomial(number, exponent, precision, residue)
        dimension, multiplicity = witness["dimension"], witness["multiplicity"]
        bound = witness["root_bound"]
        constant = witness["linear_constant"]
        unscaled_rows = []
        for degree in range(dimension):
            power = min(degree, multiplicity)
            shift = degree - power
            row = [0] * dimension
            for j in range(power + 1):
                row[j + shift] = number**max(multiplicity - degree, 0) * comb(power, j) * constant ** (power - j)
            unscaled_rows.append(row)
        rebuilt = [sum(witness["lattice_combination"][i] * unscaled_rows[i][j] for i in range(dimension)) for j in range(dimension)]
        while rebuilt and not rebuilt[-1]:
            rebuilt.pop()
        assert rebuilt == witness["polynomial"]
        exact_norm = sum(abs(value) * bound**i for i, value in enumerate(rebuilt))
        assert exact_norm == witness["weighted_l1_norm"] < (1 << exponent) ** multiplicity
        counts["independent_short_polynomial_certificates"] += 1
        expected = sorted(3**i * 5**j for i in range(exponent_three + 1) for j in range(exponent_five + 1) if (1 << exponent) <= 3**i * 5**j < (1 << (exponent + 1)) and 3**i * 5**j % modulus == residue)
        isolated = bounded_integer_roots(rebuilt, 0, bound)
        recovered = sorted(residue + modulus * root for root in isolated["roots"] if (1 << exponent) <= residue + modulus * root < (1 << (exponent + 1)) and number % (residue + modulus * root) == 0)
        assert recovered == expected and first in recovered
        for factor in expected:
            root = (factor - residue) // modulus
            assert all(_evaluate(row, root) % factor**multiplicity == 0 for row in unscaled_rows)
            assert _evaluate(rebuilt, root) == 0
        negative_residue = (residue + 2) % modulus
        negative = dyadic_residue_factors(number, exponent, negative_residue, guard)
        negative_expected = sorted(3**i * 5**j for i in range(exponent_three + 1) for j in range(exponent_five + 1) if (1 << exponent) <= 3**i * 5**j < (1 << (exponent + 1)) and 3**i * 5**j % modulus == negative_residue)
        assert negative["factors"] == negative_expected
        large_controls.append({
            "input_bits": number.bit_length(), "number": number,
            "factorization": {"3": exponent_three, "5": exponent_five},
            "dyadic_exponent": exponent, "ratio_precision": precision, "guard_bits": guard,
            "residue": residue, "interval_length_bits": exponent + 1,
            "independent_divisor_candidates": (exponent_three + 1) * (exponent_five + 1),
            "recovered_factors": recovered, "negative_residue_factors": negative_expected,
            "short_polynomial": witness, "root_isolation": isolated,
            "negative_query": negative,
        })

    complete_examples = []
    if complete:
        for number in (257 * 263, 269 * 1009):
            result = smallest_factor_by_small_roots(number)
            assert result["factor"] == _smallest_factor(number)
            assert result["lattice_queries"] > 0
            complete_examples.append({"number": number, "input_bits": number.bit_length(), **result})

    root = Path(__file__).resolve().parent.parent
    paths = ["tools/circuit_box_ratio_small_roots.py", "tools/check_circuit_box_ratio_small_roots.py"]
    return {
        "status": "passed", "seed": 2026110164283, "complete_examples_included": complete,
        "counts": counts, "large_local_controls": large_controls,
        "huge_interval_root_isolation": {"interval_endpoint_bits": 301, "polynomial_degree": len(huge_polynomial) - 1, **huge_isolation},
        "complete_examples": complete_examples,
        "proved_search_cost": "2**(n/4)*poly(n) time and polynomial bit space; classical bound, no repository score increase",
        "source_sha256": {path: hashlib.sha256((root / path).read_bytes()).hexdigest() for path in paths},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-complete", action="store_true")
    args = parser.parse_args()
    result = run_checks(not args.skip_complete)
    if args.output:
        args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": result["status"], **result["counts"], "large_control_bits": [entry["input_bits"] for entry in result["large_local_controls"]], "complete_examples": result["complete_examples"]}))


if __name__ == "__main__":
    main()
