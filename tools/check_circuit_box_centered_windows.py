#!/usr/bin/env python3
"""Exact controls for centered low-residue divisor windows."""

import argparse
from hashlib import sha256
import json
from math import comb, isqrt
from pathlib import Path
from random import Random

from circuit_box_centered_windows import (
    centered_factor_rectangle, centered_window_factors,
    smallest_factor_by_centered_windows, window_parameters,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCES = [
    "tools/circuit_box_centered_windows.py",
    "tools/check_circuit_box_centered_windows.py",
    "tools/circuit_box_ratio_small_roots.py",
]


def _smallest(number):
    return next((value for value in range(2, isqrt(number) + 1) if number % value == 0), number)


def _known_divisors(exponent_three, exponent_five):
    powers_three = [3**i for i in range(exponent_three + 1)]
    powers_five = [5**j for j in range(exponent_five + 1)]
    return [a * b for a in powers_three for b in powers_five]


def _certificate(number, result):
    witness = result["witness"]
    dimension, multiplicity = witness["dimension"], witness["multiplicity"]
    constant, bound = witness["linear_constant"], witness["root_bound"]
    rows = []
    for degree in range(dimension):
        power = min(degree, multiplicity)
        shift = degree - power
        row = [0] * dimension
        for j in range(power + 1):
            row[j + shift] = number ** max(multiplicity - degree, 0) * comb(power, j) * constant ** (power - j)
        rows.append(row)
    rebuilt = [sum(witness["lattice_combination"][i] * rows[i][j] for i in range(dimension)) for j in range(dimension)]
    while rebuilt and not rebuilt[-1]:
        rebuilt.pop()
    assert rebuilt == witness["polynomial"]
    assert sum(abs(coefficient) * bound**i for i, coefficient in enumerate(rebuilt)) == witness["weighted_l1_norm"]
    assert witness["weighted_l1_norm"] < (1 << result["factor_exponent"]) ** multiplicity
    assert (constant * (1 << result["precision"]) - result["center"]) % number == 0
    for factor in result["factors"]:
        root = (factor - result["center"]) // (1 << result["precision"])
        assert result["root_interval"][0] <= root < result["root_interval"][1]
        assert abs(root) <= bound
        assert all(sum(coefficient * root**i for i, coefficient in enumerate(row)) % factor**multiplicity == 0 for row in rows)
        assert sum(coefficient * root**i for i, coefficient in enumerate(rebuilt)) == 0


def run_checks():
    hashes = {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES}
    rng = Random(2026110182371)
    counts = {"universal_parameter_inequalities": 0, "literal_window_factor_lists": 0, "literal_rectangle_factor_lists": 0, "smallest_prime_factor_controls": 0, "independent_lattice_certificates": 0, "invalid_inputs_rejected": 0}
    tiny_cases = []
    for bits in list(range(2, 513)) + [1024, 2048, 4096, 16384]:
        dimension = 2 * bits
        required_budget = (dimension ** (2 * dimension)).bit_length()
        exponents = range(1, bits) if bits <= 512 else sorted({bits // 8, bits // 4, bits // 3, bits // 2, 3 * bits // 4, bits - 1})
        for exponent in exponents:
            precision = exponent - exponent * exponent // bits
            assert precision <= (bits + 3) // 4
            radius_exponent = exponent * exponent // bits - 1
            if radius_exponent < 1:
                continue
            multiplicity = 2 * exponent
            budget = (4 * exponent * multiplicity * dimension - dimension * (dimension - 1)
                      - 2 * bits * multiplicity * (multiplicity + 1)
                      - 2 * radius_exponent * dimension * (dimension - 1))
            assert budget >= required_budget
            counts["universal_parameter_inequalities"] += 1
            if bits < 8:
                tiny_cases.append({"bits": bits, "factor_exponent": exponent, "radius_exponent": radius_exponent, "integer_budget": budget, "required_budget": required_budget})
    assert len(tiny_cases) == 7

    negative_root_seen = False
    for _ in range(160):
        number = rng.randrange(3, 5000) | 1
        exponent = rng.randrange(number.bit_length())
        band_left, band_right = 1 << exponent, min(1 << (exponent + 1), number + 1)
        left = rng.randrange(band_left, band_right)
        right = rng.randrange(left, band_right + 1)
        parameters = window_parameters(number, left, right)
        precision = parameters["required_precision"] + rng.randrange(3)
        residue = rng.randrange(1 << precision)
        result = centered_window_factors(number, left, right, residue, precision=precision)
        expected = [value for value in range(left, right) if value % (1 << precision) == residue and number % value == 0]
        assert result["factors"] == expected
        if result["lattice_queries"]:
            _certificate(number, result)
            counts["independent_lattice_certificates"] += 1
            negative_root_seen |= any(value < result["center"] for value in expected)
        counts["literal_window_factor_lists"] += 1

    for _ in range(80):
        number = rng.randrange(3, 600) | 1
        first = sorted(rng.randrange(1, number + 30) for _ in range(2))
        second = sorted(rng.randrange(1, number + 30) for _ in range(2))
        result = centered_factor_rectangle(number, first, second, collect_pairs=True)
        expected = [(value, number // value) for value in range(first[0], first[1]) if number % value == 0 and second[0] <= number // value < second[1]]
        assert result["factor_pairs"] == expected
        assert result["count"] == len(expected)
        assert result["smallest_row"] == (expected[0][0] if expected else None)
        counts["literal_rectangle_factor_lists"] += 1

    for number in range(2, 1001):
        result = smallest_factor_by_centered_windows(number)
        assert result["factor"] == _smallest(number)
        counts["smallest_prime_factor_controls"] += 1

    local_controls = []
    for exponent_three, exponent_five, width_exponent in ((14, 18, None), (80, 55, 32), (320, 220, 150)):
        factor = 3**exponent_three
        number = factor * 5**exponent_five
        exponent = factor.bit_length() - 1
        if width_exponent is None:
            left, right = 1 << exponent, 1 << (exponent + 1)
        else:
            left = factor - (1 << (width_exponent - 2))
            right = left + (1 << width_exponent)
        precision = window_parameters(number, left, right)["required_precision"]
        residue = factor % (1 << precision)
        result = centered_window_factors(number, left, right, residue)
        expected = sorted(value for value in _known_divisors(exponent_three, exponent_five) if left <= value < right and value % (1 << precision) == residue)
        assert result["factors"] == expected and factor in expected
        assert result["lattice_queries"] == 1
        _certificate(number, result)
        counts["independent_lattice_certificates"] += 1
        negative_root_seen |= any(value < result["center"] for value in expected)
        old_precision = (number.bit_length() + 3) // 4 + 2
        if precision:
            distinct_old_classes = min(result["candidate_count"], 1 << (old_precision - precision))
        else:
            distinct_old_classes = min(right // 2 - left // 2, 1 << (old_precision - 1))
        local_controls.append({
            "number": number, "input_bits": number.bit_length(),
            "known_factorization": {"3": exponent_three, "5": exponent_five},
            "window": [left, right], "old_quarter_precision": old_precision,
            "old_nonempty_odd_residue_classes_combined": distinct_old_classes,
            "result": result,
        })
    assert negative_root_seen

    complete_examples = []
    for number in (257 * 263, 269 * 1009):
        result = smallest_factor_by_centered_windows(number)
        assert result["factor"] == _smallest(number)
        assert result["lattice_queries"] > 0
        complete_examples.append({"number": number, "input_bits": number.bit_length(), **result})

    for number, left, right, precision in ((4, 1, 2, None), (15, 0, 1, None), (15, 3, 5, None), (15, 8, 17, None), (67591, 128, 256, 0)):
        try:
            centered_window_factors(number, left, right, precision=precision)
        except ValueError:
            counts["invalid_inputs_rejected"] += 1
        else:
            raise AssertionError((number, left, right, precision))
    assert hashes == {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES}
    return {
        "status": "passed", "seed": 2026110182371, "counts": counts,
        "tiny_parameter_cases": tiny_cases, "local_controls": local_controls,
        "complete_examples": complete_examples, "source_sha256": hashes,
        "scope": "Classical centered small-root queries. Full factoring remains 2**(n/4)*poly(n); large controls are local windows with independently supplied factorizations.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_checks()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], **result["counts"], "local_bits": [item["input_bits"] for item in result["local_controls"]], "complete_examples": result["complete_examples"]}))


if __name__ == "__main__":
    main()
