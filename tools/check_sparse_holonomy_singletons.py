#!/usr/bin/env python3
"""Exact singleton Smith profiles and a sparse Gaussian kernel identity."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import check_block_smith_crossings as smith
import check_cancellation_module_automata as field


def pairing_profile(weights, total):
    values = sorted(weights, reverse=True) + [0]
    assert len(set(weights)) == len(weights) and all(weight > 0 for weight in weights)
    result = []
    for index in range(0, len(values) - 1, 2):
        order = total - values[index] - values[index + 1]
        assert order >= 0
        result.extend((order, order))
    if len(values) % 2:
        assert values[-1] == 0
        result.append(total)
    return tuple(sorted(result))


def singleton_matrix(weights, total):
    weights = [0] + list(weights)
    return [[1 << (total - (a if a == b else a + b)) for b in weights] for a in weights]


def unit_recurrence(weights):
    values = sorted(weights, reverse=True)
    unit = field.ONE
    checks = 0
    for index in range(0, len(values) - 1, 2):
        a, c = values[index:index + 2]
        difference = field.radd(unit, field.ONE)
        assert difference == field.ZERO or field.rorder(difference) > a
        da, dc = field.rat(1 | (1 << a)), field.rat(1 | (1 << c))
        denominator = field.radd(field.rmul(da, dc), field.rmul(unit, field.radd(da, dc)))
        assert field.rorder(denominator) == 0
        following = field.rdiv(field.rmul(unit, field.rmul(da, dc)), denominator)
        # An independent 2-by-2 inverse gives the same Schur scalar.
        inverse_sum = field.rdiv(field.radd(da, dc), denominator)
        schur = field.radd(unit, field.rmul(field.rmul(unit, unit), inverse_sum))
        assert following == schur
        assert field.rorder(following) == 0
        assert field.rorder(field.radd(following, field.ONE)) == c
        normalized = ((field.rshift(field.radd(unit, da), c - a), unit),
                      (unit, field.rshift(field.radd(unit, dc), a - c)))
        assert all(entry == field.ZERO or field.rorder(entry) >= 0 for row in normalized for entry in row)
        determinant = field.radd(field.rmul(normalized[0][0], normalized[1][1]), field.rmul(unit, unit))
        assert field.rorder(determinant) == 0
        unit = following
        checks += 1
    return checks


def pfaffian(matrix):
    if not matrix:
        return 1
    assert len(matrix) % 2 == 0
    answer = 0
    for partner in range(1, len(matrix)):
        indices = [j for j in range(1, len(matrix)) if j != partner]
        minor = [[matrix[a][b] for b in indices] for a in indices]
        answer ^= field.pmul(matrix[0][partner], pfaffian(minor))
    return answer


def gaussian_check(bits, radius):
    weights = [0] + [1 << j for j in range(bits)]
    modes = len(weights)
    kernel = [[int(i != j) for j in range(2 * modes)] for i in range(2 * modes)]
    for i, weight in enumerate(weights):
        kernel[i][modes + i] = kernel[modes + i][i] = 1 << weight
    masks = [a for a in range(1 << bits) if a.bit_count() <= radius]
    encoded = {}
    for mask in masks:
        selected = [j + 1 for j in range(bits) if mask >> j & 1]
        if len(selected) % 2 != radius % 2:
            selected = [0] + selected
        assert len(selected) <= radius and len(selected) % 2 == radius % 2
        encoded[mask] = selected
    assert len({tuple(indices) for indices in encoded.values()}) == len(masks)
    expected_count = sum(1 for degree in range(radius + 1) if degree % 2 == radius % 2
                         for _ in combinations(range(modes), degree))
    assert expected_count == len(masks)
    checks = 0
    for a in masks:
        for b in masks:
            indices = encoded[a] + [modes + j for j in encoded[b]]
            matrix = [[kernel[i][j] for j in indices] for i in indices]
            value = pfaffian(matrix)
            assert value == 1 << (a & b)
            assert field.rshift(field.rat(value), (1 << bits) - 1 - a - b) == field.rat(
                1 << ((1 << bits) - 1 - (a | b)))
            checks += 1
    return checks


def check(maximum_bits=8):
    records, toeplitz, recurrence_checks = [], 0, 0
    weight_sets = [tuple(1 << j for j in range(bits)) for bits in range(1, maximum_bits + 1)]
    weight_sets.extend(((1, 3, 7, 11), (2, 5, 9, 14, 20), (1, 4, 6, 13, 17, 23)))
    for weights in weight_sets:
        total = sum(weights)
        matrix = singleton_matrix(weights, total)
        actual, precision = smith.smith_profile(matrix, total)
        predicted = pairing_profile(weights, total)
        assert actual == predicted
        assert sum(actual) == (len(weights) + 1) * total - 2 * sum(weights)
        for p in sorted(set((1, 2, 3, 4, *(order + 1 for order in actual)))):
            assert smith.toeplitz_rank(matrix, p) == sum(max(0, p - order) for order in actual)
            toeplitz += 1
        recurrence_checks += unit_recurrence(weights)
        records.append({"weights": list(weights), "total": total, "profile": dict(sorted(Counter(actual).items())),
                        "orders": list(actual), "precision": precision})
    gaussian = sum(gaussian_check(bits, radius) for bits, radius in ((2, 1), (3, 2), (4, 2), (4, 3), (5, 3)))
    return {"status": "pass", "field": 2,
            "scope": "All-size weighted singleton Smith theorem and Gaussian encoding; no higher-radius profile theorem or improved factoring exponent.",
            "records": records, "independent_toeplitz_ranks": toeplitz,
            "exact_unit_recurrence_steps": recurrence_checks, "literal_gaussian_entries": gaussian,
            "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "smith_checker_sha256": sha256(Path(smith.__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-bits", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = check(args.maximum_bits)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
