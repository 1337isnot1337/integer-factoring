#!/usr/bin/env python3
"""Polynomial-time modular angle factors, with independent exact checks.

Anthropic's original operators are pinned in check_partition_channels.py.
The polynomial-time Pfaffian method follows PR22's modular-pfaffian.md;
the three-block angle matrix is derived in homogeneous-central-modular-symbols.md.
No stronger factoring exponent is asserted.
"""

import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import random

from check_homogeneous_central_symbols import (
    angle_symbol_two, binomial_two, multiply_two, reference_angle_two,
)
from check_modular_crossing import regular_partitions


def divide_exact(numerator: int, denominator: int) -> int:
    if not denominator:
        raise ZeroDivisionError("zero polynomial denominator")
    quotient = 0
    while numerator and numerator.bit_length() >= denominator.bit_length():
        shift = numerator.bit_length() - denominator.bit_length()
        quotient ^= 1 << shift
        numerator ^= denominator << shift
    if numerator:
        raise ValueError("polynomial quotient is not exact")
    return quotient


def determinant_polynomial(matrix: list[list[int]]) -> int:
    """Bareiss elimination over F2[Z], with row pivots and exact divisions."""
    n = len(matrix)
    if not n:
        return 1
    if any(len(row) != n for row in matrix):
        raise ValueError("square matrix required")
    rows = [row[:] for row in matrix]
    previous = 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if rows[i][k]), None)
        if pivot_row is None:
            return 0
        rows[k], rows[pivot_row] = rows[pivot_row], rows[k]
        pivot = rows[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = multiply_two(pivot, rows[i][j]) ^ multiply_two(rows[i][k], rows[k][j])
                rows[i][j] = divide_exact(numerator, previous)
            rows[i][k] = 0
        previous = pivot
    return rows[-1][-1]


def pfaffian_polynomial(matrix: list[list[int]]) -> int:
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("square matrix required")
    if any(matrix[i][i] for i in range(n)) or any(
        matrix[i][j] != matrix[j][i] for i in range(n) for j in range(i)
    ):
        raise ValueError("an alternating characteristic-two matrix is required")
    if n % 2:
        return 0
    determinant = determinant_polynomial(matrix)
    answer = 0
    while determinant:
        degree = (determinant & -determinant).bit_length() - 1
        if degree % 2:
            raise AssertionError("alternating determinant is not a polynomial square")
        answer ^= 1 << (degree // 2)
        determinant &= determinant - 1
    return answer


def angle_matrix(parts: tuple[int, ...], output: tuple[int, ...], input_: tuple[int, ...]) -> list[list[int]]:
    """Blocks O, U, Y; only O-state and equal U-Y pairs are allowed."""
    no, nu, ny = len(parts), len(output), len(input_)
    n = no + nu + ny
    matrix = [[0] * n for _ in range(n)]
    for i, m in enumerate(parts):
        a = (m - 1) // 2
        for j, k in enumerate(output):
            exponent = 2 * (a & ((k - 1) // 2))
            matrix[i][no + j] = matrix[no + j][i] = 1 << exponent
        for j, k in enumerate(input_):
            exponent = 2 * (a & ~((k - 1) // 2))
            matrix[i][no + nu + j] = matrix[no + nu + j][i] = 1 << exponent
    for i, u in enumerate(output):
        for j, y in enumerate(input_):
            matrix[no + i][no + nu + j] = matrix[no + nu + j][no + i] = int(u == y)
    return matrix


def angle_pfaffian(parts: tuple[int, ...], output: tuple[int, ...], input_: tuple[int, ...]) -> int:
    for label in (parts, output, input_):
        if len(set(label)) != len(label) or any(k <= 0 or k % 2 == 0 for k in label):
            raise ValueError("distinct positive odd parts required")
    polynomial = pfaffian_polynomial(angle_matrix(parts, output, input_))
    return multiply_two(binomial_two(sum(parts) + len(parts)), polynomial)


def matching_pfaffian(matrix: list[list[int]]) -> int:
    """Independent perfect-matching recursion, only for the finite checker."""
    @lru_cache(None)
    def visit(mask):
        if not mask:
            return 1
        first = mask & -mask
        i = first.bit_length() - 1
        rest = mask ^ first
        answer = 0
        choices = rest
        while choices:
            bit = choices & -choices
            j = bit.bit_length() - 1
            answer ^= multiply_two(matrix[i][j], visit(rest ^ bit))
            choices ^= bit
        return answer
    return visit((1 << len(matrix)) - 1) if len(matrix) % 2 == 0 else 0


def verify() -> dict:
    rng = random.Random(40120261030)
    matrix_cases = 0
    for n in range(0, 13, 2):
        for trial in range(10):
            matrix = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i):
                    value = rng.getrandbits(5) if trial < 8 else int(trial == 9)
                    matrix[i][j] = matrix[j][i] = value
            assert pfaffian_polynomial(matrix) == matching_pfaffian(matrix)
            matrix_cases += 1

    triples = []
    for h in range(0, 13):
        for parts in regular_partitions(2, h):
            for ei in range(0, 15):
                inputs = regular_partitions(2, ei)
                for eo in range(0, 15):
                    outputs = regular_partitions(2, eo)
                    if inputs and outputs:
                        triples.append((parts, rng.choice(outputs), rng.choice(inputs)))
    rng.shuffle(triples)
    triples = triples[:1200]
    reference_cases = 0
    for parts, output, input_ in triples:
        actual = angle_pfaffian(parts, output, input_)
        assert actual == angle_symbol_two(parts, output, input_)
        assert actual.bit_length() <= 2 * sum(parts) + 1
        if sum(parts) <= 8:
            assert actual == reference_angle_two(parts, output, input_)
            reference_cases += 1

    division_cases = 0
    for _ in range(300):
        a, b = rng.getrandbits(35), rng.getrandbits(20) or 1
        assert divide_exact(multiply_two(a, b), b) == a
        division_cases += 1
    try:
        divide_exact(1, 2)
        raise AssertionError("nonexact division was accepted")
    except ValueError:
        pass
    try:
        pfaffian_polynomial([[1, 0], [0, 0]])
        raise AssertionError("nonalternating matrix was accepted")
    except ValueError:
        pass
    # Removing the U-Y block incorrectly drops an untouched state.
    assert angle_pfaffian((), (1,), (1,)) == 1
    broken = angle_matrix((), (1,), (1,))
    broken[0][1] = broken[1][0] = 0
    assert pfaffian_polynomial(broken) == 0
    return {
        "scope": "Exact polynomial angle-factor construction; no improved factoring exponent",
        "alternating_matrix_cases": matrix_cases,
        "angle_triples": len(triples),
        "original_integral_reference_triples": reference_cases,
        "exact_polynomial_divisions": division_cases,
        "negative_controls": ["nonexact division", "nonalternating matrix", "omitted unchanged-state edge"],
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    args = parser.parse_args()
    evidence = verify()
    text = json.dumps(evidence, sort_keys=True, indent=2) + '\n'
    if args.output:
        Path(args.output).write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
