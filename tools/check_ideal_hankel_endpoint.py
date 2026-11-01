#!/usr/bin/env python3
"""Exact binary checks for the normalized ideal-Hankel endpoint theorem."""

import argparse
import json
from pathlib import Path
from random import Random


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            bit = row.bit_length() - 1
            if bit in pivots:
                row ^= pivots[bit]
            else:
                pivots[bit] = row
                break
    return len(pivots)


def identity(n):
    return tuple(1 << i for i in range(n))


def multiply(left, right):
    result = []
    for row in left:
        value = 0
        while row:
            bit = row & -row
            row ^= bit
            value ^= right[bit.bit_length() - 1]
        result.append(value)
    return tuple(result)


def inverse(matrix):
    n = len(matrix)
    rows = [row | (1 << (n + i)) for i, row in enumerate(matrix)]
    for column in range(n):
        pivot = next((i for i in range(column, n) if rows[i] >> column & 1), None)
        if pivot is None:
            raise ValueError("singular matrix")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for i in range(n):
            if i != column and rows[i] >> column & 1:
                rows[i] ^= rows[column]
    assert tuple(row & ((1 << n) - 1) for row in rows) == identity(n)
    return tuple(row >> n for row in rows)


def tensor(left, right):
    n = len(right)
    result = []
    for a in left:
        for b in right:
            value = 0
            row = a
            while row:
                bit = row & -row
                row ^= bit
                value ^= b << (n * (bit.bit_length() - 1))
            result.append(value)
    return tuple(result)


def swap(d):
    return tuple(1 << (y * d + x) for x in range(d) for y in range(d))


def opposite(matrix, d):
    result = [0] * (d * d)
    for row_index, row in enumerate(matrix):
        x, y = divmod(row_index, d)
        while row:
            bit = row & -row
            row ^= bit
            v, u = divmod(bit.bit_length() - 1, d)
            result[x * d + u] ^= 1 << (v * d + y)
    return tuple(result)


def conjugate(matrix, basis):
    return multiply(multiply(basis, matrix), inverse(basis))


def hankel_coefficients(n):
    return tuple(tuple(1 << (degree - i) if 0 <= degree - i < n else 0
                       for i in range(n)) for degree in range(n))


def smith_multiplicities(coefficients, maximum):
    n = len(coefficients[0])
    previous_rank = previous_count = 0
    result = []
    for order in range(maximum + 1):
        rows = []
        for block in range(order + 1):
            for row in range(n):
                value = 0
                for earlier in range(max(0, block - len(coefficients) + 1), block + 1):
                    value ^= coefficients[block - earlier][row] << (n * earlier)
                rows.append(value)
        current_rank = rank(rows)
        current_count = current_rank - previous_rank
        result.append(current_count - previous_count)
        previous_rank, previous_count = current_rank, current_count
        if current_count == n:
            return result
    raise AssertionError("valuation cutoff was too small")


def random_basis(n, rng):
    rows = list(identity(n))
    for _ in range(8 * n):
        i, j = rng.sample(range(n), 2)
        rows[i] ^= rows[j]
    return tuple(rows)


def verify():
    rng = Random(202610300742)
    normalization_cases = smith_cases = twisted_swap_cases = 0
    for n in range(2, 26):
        h1 = tuple((1 << (n - i)) - 1 for i in range(n))
        reversal = tuple(1 << (n - i - 1) for i in range(n))
        jordan = tuple((1 << i) | ((1 << (i + 1)) if i + 1 < n else 0) for i in range(n))
        assert multiply(reversal, inverse(h1)) == jordan
        for _ in range(3):
            basis = random_basis(n, rng)
            right = multiply(inverse(h1), inverse(basis))
            assert multiply(multiply(basis, h1), right) == identity(n)
            assert multiply(multiply(basis, reversal), right) == conjugate(jordan, basis)
            normalization_cases += 1
        if n in (4, 9, 16, 25):
            valuations = smith_multiplicities(hankel_coefficients(n), n)
            assert valuations == [1] + [0] * (n - 1) + [n - 1]
            smith_cases += 1

    for d in range(2, 9):
        for _ in range(40):
            p = random_basis(d, rng)
            unipotent = tuple((1 << i) | sum(rng.randrange(2) << j for j in range(i + 1, d))
                              for i in range(d))
            q = multiply(inverse(p), unipotent)
            leading = multiply(tensor(p, q), swap(d))
            assert rank(opposite(leading, d)) == 1
            assert multiply(leading, leading) == tensor(multiply(p, q), multiply(q, p))
            powered = leading
            for _ in range(1 + (d - 1).bit_length()):
                powered = multiply(powered, powered)
            assert powered == identity(d * d)
            twisted_swap_cases += 1

    d = 2
    q = (3, 2)
    leading = multiply(tensor(identity(d), q), swap(d))
    nilpotent = tuple(a ^ b for a, b in zip(leading, identity(4)))
    current = identity(4)
    for exponent in range(5):
        assert rank(current) == 4 - exponent
        current = multiply(current, nilpotent)
    basis_columns = (1, 7, 4, 8)
    basis = tuple(sum(((column >> i) & 1) << j for j, column in enumerate(basis_columns))
                  for i in range(4))
    h1 = tuple((1 << (4 - i)) - 1 for i in range(4))
    reversal = tuple(1 << (3 - i) for i in range(4))
    assert conjugate(multiply(reversal, inverse(h1)), basis) == leading
    assert rank(opposite(leading, 2)) == 1

    dimensions = 0
    for d in range(3, 10001):
        regular_order = 1 << (d * d - 1).bit_length()
        twisted_upper = 1 << (1 + (d - 1).bit_length())
        assert regular_order > twisted_upper
        dimensions += 1
    return {
        "status": "pass",
        "constant_normalization_cases": normalization_cases,
        "exact_smith_profiles": smith_cases,
        "twisted_swap_order_cases": twisted_swap_cases,
        "order_separation_dimensions": dimensions,
        "two_label_endpoint": "rank-one opposite cut and one regular unipotent block of size four",
        "scope": "The specified ideal Hankel curve under constant normalization, not general Smith profiles.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    evidence = verify()
    rendered = json.dumps(evidence, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
