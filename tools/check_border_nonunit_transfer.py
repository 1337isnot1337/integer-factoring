#!/usr/bin/env python3
"""Finite transcription checks for border extraction and chain-ring budgets.

This does not search for an improved seed. Dense finite matrices here
test identities and rank estimates; they do not implement the final
asymptotic machine or measure its memory usage.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import product
import json
from math import factorial
from pathlib import Path
import random


def field_rank(matrix: list[list[int]], p: int) -> int:
    a = [[entry % p for entry in row] for row in matrix]
    if not a:
        return 0
    r = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inverse = pow(a[r][col], -1, p)
        a[r] = [(entry * inverse) % p for entry in a[r]]
        for i in range(r + 1, len(a)):
            scale = a[i][col]
            if scale:
                a[i] = [(x - scale * y) % p for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def chain_inner_rank(matrix: list[list[int]], p: int, h: int) -> int:
    """Number of nonzero Smith entries over Z/p^h, via valuation pivots."""
    modulus = p**h
    a = [[entry % modulus for entry in row] for row in matrix]
    rows = len(a)
    cols = len(a[0])

    def valuation(entry: int) -> int:
        if not entry:
            return h
        v = 0
        while entry % p == 0:
            entry //= p
            v += 1
        return v

    r = 0
    while r < min(rows, cols):
        v, i, j = min((valuation(a[i][j]), i, j)
                      for i in range(r, rows) for j in range(r, cols))
        if v == h:
            break
        a[r], a[i] = a[i], a[r]
        for row in a:
            row[r], row[j] = row[j], row[r]
        pivot_power = p**v
        inverse_unit = pow(a[r][r] // pivot_power, -1, modulus)
        a[r] = [(entry * inverse_unit) % modulus for entry in a[r]]
        assert a[r][r] == pivot_power
        for i in range(r + 1, rows):
            assert a[i][r] % pivot_power == 0
            scale = a[i][r] // pivot_power
            a[i] = [(x - scale * y) % modulus for x, y in zip(a[i], a[r])]
        for j in range(r + 1, cols):
            assert a[r][j] % pivot_power == 0
            scale = a[r][j] // pivot_power
            for i in range(rows):
                a[i][j] = (a[i][j] - scale * a[i][r]) % modulus
        r += 1
    return r


def coefficient_product(polynomials: list[list[int]], degree: int, p: int) -> int:
    coefficients = [1] + [0] * degree
    for polynomial in polynomials:
        new = [0] * (degree + 1)
        for i, x in enumerate(coefficients):
            if x:
                for j, y in enumerate(polynomial[:degree + 1 - i]):
                    new[i + j] = (new[i + j] + x * y) % p
        coefficients = new
    return coefficients[degree]


def tensor_matrix(entries: dict[tuple[int, int, int, int], int], d: int,
                  cut: str) -> list[list[int]]:
    result = [[0] * (d*d) for _ in range(d*d)]
    for (x, y, u, v), value in entries.items():
        if cut == "row":
            result[x*d+y][u*d+v] = value
        else:
            result[x*d+u][y*d+v] = value
    return result


def border_case(p: int, k: int) -> dict[str, int]:
    # A_0=fJ, A_1=zE-fJ-z^2K, A_2=z^2K, with f=1+z+z^2.
    # J and K are fully separable; A_1 has the safe full-rank bound 4.
    # Their sum is zE. Constants cancel outside the crossing support.
    d = 2
    channels = [1, d*d, 1]
    local = []
    for i in range(3):
        tensor = {}
        for indices in product(range(d), repeat=4):
            x, y, u, v = indices
            crossing = int(x == v and y == u)
            special = int(x == 0 and v == 1)
            if i == 0:
                polynomial = [1, 1, 1]
            elif i == 1:
                polynomial = [-1, crossing - 1, -1 - special]
            else:
                polynomial = [0, 0, special]
            tensor[indices] = [c % p for c in polynomial]
        local.append(tensor)

    size = d**k
    digits = [tuple((number >> (k - 1 - j)) & 1 for j in range(k))
              for number in range(size)]
    pieces: dict[tuple[int, ...], dict[tuple[int, int, int, int], int]] = defaultdict(dict)
    for word in product(range(3), repeat=k):
        counts = tuple(word.count(i) for i in range(3))
        output = pieces[counts]
        for indices in product(range(size), repeat=4):
            polynomial_entries = [local[word[j]][tuple(digits[a][j] for a in indices)]
                                  for j in range(k)]
            coefficient = coefficient_product(polynomial_entries, k, p)
            if coefficient:
                output[indices] = (output.get(indices, 0) + coefficient) % p

    for indices in product(range(size), repeat=4):
        x, y, u, v = indices
        actual = sum(piece.get(indices, 0) for piece in pieces.values()) % p
        assert actual == int(x == v and y == u), (p, k, indices)

    nonzero_pieces = 0
    for counts, entries in pieces.items():
        denominator = 1
        channel_product = 1
        for count, channels_i in zip(counts, channels):
            denominator *= factorial(count)
            channel_product *= channels_i**count
        multinomial = factorial(k) // denominator
        limit = (k + 1) * multinomial * channel_product
        row_rank = field_rank(tensor_matrix(entries, size, "row"), p)
        col_rank = field_rank(tensor_matrix(entries, size, "col"), p)
        assert row_rank <= limit and col_rank <= limit
        nonzero_pieces += row_rank > 0
    return {"p": p, "k": k, "labels": size, "types": len(pieces),
            "nonzero_pieces": nonzero_pieces, "required_precision": k + 1}


def chain_budget_checks() -> int:
    rng = random.Random(20261030)
    checked = 0
    for p, h, d in product([2, 3], [2, 3], [2, 3]):
        modulus = p**h
        for t in range(h):
            for c in range(3):
                for _ in range(12):
                    factors = [([[rng.randrange(modulus) for _ in range(d)] for _ in range(d)],
                                [[rng.randrange(modulus) for _ in range(d)] for _ in range(d)])
                               for _ in range(c)]
                    a = [[0] * (d*d) for _ in range(d*d)]
                    for x, y, u, v in product(range(d), repeat=4):
                        b = sum(left[x][u] * right[y][v] for left, right in factors)
                        perturbation = rng.randrange(modulus)
                        target = p**t * (int(x == v and y == u) + p * perturbation)
                        a[x*d+y][u*d+v] = (target - b) % modulus
                    r = chain_inner_rank(a, p, h)
                    for s in range(1, d + 1):
                        assert r + s*s*c >= s*d, (p, h, d, t, c, r, s)
                        checked += 1
    return checked


def digit_rank_checks() -> list[dict[str, int]]:
    results = []
    for p, e in [(2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (5, 1)]:
        size = p**e
        matrix = [[a+b for b in range(size)] for a in range(size)]
        before = chain_inner_rank(matrix, p, e + 1)
        after = field_rank([[(entry // size) % p for entry in row] for row in matrix], p)
        assert before <= 2 and after == size - 1
        results.append({"p": p, "precision_digit": e, "labels": size,
                        "ring_channels": before, "digit_rank": after})
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "status": "passed",
        "scope": "Finite identity and rank checks; no improved seed or measured asymptotic space bound.",
        "border_cases": [border_case(p, k) for p, k in
                         [(2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (5, 1), (5, 2)]],
        "chain_weighted_budget_inequalities": chain_budget_checks(),
        "digit_rank_cases": digit_rank_checks(),
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
