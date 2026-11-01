#!/usr/bin/env python3
"""Independent finite type-amplification checks for the gap-free theorem."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import factorial, isqrt
from pathlib import Path
import random


def cross(x, y, u, v):
    return int(x == v and y == u)


def index(d, x, y, u, v):
    return ((x * d + y) * d + u) * d + v


def tensor(d, function, p):
    return [function(x, y, u, v) % p for x, y, u, v in product(range(d), repeat=4)]


def rank(matrix, p):
    if p == 2:
        basis = {}
        for row in matrix:
            word = sum((entry & 1) << i for i, entry in enumerate(row))
            while word:
                pivot = word.bit_length() - 1
                if pivot not in basis:
                    basis[pivot] = word
                    break
                word ^= basis[pivot]
        return len(basis)
    rows = [row[:] for row in matrix]
    answer = 0
    columns = len(rows[0]) if rows else 0
    for column in range(columns):
        pivot = next((r for r in range(answer, len(rows)) if rows[r][column] % p), None)
        if pivot is None:
            continue
        rows[answer], rows[pivot] = rows[pivot], rows[answer]
        inverse = pow(rows[answer][column], -1, p)
        rows[answer] = [entry * inverse % p for entry in rows[answer]]
        for r in range(answer + 1, len(rows)):
            multiplier = rows[r][column]
            if multiplier:
                rows[r] = [(x - multiplier * y) % p for x, y in zip(rows[r], rows[answer])]
        answer += 1
    return answer


def cut_ranks(values, d, p):
    row = [[0] * (d * d) for _ in range(d * d)]
    col = [[0] * (d * d) for _ in range(d * d)]
    for x, y, u, v in product(range(d), repeat=4):
        value = values[index(d, x, y, u, v)]
        row[x * d + y][u * d + v] = value
        col[x * d + u][y * d + v] = value
    return rank(row, p), rank(col, p)


def width(value):
    return (value - 1).bit_length() if value else 0


def seed_quality(seed, a, p):
    d = 1 << a
    expected = tensor(d, cross, p)
    assert [sum(values[j] for values in seed) % p for j in range(d**4)] == expected
    ranks = [cut_ranks(values, d, p) for values in seed]
    assert all(r and c for r, c in ranks)
    w = max(width(r) + width(c) for r, c in ranks)
    quality = Fraction(w + 2 * width(len(seed)), 4 * a)
    assert quality >= Fraction(1, 2)
    return quality, ranks, w


def check_amplification(seed, a, target, p):
    quality, seed_ranks, w = seed_quality(seed, a, p)
    root_scale = 1 << 40
    # Certified rational lower bound for the possibly irrational sum.
    # Finite examples have ample slack unless the single term is an
    # exact square. No floating-point logarithms or roots are used.
    root_numerators = [isqrt(r * c * root_scale**2) for r, c in seed_ranks]
    rank_sum_numerator = sum(root_numerators)
    k, remainder = divmod(target, a)
    count = len(seed)
    types = sorted({tuple(word.count(i) for i in range(count)) for word in product(range(count), repeat=k)})
    d = 1 << target
    seed_d = 1 << a
    arrays = {nu: [0] * d**4 for nu in types}
    coordinate_checks = 0
    for external in product(range(d), repeat=4):
        rest = tuple(value >> (a * k) for value in external)
        rest_value = cross(*rest)
        blocks = [tuple((value >> (a * t)) & (seed_d - 1) for value in external) for t in range(k)]

        # Coefficient dynamic program.
        states = {(0,) * count: 1}
        for block in blocks:
            next_states = defaultdict(int)
            for nu, coefficient in states.items():
                for piece in range(count):
                    new_type = tuple(nu[i] + (i == piece) for i in range(count))
                    next_states[new_type] += coefficient * seed[piece][index(seed_d, *block)]
            states = {nu: coefficient % p for nu, coefficient in next_states.items()}

        # Independent enumeration of tensor-product piece words.
        direct = defaultdict(int)
        for word in product(range(count), repeat=k):
            coefficient = rest_value
            for piece, block in zip(word, blocks):
                coefficient *= seed[piece][index(seed_d, *block)]
            direct[tuple(word.count(i) for i in range(count))] += coefficient
        for nu in types:
            coefficient = rest_value * states.get(nu, 0) % p
            assert coefficient == direct[nu] % p
            arrays[nu][index(d, *external)] = coefficient
            coordinate_checks += 1
        assert sum(direct.values()) % p == cross(*external)

    nonzero_types = 0
    for nu, values in arrays.items():
        row_rank, col_rank = cut_ranks(values, d, p)
        if not row_rank:
            assert not col_rank
            continue
        nonzero_types += 1
        multinomial = factorial(k)
        for multiplicity in nu:
            multinomial //= factorial(multiplicity)
        row_bound = col_bound = multinomial * (1 << (2 * remainder))
        for multiplicity, (r, c) in zip(nu, seed_ranks):
            row_bound *= r**multiplicity
            col_bound *= c**multiplicity
        assert row_rank <= row_bound and col_rank <= col_bound
        type_weight_square = multinomial**2
        for multiplicity, (r, c) in zip(nu, seed_ranks):
            type_weight_square *= (r * c)**multiplicity
        assert type_weight_square * root_scale**(2 * k) <= rank_sum_numerator**(2 * k)
        assert row_bound * col_bound == type_weight_square * (1 << (4 * remainder))
        assert (1 << (width(row_rank) + width(col_rank))) * root_scale**(2 * k) <= 4 * rank_sum_numerator**(2 * k) * (1 << (4 * remainder))
        assert width(row_rank) <= 2 * target and width(col_rank) <= 2 * target
        exact_bound = k * (w + 2 * width(count)) + 4 * remainder + 2
        assert width(row_rank) + width(col_rank) <= exact_bound
        if quality <= 1:
            assert exact_bound <= 4 * quality * target + 4 * a + 2
    assert len(types) <= (k + 1) ** count
    return coordinate_checks, nonzero_types


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    e = tensor(2, cross, 2)
    split = [tensor(2, lambda x, y, u, v, b=b: cross(x, y, u, v) * (x == b), 2) for b in range(2)]
    r = tensor(2, lambda x, y, u, v: (x == y) * (u == v), 2)
    c = tensor(2, lambda x, y, u, v: (x == u) * (y == v), 2)
    h = [(x + y + z) % 2 for x, y, z in zip(e, r, c)]
    rng = random.Random(20261030)
    random_piece = [rng.randrange(3) for _ in range(16)]
    e3 = tensor(2, cross, 3)
    random_split = [random_piece, [(x - y) % 3 for x, y in zip(e3, random_piece)]]
    cases = [
        ([e], 1, 3, 2),
        (split, 1, 4, 2),
        ([r, c, h], 1, 3, 2),
        (random_split, 1, 3, 3),
        ([tensor(4, cross, 2)], 2, 3, 2),
    ]
    coordinates = nonzero_types = 0
    for seed, a, target, p in cases:
        n, m = check_amplification(seed, a, target, p)
        coordinates += n
        nonzero_types += m

    entropy_checks = 0
    for target in range(1, 65537):
        cap = isqrt(target)
        # ceil(log2(lambda+1)) overestimates the actual type entropy;
        # even this integer envelope is <= 2lambda.
        assert cap * target.bit_length() <= 2 * target
        assert cap * cap <= target < (cap + 1) ** 2
        entropy_checks += 1
    result = {
        "status": "passed",
        "scope": "Finite type identities and resource envelopes; no new strong seed or factoring exponent.",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "fields": [2, 3],
        "seed_amplification_cases": len(cases),
        "independent_word_vs_coefficient_entries": coordinates,
        "nonzero_type_rank_pairs_checked": nonzero_types,
        "exact_rank_sum_type_certificates": nonzero_types,
        "exact_integer_entropy_envelopes": entropy_checks,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
