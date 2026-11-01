#!/usr/bin/env python3
"""Exact affine-flat supports, corrections, and independent profile checks."""

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import random

from check_modular_recombination import CROSSING_SYMMETRIES, crossing, encode, ranks
from hamming_pfaffian import profile
from hamming_restrictions import (
    affine_histogram,
    affine_profile,
    histogram_profile,
    word_histogram,
)
from hamming_tensor_power import local_terms

PAIRS = tuple(combinations(range(4), 2))


def equal_bit(pattern, i, j):
    return ((pattern >> (3 - i)) & 1) == ((pattern >> (3 - j)) & 1)


def hyperplanes():
    return {
        frozenset(p for p in range(16) if (p & normal).bit_count() % 2 == value)
        for normal in range(1, 16)
        for value in (0, 1)
    }


def generic_pair(pair):
    values = [0] * 4
    values[pair[0]] = values[pair[1]] = 2
    remaining = [i for i in range(4) if i not in pair]
    values[remaining[0]], values[remaining[1]] = 4, 8
    return tuple(values)


def permute_pattern(pattern, permutation):
    return sum(
        ((pattern >> (3 - old)) & 1) << (3 - j) for j, old in enumerate(permutation)
    )


def support_checks():
    planes, found = hyperplanes(), []
    assert len(planes) == 30
    for j, (local, scalar) in enumerate(local_terms()):
        assert scalar and max(v.bit_length() - 1 for v in local) <= 4
        assert scalar.bit_length() - 1 + local[0].bit_length() - 1 <= 12
        support = {i for i, value in enumerate(local) if value}
        assert len(support) <= 8
        found += [(j, sorted(flat)) for flat in planes if flat <= support]
    assert found == [
        (0, [0, 2, 5, 7, 8, 10, 13, 15]),
        (5, [0, 2, 5, 7, 8, 10, 13, 15]),
        (13, [0, 3, 4, 7, 8, 11, 12, 15]),
        (14, [0, 3, 4, 7, 8, 11, 12, 15]),
    ]
    checked = 0
    nonzero = 0
    for r in range(3, 9):
        for permutation in CROSSING_SYMMETRIES:
            generic = [
                affine_profile(r, tuple(generic_pair(pair)[i] for i in permutation))
                for pair in PAIRS
            ]
            assert all(g.bit_count() % 2 == 0 for g in generic)
            nonzero += sum(bool(g) for g in generic)
            for flat in planes:
                changed = {permute_pattern(p, permutation) for p in flat}
                histogram = tuple(
                    (1 << (r - 3)) if p in changed else 0 for p in range(16)
                )
                actual = histogram_profile(histogram)
                predicted = 0
                for pair, value in zip(PAIRS, generic):
                    if all(equal_bit(p, *pair) for p in flat):
                        predicted ^= value
                assert actual == predicted
                checked += 1
        if r >= 4:
            assert histogram_profile((1 << (r - 4),) * 16) == 0
    assert nonzero
    return found, checked


def source_checks():
    rng = random.Random(20261030)
    checked = 0
    for r in range(1, 6):
        m = 1 << r
        for _ in range(40):
            indices = tuple(rng.randrange(1 << (r + 1)) for _ in range(4))
            words = tuple(
                sum(
                    (((index >> 1) & point).bit_count() % 2 ^ (index & 1)) << point
                    for point in range(m)
                )
                for index in indices
            )
            assert affine_histogram(r, indices) == word_histogram(m, words)
            h = 1 << (m + 2)
            labels = [(2 * word + 1, h - 2 * word - 1) for word in words]
            x, y, u, v = labels
            assert affine_profile(r, indices) == profile(x, v, u, y)
            checked += 1
    return checked


def binary_rank(values):
    basis = {}
    for value in values:
        while value:
            p = value.bit_length() - 1
            if p in basis:
                value ^= basis[p]
            else:
                basis[p] = value
                break
    return len(basis)


def complete_code_checks():
    r, d = 3, 16
    generic = [affine_profile(r, generic_pair(pair)) for pair in PAIRS]
    equality, forbidden, tensors = [0] * 6, 0, {}
    frames = {}
    for a in range(1, 1 << r):
        for b in range(1, 1 << r):
            if a != b:
                frames.setdefault(tuple(sorted((a, b, a ^ b))), (a, b))
    assert len(frames) == (2**r - 1) * (2**r - 2) // 6
    for x, y, u, v in product(range(d), repeat=4):
        indices = x, y, u, v
        bit = 1 << ((x * d + y) * d * d + u * d + v)
        if binary_rank([i >> 1 for i in indices]) >= 3:
            forbidden |= bit
        for j, pair in enumerate(PAIRS):
            if indices[pair[0]] == indices[pair[1]]:
                equality[j] |= bit
        polynomial = affine_profile(r, indices)
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        while polynomial:
            low = polynomial & -polynomial
            tensors[low.bit_length() - 1] = tensors.get(low.bit_length() - 1, 0) | bit
            polynomial ^= low
    total = 0
    for vector in tensors.values():
        total ^= vector
    assert total == encode(crossing(d), d)
    for degree, vector in tensors.items():
        corrected = vector
        for g, eq in zip(generic, equality):
            if g >> degree & 1:
                corrected ^= eq
        assert corrected & forbidden == 0
        values = set()
        for a, b in frames.values():
            x, y = 2 * a, 2 * b
            pos = (x * d + y) * d * d + y * d + x
            values.add((corrected >> pos) & 1)
        assert len(values) == 1
    for eq in equality:
        mask = (1 << (d * d)) - 1
        rows = [(eq >> (i * d * d)) & mask for i in range(d * d)]
        assert max(ranks(rows, d)) <= d
    # The uncorrected profile really does have rank-three support.
    assert any(vector & forbidden for vector in tensors.values())
    return {
        "crossing_entries": d**4,
        "coefficient_matrices": len(tensors),
        "independent_pair_blocks": len(frames),
        "equality_rank_checks": len(equality),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    found, flat_checks = support_checks()
    direct = source_checks()
    complete = complete_code_checks()
    for r in range(9, 65):
        d = 1 << (r + 1)
        count = ((1 << r) - 1) * ((1 << r) - 2) // 6
        assert count - 6 * d >= d * d // 48
        assert 24 * (count - 6 * d) == d * d - 150 * d + 8
    paths = [Path(__file__), Path(__file__).with_name("hamming_restrictions.py")]
    result = {
        "status": "affine-code construction forces rho>=1; no factoring lower bound",
        "surviving_local_hyperplanes": found,
        "rank_three_symmetry_flat_checks": flat_checks,
        "independent_source_profile_checks": direct,
        "complete_small_code": complete,
        "raw_coefficient_position_bound": 5 * 13 * 65,
        "source_sha256": {p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
