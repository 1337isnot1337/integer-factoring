#!/usr/bin/env python3
"""Exact checks of the PBW-length degeneration in both crossing cuts."""

from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from check_noncommutative_crossing_hopf import (
    add, coefficients, column_ranks, compose, maps, multiply,
    parameter_translate, scale, variable,
)
from check_noncommutative_crossing_orbit import binary_rank
from check_noncommutative_crossing_pencils import generator_maps
from check_noncommutative_crossing_ut import RestrictedUpperTriangular


@lru_cache(None)
def abelian_matrices(weights):
    width = len(weights)
    images = []
    for side in range(2):
        for j, weight in enumerate(weights):
            images.append({1 << (side * width + j): 1 << weight,
                           1 << ((1 - side) * width + j): 1 ^ (1 << weight)})
    return coefficients([parameter_translate(image) for image in images])


def both_ranks(matrix, width):
    return binary_rank(matrix), column_ranks([matrix], width)[0]


def degeneration_check(original, abelian, width):
    dimension = 1 << (2 * width)
    count = max(len(original), len(abelian))
    zero = [0] * dimension
    positive_order_entries = 0
    for degree in range(count):
        source = original[degree] if degree < len(original) else zero
        target = abelian[degree] if degree < len(abelian) else zero
        leading = []
        for output, row in enumerate(source):
            diagonal = 0
            output_length = output.bit_count()
            while row:
                bit = row & -row
                input = bit.bit_length() - 1
                valuation = input.bit_count() - output_length
                assert valuation >= 0
                if valuation == 0:
                    diagonal ^= bit
                else:
                    positive_order_entries += 1
                row ^= bit
            leading.append(diagonal)
        assert leading == target
    return positive_order_entries


def combine(matrices, mask):
    result = [0] * len(matrices[0])
    for degree, matrix in enumerate(matrices):
        if mask >> degree & 1:
            result = [a ^ b for a, b in zip(result, matrix)]
    return result


def check(name, weights, images, exhaustive=False):
    original = coefficients([parameter_translate(image) for image in images])
    abelian = abelian_matrices(tuple(weights))
    width = len(weights)
    positive = degeneration_check(original, abelian, width)
    coefficient_records = []
    for i, matrix in enumerate(original):
        expected = abelian[i] if i < len(abelian) else [0] * len(matrix)
        r, c = both_ranks(matrix, width)
        ar, ac = both_ranks(expected, width)
        assert r >= ar and c >= ac
        coefficient_records.append([i, r, c, ar, ac])
    count = len(original)
    masks = range(1, 1 << count) if exhaustive else sorted(
        {1, (1 << count) - 1, (1 << (count - 1)) | 1, 48 & ((1 << count) - 1)}
        | {((1 << i) | (1 << (i + 1))) for i in range(count - 1)}
        | {sum(1 << i for i in range(parity, count, 2)) for parity in (0, 1)})
    recombinations = []
    for mask in masks:
        source, target = combine(original, mask), combine(abelian, mask)
        r, c = both_ranks(source, width)
        ar, ac = both_ranks(target, width)
        assert r >= ar and c >= ac
        recombinations.append([mask, r, c, ar, ac])
    return {"name": name, "label_bits": width, "labels": 1 << width,
            "primitive_weights": list(weights),
            "coefficientwise_leading_tensor_equals_abelian": True,
            "positive_PBW_valuation_entries": positive,
            "coefficient_ranks_and_abelian_ranks": coefficient_records,
            "recombination_ranks_and_abelian_ranks": recombinations,
            "all_binary_recombinations": exhaustive}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    identity, one, _, _ = maps(1, 1, (1,))
    _, two, _, _ = maps(2, 1, (1, 1))
    field, _ = generator_maps(4, 2, ((0, 2, 1), (0, 3, 2), (1, 2, 2), (1, 3, 3)))
    jordan, _ = generator_maps(4, 2, ((0, 2, 1), (1, 3, 1), (0, 3, 2)))
    gauge = list(identity)
    gauge[2] = add(gauge[2], scale(multiply(variable(0), variable(1)), 2))
    gauge[5] = add(gauge[5], scale(multiply(variable(3), variable(4)), 4))
    gauged = compose(gauge, compose(one, gauge))
    upper = RestrictedUpperTriangular(4)
    cases = [
        check("one_Heisenberg_pair", (1, 1, 2), one, True),
        check("two_pairs_shared_center", (1, 1, 1, 1, 2), two),
        check("irreducible_F4_pencil", (1, 1, 1, 1, 2, 2), field),
        check("nilpotent_Jordan_pencil", (1, 1, 1, 1, 2, 2), jordan),
        check("parameter_dependent_local_gauge", (1, 1, 2), gauged),
        check("restricted_UT4", upper.weights, upper.generator_images()),
    ]
    result = {
        "status": "all coefficient leading tensors and both-cut rank domination checks pass",
        "scope": "finite checks supplement the ordinary all-size filtered proof",
        "cases": cases,
        "verified_remote_sources": [
            {"commit": "b49b9e19480c65d139ef15beecacb92064c1a960",
             "path": "research/full-conservation-rank-sum-barrier.md",
             "blob_sha": "1b3e426458157f51a3833476d24e9d4bb3d53803"},
            {"commit": "023e4fc4482d72f25738f6b1d0f9cd31583baed1",
             "path": "research/modular-binary-word-barrier.md",
             "blob_sha": "8b945e745cc3611e378cf126d26ecf0a6e9f4e02"},
        ],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
