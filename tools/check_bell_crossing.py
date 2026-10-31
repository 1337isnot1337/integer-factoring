#!/usr/bin/env python3
"""Replay physical Bell cuts and exact ternary pencil ranks; no Lean work."""

import argparse
from itertools import product
import json
from math import log2
import random

from bell_crossing import (
    bell_basis, common_determinant_floor, direct_pencil_profiles, direct_transform,
    entropy_exceeds, matmul, pencil_profiles, physical_matrix, rank_mod, transform,
)


def check() -> dict:
    rng = random.Random(20261031)
    local_checks = 0
    physical_checks = 0
    rank_checks = 0
    for bits in (1, 2, 3):
        q = 1 << (2 * bits)
        for prime in (3, 5, 7):
            one = (1,) * q
            assert transform(one, bits, prime) == one
            vectors = [one] + [tuple(int(j == i) for j in range(q))
                               for i in range(q)]
            vectors += [tuple(rng.randrange(prime) for _ in range(q))
                        for _ in range(4)]
            basis = bell_basis(bits, prime)
            for values in vectors:
                other = transform(values, bits, prime)
                assert transform(other, bits, prime) == tuple(x % prime for x in values)
                local_checks += 1
                for column, diagonal in ((False, values), (True, other)):
                    matrix = physical_matrix(values, bits, prime, column)
                    expected = [[x * diagonal[j] % prime for j, x in enumerate(row)]
                                for row in basis]
                    assert matmul(matrix, basis, prime) == expected
                    physical_checks += 1
                    if bits <= 2 or values in vectors[-2:]:
                        assert rank_mod(matrix, prime) == sum(x % prime != 0 for x in diagonal)
                        rank_checks += 1
            for column in (False, True):
                assert physical_matrix(one, bits, prime, column) == [
                    [int(i == j) for j in range(q)] for i in range(q)]

    literal_profiles = 0
    for values in product(range(3), repeat=4):
        for length in range(1, 5):
            assert pencil_profiles(values, 1, length) == direct_pencil_profiles(values, 1, length)
            literal_profiles += 1
    direct_transform_checks = 0
    for bits in range(1, 6):
        for prime in (3, 5, 7):
            values = tuple(rng.randrange(prime) for _ in range(1 << (2 * bits)))
            assert transform(values, bits, prime) == direct_transform(values, bits, prime)
            direct_transform_checks += 1
    # Negative control: reducing integer word multiplicities mod p is wrong.
    assert pencil_profiles((1, 1, 1, 1), 1, 3)[0] == (0, 0, 0, 64)
    assert pencil_profiles((0, 0, 0, 0), 1, 3)[0] == (64, 0, 0, 0)
    numerator, denominator = 169, 239
    assert 2 * numerator ** 2 > denominator ** 2
    floors = []
    for q in (4, 16, 64):
        exceeds = entropy_exceeds(q, numerator, denominator)
        assert exceeds == (q in (4, 16))
        gamma = (log2(2 * q) - (q + 1) / (2 * q) * log2(q + 1)) / log2(q)
        floors.append({"symbols": q, "leading_rho_floor_approx": gamma,
                       "exceeds_169_over_239_exactly": exceeds,
                       "length_2N_floor": common_determinant_floor(q, 1)})
    return {"status": "passed", "local_transform_checks": local_checks,
            "physical_eigenbasis_checks": physical_checks,
            "independent_elimination_ranks": rank_checks,
            "literal_tensor_word_profile_checks": literal_profiles,
            "fast_against_direct_transform_checks": direct_transform_checks,
            "small_seed_determinant_floors": floors,
            "factoring_exponent": "1/2; no improvement established"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = check()
    if args.output:
        with open(args.output, "w") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
