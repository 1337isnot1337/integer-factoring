"""Exact focused checks for the complete-binary-core Schur identity.

This checks an identity and its separated-support rank consequence.  It
makes no asymptotic rank or factoring-runtime claim.
"""
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

from modular_pfaffian import final_piece_profile, polynomial_product
from check_modular_pfaffian import binary_rank


def polynomial_power(value, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = polynomial_product(result, value)
        value = polynomial_product(value, value)
        exponent >>= 1
    return result


def lift(label, core):
    return tuple(sorted(tuple(range(1, 2 * core, 2)) +
                        tuple(2 * core + t for t in label), reverse=True))


def nonzero_bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def matrices_from_profiles(profiles, dimension):
    matrices = {}
    for (x, y, u, v), polynomial in profiles.items():
        for exponent in nonzero_bits(polynomial):
            matrices.setdefault(exponent, [0] * (dimension * dimension))[
                x * dimension + y] |= 1 << (u * dimension + v)
    return matrices


def rank_pairs(matrices, d):
    pairs = set()
    for matrix in matrices.values():
        column = [sum(((matrix[x*d+y] >> (u*d+v)) & 1) << (y*d+v)
                      for y in range(d) for v in range(d))
                  for x in range(d) for u in range(d)]
        pairs.add((binary_rank(matrix), binary_rank(column)))
    return sorted(pairs)


def one_case(name, core, labels, expect_overlap, expect_separation):
    assert core >= 1 and core & (core - 1) == 0
    h0 = sum(labels[0])
    r = len(labels[0])
    assert all(sum(label) == h0 and len(label) == r for label in labels)
    assert all(0 < t < 2 * core and t % 2 for label in labels for t in label)
    transformed = tuple(lift(label, core) for label in labels)
    h = core * core + 2 * core * r + h0
    assert all(sum(label) == h for label in transformed)
    d = len(labels)
    kernel = 1 | (1 << (2 * core)) | (1 << (4 * core))
    multiplier = polynomial_power(kernel, core + 2 * r)
    assert multiplier & 1
    assert multiplier.bit_count() % 2 == 1
    low_profiles = {}
    high_profiles = {}
    comparisons = 0
    for x, y, u, v in product(range(d), repeat=4):
        low = final_piece_profile(labels[x], labels[v], labels[u], labels[y])
        high = final_piece_profile(transformed[x], transformed[v],
                                   transformed[u], transformed[y])
        assert high == polynomial_product(multiplier, low)
        assert high.bit_count() % 2 == int(x == v and y == u)
        low_profiles[x, y, u, v] = low
        high_profiles[x, y, u, v] = high
        comparisons += 1
    low_matrices = matrices_from_profiles(low_profiles, d)
    high_matrices = matrices_from_profiles(high_profiles, d)
    support_counts = {}
    for shift in nonzero_bits(multiplier):
        for exponent in low_matrices:
            position = shift + exponent
            support_counts[position] = support_counts.get(position, 0) + 1
    overlap_positions = sum(count > 1 for count in support_counts.values())
    assert bool(overlap_positions) == expect_overlap
    assert (core > 2 * h0) == expect_separation
    low_pairs = rank_pairs(low_matrices, d)
    high_pairs = rank_pairs(high_matrices, d)
    if expect_separation:
        assert not overlap_positions
        assert low_pairs == high_pairs
        originals = {tuple(matrix) for matrix in low_matrices.values()}
        assert all(tuple(matrix) in originals for matrix in high_matrices.values())
        assert {tuple(matrix) for matrix in high_matrices.values()} == originals
    result = {
        "name": name, "core_size": core, "active_cardinality": r,
        "low_weight": h0, "lifted_weight": h, "labels": [list(t) for t in labels],
        "dimension": d, "profile_comparisons": comparisons,
        "overlapping_coefficient_positions": overlap_positions,
        "separation_condition": expect_separation,
        "low_rank_pairs": low_pairs, "lifted_rank_pairs": high_pairs,
        "low_nonzero_pieces": len(low_matrices),
        "lifted_nonzero_pieces": len(high_matrices),
    }
    polynomial_product.cache_clear()
    return result


def main():
    cases = [
        one_case("small_core", 2, ((1,),), False, False),
        one_case("overlapping_supports", 8,
                 ((7, 5, 3), (9, 5, 1), (11, 3, 1)), True, False),
        one_case("separated_supports", 32, ((5, 3), (7, 1)), False, True),
    ]
    print(json.dumps({
        "cases": cases,
        "profile_comparisons": sum(c["profile_comparisons"] for c in cases),
        "implementation_sha256": sha256(Path(__file__).with_name(
            "modular_pfaffian.py").read_bytes()).hexdigest(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "status": "all exact checks passed",
        "scope": "profile convolution and separated-support rank equality; no new factoring exponent",
    }, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
