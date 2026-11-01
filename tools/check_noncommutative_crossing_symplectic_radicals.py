#!/usr/bin/env python3
"""Closed sparse bases, Gram ranks, and half-orientation symplectic radicals."""

from itertools import combinations, product
from math import comb, log2
from pathlib import Path
import argparse
import hashlib
import json

from check_modular_crossing import rank_two
from explore_symplectic_boundary_codes import boundary_code, half_orientation_code


def choose(n, k):
    return comb(n, k) if n >= 0 and 0 <= k <= n else 0


def toggle(result, value):
    if value in result:
        result.remove(value)
    else:
        result.add(value)


def omega(poly, pairs):
    result = set()
    for source in poly:
        for i in range(pairs):
            mask = 3 << (2 * i)
            if not source & mask:
                toggle(result, source | mask)
    return result


def contraction(poly, pairs):
    result = set()
    for source in poly:
        for i in range(pairs):
            mask = 3 << (2 * i)
            if source & mask == mask:
                toggle(result, source ^ mask)
    return result


def mode_complement(mask, pairs):
    result = 0
    for i in range(pairs):
        value = (mask >> (2 * i)) & 3
        result |= (0, 2, 1, 3)[value] << (2 * i)
    return result


def dimensions(pairs, length):
    boundary = gram = radical = 0
    for singles in range(pairs + 1):
        if (length - singles) & 1:
            continue
        even = pairs - singles
        filled = (length - singles) // 2
        if not even:
            continue
        multiplicity = choose(pairs, singles) * (1 << singles)
        b = choose(even - 1, filled - 1)
        g = b if even & 1 else choose(even - 2, filled - 2)
        r = 0 if even & 1 else choose(even - 2, filled - 1)
        assert b == g + r
        boundary += multiplicity * b
        gram += multiplicity * g
        radical += multiplicity * r
    fixed = choose(pairs - 2, length // 2 - 1) if not (pairs | length) & 1 else 0
    assert (radical - fixed) % 2 == 0
    return boundary, gram, radical, fixed, (radical - fixed) // 2


def closed_radical(pairs, length):
    ambient = tuple(sum(1 << i for i in subset)
                    for subset in combinations(range(2 * pairs), length))
    positions = {mask: index for index, mask in enumerate(ambient)}
    vectors, pivots, sectors, half = [], [], [], []
    for singles in range(pairs + 1):
        if (length - singles) & 1:
            continue
        even = pairs - singles
        filled = (length - singles) // 2
        if even < 2 or even & 1 or not 1 <= filled < even:
            continue
        for odd_sites in combinations(range(pairs), singles):
            odd_set = set(odd_sites)
            even_sites = tuple(i for i in range(pairs) if i not in odd_set)
            p, q = even_sites[0], even_sites[-1]
            for orientation in product((1, 2), repeat=singles):
                odd_mask = sum(value << (2 * i) for i, value in zip(odd_sites, orientation))
                for selected in combinations(even_sites[1:-1], filled - 1):
                    seed = odd_mask | (3 << (2 * p))
                    pivot = odd_mask | (3 << (2 * q))
                    for i in selected:
                        seed |= 3 << (2 * i)
                        pivot |= 3 << (2 * i)
                    polynomial = omega(contraction({seed}, pairs), pairs)
                    # A second construction uses the closed Johnson-neighbor formula.
                    occupied = {p, *selected}
                    explicit = {seed} if filled & 1 else set()
                    for removed in occupied:
                        for added in set(even_sites) - occupied:
                            explicit.add(seed ^ (3 << (2 * removed)) ^ (3 << (2 * added)))
                    assert polynomial == explicit
                    assert len(polynomial) == filled * (even - filled) + (filled & 1)
                    assert not omega(polynomial, pairs)
                    assert not contraction(polynomial, pairs)
                    vector = sum(1 << positions[mask] for mask in polynomial)
                    assert vector.bit_length() - 1 == positions[pivot]
                    vectors.append(vector)
                    pivots.append(positions[pivot])
                    sectors.append((odd_sites, orientation))
                    if odd_sites and orientation[0] == 1:
                        half.append(len(vectors) - 1)
    for i, value in enumerate(vectors):
        assert all(((value >> pivot) & 1) == int(i == j) for j, pivot in enumerate(pivots))
    groups = {}
    for i, sector in enumerate(sectors):
        groups.setdefault(sector, []).append(i)
    for indices in groups.values():
        assert all((vectors[i] & vectors[j]).bit_count() % 2 == 0
                   for i, j in product(indices, repeat=2))
    for i in half:
        shifted = sum(1 << positions[mode_complement(ambient[j], pairs)]
                      for j in range(len(ambient)) if vectors[i] >> j & 1)
        assert all(not (shifted >> pivots[j] & 1) for j in half)
    return ambient, vectors, pivots, half


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = []
    for pairs in range(2, 8):
        for length in range(2, pairs + 1):
            boundary, gram, radical, fixed, half_dimension = dimensions(pairs, length)
            ambient, vectors, pivots, half = closed_radical(pairs, length)
            assert len(vectors) == radical and len(half) == half_dimension
            assert rank_two(vectors) == radical
            record = {"pairs": pairs, "length": length, "ambient_dimension": len(ambient),
                      "boundary_dimension": boundary, "boundary_Gram_rank": gram,
                      "radical_dimension": radical, "J_fixed_radical_dimension": fixed,
                      "half_orientation_dimension": half_dimension,
                      "maximum_input_terms": max(map(int.bit_count, vectors), default=0)}
            if pairs <= 6:
                other = boundary_code(pairs, length, radical=True)
                assert other["original_boundary_dimension"] == boundary
                assert other["original_boundary_Gram_rank"] == gram
                assert set(zip(other["output_indices"], other["input_vectors"])) == set(zip(pivots, vectors))
                selected = half_orientation_code(other)
                assert set(zip(selected["output_indices"], selected["input_vectors"])) == {
                    (pivots[i], vectors[i]) for i in half}
                record["independent_RREF_constructor_agrees_exactly"] = True
            cases.append(record)
    asymptotic = []
    for pairs in (16, 32, 64, 128, 256, 512):
        length = round(451 * pairs / 1000)
        if (pairs - length) & 1:
            length += 1
        _, _, radical, fixed, half_dimension = dimensions(pairs, length)
        ambient = comb(2 * pairs, length)
        alpha = length / (2 * pairs)
        limit = alpha * alpha * (1 - alpha) ** 2 / (alpha * alpha + (1 - alpha) ** 2) ** 2
        asymptotic.append({"pairs": pairs, "length": length, "alpha": alpha,
                           "ambient_dimension": str(ambient), "radical_dimension": str(radical),
                           "half_orientation_dimension": str(half_dimension),
                           "radical_fraction": radical / ambient,
                           "asymptotic_radical_fraction_at_same_alpha": limit,
                           "half_fraction": half_dimension / ambient,
                           "label_bits_per_physical_mode": log2(half_dimension) / (2 * pairs)})
    result = {
        "status": "All closed bases, coordinate projections, Gram radicals, J selectors, and independent RREF comparisons passed.",
        "scope": "All-size ordinary proof with exact finite audits; no improved rank-sum exponent claimed.",
        "finite_cases": cases,
        "entropy_controls": asymptotic,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()

