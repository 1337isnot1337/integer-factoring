#!/usr/bin/env python3
"""Exact checks for homogeneous permutation labels; no asymptotic rank claim."""

import argparse
from functools import lru_cache
from itertools import permutations, product
import json
from math import factorial, isqrt, prod
from pathlib import Path
import random

from fermionic_piece_polynomial import (
    monomial_piece_polynomial,
    polynomial_product,
)


def rank_two(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = row
                break
            row ^= basis[pivot]
    return len(basis)


def inverse(permutation):
    result = [0] * len(permutation)
    for i, value in enumerate(permutation):
        result[value] = i
    return tuple(result)


def compose(first, second):
    return tuple(first[i] for i in second)


def gram_kernel(permutation):
    size = len(permutation)
    used = set()
    nontrivial = []
    for start in range(size):
        if start in used:
            continue
        length, current = 0, start
        while current not in used:
            used.add(current)
            length += 1
            current = permutation[current]
        if length > 1:
            nontrivial.append(length)
    if not nontrivial:
        return int(size % 2 == 0)
    return int(len(nontrivial) == 1 and nontrivial[0] % 2 == 1)


def direct_gram_identity(permutation):
    size = len(permutation)
    rows = [
        sum(int(i != j and i != permutation[j]) << j for j in range(size))
        for i in range(size)
    ]
    return int(rank_two(rows) == size)


@lru_cache(None)
def partitions(size, maximum=None):
    if maximum is None:
        maximum = size
    if not size:
        return ((),)
    return tuple(
        (first,) + tail
        for first in range(min(size, maximum), 0, -1)
        for tail in partitions(size - first, first)
    )


def hooks(shape):
    return tuple(
        row - j + sum(other > j for other in shape[i + 1:])
        for i, row in enumerate(shape)
        for j in range(row)
    )


@lru_cache(None)
def dimension(shape):
    return factorial(sum(shape)) // prod(hooks(shape))


@lru_cache(None)
def single_cycle_character(shape, length):
    """Murnaghan--Nakayama for cycle type (length, 1, ..., 1)."""
    answer = 0
    for smaller in partitions(sum(shape) - length):
        if len(smaller) > len(shape):
            continue
        if any(value > shape[i] for i, value in enumerate(smaller)):
            continue
        cells = {
            (i, j)
            for i, row in enumerate(shape)
            for j in range(smaller[i] if i < len(smaller) else 0, row)
        }
        if not cells:
            continue
        if any(
            (i + 1, j) in cells and (i, j + 1) in cells
            and (i + 1, j + 1) in cells
            for i, j in cells
        ):
            continue
        reached = {next(iter(cells))}
        pending = list(reached)
        while pending:
            i, j = pending.pop()
            for cell in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
                if cell in cells and cell not in reached:
                    reached.add(cell)
                    pending.append(cell)
        if len(reached) == len(cells):
            sign = (-1) ** (len({i for i, _ in cells}) - 1)
            answer += sign * dimension(smaller)
    return answer


def central_eigenvalue(shape):
    size = sum(shape)
    numerator = sum(
        factorial(size) // (length * factorial(size - length))
        * single_cycle_character(shape, length)
        for length in range(3, size + 1, 2)
    )
    value, remainder = divmod(numerator, dimension(shape))
    assert remainder == 0
    return value + int(size % 2 == 0)


def balanced_shape(size):
    side = isqrt(size)
    side -= side % 2
    assert size >= 4 and size % 2 == 0 and side > 0
    return (side + size - side * side,) + (side,) * (side - 1)


def edge_masks(permutation):
    size = len(permutation)
    return tuple((1 << i) | (1 << (size + value))
                 for i, value in enumerate(permutation))


def xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def lifted_polynomial(operator_a, operator_b, output, input_, size):
    """Homogeneous six-block lift, packed with carry-free degree radices."""
    radix = 4 * size + 1
    places = [radix * 5**i for i in range(2 * size)]

    def monomial(mask):
        return 1 << sum(places[i] for i in range(2 * size) if mask >> i & 1)

    def binomial(mask):
        answer = 1
        for i in range(2 * size):
            if mask >> i & 1:
                answer ^= answer << places[i]
        return answer

    vertices = (
        [("U", e) for e in output] + [("B", e) for e in operator_b]
        + [("A", e) for e in operator_a] + [("Y", e) for e in input_]
    )
    count = len(vertices)
    edges = [[0] * count for _ in range(count)]
    for i, (left, e) in enumerate(vertices):
        for j in range(i + 1, count):
            right, f = vertices[j]
            groups = left, right
            value = 0
            if groups == ("U", "B") and e == f:
                value = polynomial_product(3, binomial(f))
            elif groups == ("U", "A") and e == f:
                value = polynomial_product(6, polynomial_product(monomial(e), binomial(f)))
            elif groups == ("U", "Y") and e == f:
                value = monomial(f) << 1
            elif groups == ("B", "A") and not e & f:
                value = monomial(f) << 1
            elif groups == ("B", "Y") and not e & f:
                value = polynomial_product(6, polynomial_product(monomial(f), binomial(e)))
            elif groups == ("A", "Y") and not e & f:
                value = polynomial_product(3, binomial(e))
            edges[i][j] = value

    @lru_cache(None)
    def matching(remaining):
        if not remaining:
            return 1
        first_bit = remaining & -remaining
        first = first_bit.bit_length() - 1
        rest = remaining ^ first_bit
        answer = 0
        candidates = rest
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            other = bit.bit_length() - 1
            if edges[first][other]:
                answer ^= polynomial_product(edges[first][other], matching(rest ^ bit))
        return answer

    return matching((1 << count) - 1)


def specialize(polynomial, size):
    answer = 0
    monomials = 0
    radix = 4 * size + 1
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        code = bit.bit_length() - 1
        degree_l, code = code % radix, code // radix
        assert degree_l % 2 == 0
        exponent = degree_l
        for i in range(2 * size):
            degree = code % 5
            code //= 5
            exponent += degree << (i + 1)
        assert code == 0
        answer ^= 1 << exponent
        monomials += 1
    return answer, monomials


def all_group_matchings(size):
    groups = tuple("U" * size + "B" * size + "A" * size + "Y" * size)
    total = 0

    def walk(remaining, counts):
        nonlocal total
        if not remaining:
            assert counts["UB"] == counts["AY"]
            assert counts["UA"] == counts["BY"]
            assert counts["UY"] == counts["BA"]
            b, c, d = counts["UB"], counts["UA"], counts["UY"]
            assert b + c + d == size
            from math import comb
            for k in range(2 * size + 1):
                shift = k - (size - b)
                value = comb(size - d, shift) % 2 if 0 <= shift <= size - d else 0
                if value:
                    assert size - b <= k <= size + c
                if k == size:
                    assert value == int(not b & c)
            total += 1
            return
        first, *rest = remaining
        for index, other in enumerate(rest):
            if groups[first] == groups[other]:
                continue
            pair = groups[first] + groups[other]
            counts[pair] += 1
            walk(rest[:index] + rest[index + 1:], counts)
            counts[pair] -= 1

    walk(list(range(4 * size)), {key: 0 for key in ("UB", "UA", "UY", "BA", "BY", "AY")})
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    rng = random.Random(19374021)
    counts = {}

    determinant_checks = 0
    for size in range(2, 9):
        for permutation in permutations(range(size)):
            assert direct_gram_identity(permutation) == gram_kernel(permutation)
            determinant_checks += 1
    counts["direct_disjointness_gram_checks"] = determinant_checks

    ranks = []
    for size in range(2, 7):
        labels = tuple(permutations(range(size)))
        indices = {p: i for i, p in enumerate(labels)}
        support = tuple(p for p in labels if gram_kernel(p))
        rows = [sum(1 << indices[compose(label, p)] for p in support) for label in labels]
        rank = rank_two(rows)
        odd_multiplicity = sum(
            dimension(shape) ** 2 for shape in partitions(size)
            if central_eigenvalue(shape) % 2
        )
        assert rank >= odd_multiplicity
        ranks.append({"a": size, "dimension": len(labels), "rank": rank,
                      "odd_central_multiplicity": odd_multiplicity})
    counts["small_gram_ranks"] = ranks

    balanced_characters = 0
    all_characters = 0
    for size in range(2, 15):
        assert sum(dimension(shape) ** 2 for shape in partitions(size)) == factorial(size)
        for shape in partitions(size):
            value = central_eigenvalue(shape)
            all_characters += 1
            even_contents = sum((j - i) % 2 == 0
                                for i, row in enumerate(shape) for j in range(row))
            if size % 2 == 0 and even_contents * 2 == size:
                assert value % 2 == 1
                balanced_characters += 1
    counts["integer_central_characters"] = all_characters
    counts["balanced_content_parity_checks"] = balanced_characters

    hook_shapes = 0
    for size in range(4, 1001, 2):
        shape = balanced_shape(size)
        assert sum(shape) == size and all(row % 2 == 0 for row in shape)
        largest = shape[0] + len(shape) - 1
        assert largest * largest <= 49 * size
        if size <= 28:
            assert central_eigenvalue(shape) % 2 == 1
            hook_product = prod(hooks(shape))
            assert hook_product * hook_product <= (49 * size) ** size
        hook_shapes += 1
    counts["balanced_hook_shapes"] = hook_shapes

    support_eight = {
        p for p in permutations(range(8)) if gram_kernel(p)
    }
    seven_cycle = tuple(list(range(1, 7)) + [0, 7])
    square_coefficient = sum(
        compose(inverse(p), seven_cycle) in support_eight
        for p in support_eight
    ) % 2
    assert square_coefficient == 1
    bad_eigenvalue = central_eigenvalue((4, 3, 2, 1))
    assert bad_eigenvalue == -176
    counts["negative_controls"] = {
        "a8_square_coefficient_at_7_cycle": square_coefficient,
        "a10_staircase_central_eigenvalue": bad_eigenvalue,
    }

    lift_checks = 0
    lifted_monomials = 0
    for size, sample_count in ((2, None), (3, 32), (4, 4)):
        labels = tuple(permutations(range(size)))
        if sample_count is None:
            cases = product(labels, repeat=4)
        else:
            cases = [tuple(rng.choice(labels) for _ in range(4))
                     for _ in range(sample_count)]
        for a, b, u, y in cases:
            masks = tuple(edge_masks(p) for p in (a, b, u, y))
            polynomial = lifted_polynomial(*masks, size)
            specialized, terms = specialize(polynomial, size)
            modes = tuple(tuple(sorted(1 + 2 * e for e in label)) for label in masks)
            direct, shift = monomial_piece_polynomial(*modes)
            h = size + 2 * (2 ** (2 * size) - 1)
            assert shift == -2 * h and specialized == direct
            assert polynomial.bit_count() % 2 == (
                int(u == y) * gram_kernel(compose(inverse(a), b))
            )
            lift_checks += 1
            lifted_monomials += terms
            polynomial_product.cache_clear()
    counts["independent_homogeneous_lift_checks"] = lift_checks
    counts["multigraded_monomials_checked"] = lifted_monomials
    counts["group_matching_count_constraints"] = sum(all_group_matchings(a) for a in range(1, 4))
    counts["scope"] = "large effective homogeneous label space; no improved joint-width bound"

    encoded = json.dumps(counts, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
