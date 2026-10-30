#!/usr/bin/env python3
"""Exact quotient-symbol checks for the normalized modular crossing.

The reference route uses the original integral dual polynomial and operator.
The new route uses a Boolean-zeta monomial determinant over F_2.  Finite
checks are implementation evidence, not an asymptotic factoring theorem.
"""

from functools import lru_cache
from itertools import combinations
from math import comb
import argparse
import json

from check_partition_channels import dual_polynomial, operator
from check_modular_crossing import (
    column_bits, determinant_operator_two, rank_two, regular,
    regular_partitions,
)


def multiply_two(left: int, right: int) -> int:
    """Multiply F_2[Z] polynomials encoded by their coefficient bits."""
    result = 0
    while right:
        low = right & -right
        result ^= left << (low.bit_length() - 1)
        right ^= low
    return result


@lru_cache(None)
def binomial_two(exponent: int) -> int:
    """The complete polynomial (1+Z)^exponent, without integer expansion."""
    result = 1
    shift = 1
    while exponent:
        if exponent & 1:
            result ^= result << shift
        exponent >>= 1
        shift <<= 1
    return result


def determinant_monomials_two(exponents: list[list[int]]) -> int:
    """det[Z^e_ij] over F_2[Z], by subset dynamic programming."""
    size = len(exponents)
    values = {0: 1}
    for row in exponents:
        next_values = {}
        for mask, value in values.items():
            for column in range(size):
                if not (mask >> column) & 1:
                    target = mask | (1 << column)
                    next_values[target] = next_values.get(target, 0) ^ (
                        value << row[column])
        values = {mask: value for mask, value in next_values.items() if value}
    return values.get((1 << size) - 1, 0)


@lru_cache(None)
def angle_symbol_two(parts: tuple[int, ...], output: tuple[int, ...],
                     input_: tuple[int, ...]) -> int:
    """All source angle coefficients for one regular label and two states."""
    assert regular(parts, 2) and regular(output, 2) and regular(input_, 2)
    length = len(parts)
    unchanged_count = len(output) + len(input_) - length
    if unchanged_count < 0 or unchanged_count % 2:
        return 0
    unchanged_count //= 2
    common = tuple(sorted(set(output) & set(input_)))
    if unchanged_count > len(common):
        return 0
    determinants = 0
    for unchanged in combinations(common, unchanged_count):
        unchanged = set(unchanged)
        columns = [(True, (k - 1) // 2) for k in output if k not in unchanged]
        columns += [(False, (k - 1) // 2) for k in input_ if k not in unchanged]
        exponents = [
            [2 * (a & (c if creation else ~c)) for creation, c in columns]
            for a in ((m - 1) // 2 for m in parts)
        ]
        determinants ^= determinant_monomials_two(exponents)
    return multiply_two(binomial_two(sum(parts) + length), determinants)


def reference_angle_two(parts: tuple[int, ...], output: tuple[int, ...],
                        input_: tuple[int, ...]) -> int:
    """Independent original integer expansion, including discarded monomials."""
    result = 0
    weight = sum(parts)
    for term, coefficient in dual_polynomial(parts):
        if (coefficient * operator(term, output, input_)) % 2:
            # Deliberately use integer binomial coefficients on this route.
            result ^= sum((comb(weight + sum(term), i) % 2) << i
                          for i in range(weight + sum(term) + 1))
    return result


def source_zeta_two(parts: tuple[int, ...]) -> frozenset[tuple[int, ...]]:
    """Reduce the original, integral source dual in the parameter quotient."""
    return frozenset(term for term, coefficient in dual_polynomial(parts)
                     if coefficient % 2 and regular(term, 2))


def boolean_zeta_two(parts: tuple[int, ...]) -> frozenset[tuple[int, ...]]:
    terms = {()}
    for part in parts:
        a = (part - 1) // 2
        choices = tuple(2 * b + 1 for b in range(a + 1) if b & ~a == 0)
        next_terms = set()
        for term in terms:
            for choice in choices:
                if choice in term:
                    continue
                target = tuple(sorted(term + (choice,)))
                if target in next_terms:
                    next_terms.remove(target)
                else:
                    next_terms.add(target)
        terms = next_terms
    return frozenset(terms)


def normalized_central_two(weight: int, labels=None) -> list[int]:
    """Complete central piece, using the original filtered source duals."""
    if labels is None:
        labels = regular_partitions(2, weight)
    dimension = len(labels)
    result = [0] * (dimension * dimension)
    for energy in range(0, 2 * weight + 1, 2):
        for intermediate in regular_partitions(2, energy):
            right = sum(determinant_operator_two(v, u, intermediate)
                        << (vi * dimension + ui)
                        for vi, v in enumerate(labels)
                        for ui, u in enumerate(labels))
            if not right:
                continue
            for xi, x in enumerate(labels):
                for yi, y in enumerate(labels):
                    if (angle_symbol_two(x, intermediate, y) >> energy) & 1:
                        result[xi * dimension + yi] ^= right
    return result


def verify(max_weight: int = 9) -> dict:
    states = tuple(x for h in range(max_weight + 1)
                   for x in regular_partitions(2, h))
    labels = tuple(x for x in states if x)
    angle_checks = adjoint_checks = zeta_checks = 0
    for parts in labels:
        assert source_zeta_two(parts) == boolean_zeta_two(parts)
        zeta_checks += 1
        h = sum(parts)
        for output in states:
            for input_ in states:
                actual = angle_symbol_two(parts, output, input_)
                expected = reference_angle_two(parts, output, input_)
                assert actual == expected, (parts, output, input_, actual, expected)
                assert actual.bit_length() <= 2 * h + 1
                assert all(not ((actual >> i) & 1) for i in range(1, 2 * h + 1, 2))
                reverse = angle_symbol_two(parts, input_, output)
                reflected = sum(((reverse >> (2 * h - i)) & 1) << i
                                for i in range(2 * h + 1))
                assert actual == reflected
                assert ((actual >> (2 * h)) & 1) == operator(parts, output, input_) % 2
                assert (actual & 1) == operator(parts, input_, output) % 2
                angle_checks += 1
                adjoint_checks += 2

    ideal_checks = 0
    for prime in (2, 3, 5):
        regular_states = tuple(x for h in range(8)
                               for x in regular_partitions(prime, h))
        parameter_terms = [(prime * j,) for j in range(1, 4)]
        parameter_terms += [(m,) * prime for m in range(1, 4)]
        parameter_terms += [tuple(sorted(term + (1,))) for term in parameter_terms]
        for parts in parameter_terms:
            for output in regular_states:
                for input_ in regular_states:
                    assert operator(parts, output, input_) % prime == 0
                    ideal_checks += 1

    # Compare complete piece entries through a second, pre-existing evaluator.
    from check_modular_crossing import normalized_piece
    ranks = []
    for weight in (6, 8):
        basis = regular_partitions(2, weight)
        fast = normalized_central_two(weight, basis)
        original = normalized_piece(2, weight, 0, basis)
        packed = [sum(value << j for j, value in enumerate(row)) for row in original]
        assert fast == packed
        ranks.append({"weight": weight, "labels": len(basis),
                      "row_rank": rank_two(fast),
                      "column_rank": rank_two(column_bits(fast, len(basis)))})
    return {
        "status": "pass",
        "scope": "finite exact quotient, angle-symbol, transpose, and complete-piece checks; no improved exponent",
        "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
        "checks": {"operator_parameter_ideal_entries": ideal_checks,
                   "full_angle_polynomials": angle_checks,
                   "endpoint_transpose_entries": adjoint_checks,
                   "source_boolean_zeta_expansions": zeta_checks},
        "normalized_central_pieces": ranks,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-weight", type=int, default=9)
    parser.add_argument("--central", type=int, help="record one additional complete central piece")
    parser.add_argument("--output")
    args = parser.parse_args()
    report = verify(args.max_weight)
    if args.central:
        basis = regular_partitions(2, args.central)
        matrix = normalized_central_two(args.central, basis)
        report["additional_central_piece"] = {
            "weight": args.central, "labels": len(basis),
            "row_rank": rank_two(matrix),
            "column_rank": rank_two(column_bits(matrix, len(basis))),
            "dual_family": "original filtered source duals reduced in the parameter quotient",
        }
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(encoded)
    print(encoded, end="")
