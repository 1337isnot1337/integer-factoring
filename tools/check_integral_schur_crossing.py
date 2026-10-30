#!/usr/bin/env python3
"""Exact checks of the integral Schur lattice and its crossing families.

All rational power-sum calculations are completed before reduction modulo
two.  In the complete-piece calculation, a common denominator is retained
inside the modulus, so no even denominator is inverted.  These are finite
checks of the accompanying proof, not asymptotic rank certificates.
"""

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
import argparse
import json

from check_partition_channels import (
    column_matrix, dual_polynomial, operator, partition_denominator,
    partitions, rank_mod,
)
from check_modular_crossing import inverse_two, modular_operator


@lru_cache(None)
def complete_symmetric(weight):
    if weight < 0:
        return ()
    return tuple((mu, Fraction(1, partition_denominator(mu)))
                 for mu in partitions(weight))


def multiply_power_sums(left, right):
    result = defaultdict(Fraction)
    for a, c in left:
        for b, d in right:
            result[tuple(sorted(a + b))] += c * d
    return tuple((mu, c) for mu, c in result.items() if c)


@lru_cache(None)
def schur(mu):
    """Jacobi--Trudi, independent of the boson operator recursion."""
    descending = mu[::-1]
    length = len(mu)
    partial = {0: (((), Fraction(1)),)}
    for row in range(length):
        next_partial = {}
        for mask, polynomial in partial.items():
            for column in range(length):
                if mask >> column & 1:
                    continue
                weight = descending[row] - row + column
                if weight < 0:
                    continue
                sign = (-1) ** ((mask >> (column + 1)).bit_count())
                new_mask = mask | (1 << column)
                out = defaultdict(Fraction, next_partial.get(new_mask, ()))
                for part, c in multiply_power_sums(
                        polynomial, complete_symmetric(weight)):
                    out[part] += sign * c
                next_partial[new_mask] = tuple((p, c) for p, c in out.items() if c)
        partial = next_partial
    return partial[(1 << length) - 1]


@lru_cache(None)
def character_table(weight):
    basis = partitions(weight)
    rows = []
    for mu in basis:
        coefficients = dict(schur(mu))
        row = []
        for nu in basis:
            value = coefficients.get(nu, 0) * partition_denominator(nu)
            assert Fraction(value).denominator == 1
            row.append(int(value))
        rows.append(tuple(row))
    for i in range(len(basis)):
        for j in range(len(basis)):
            value = sum(Fraction(rows[i][k] * rows[j][k], partition_denominator(mu))
                        for k, mu in enumerate(basis))
            assert value == int(i == j)
    return tuple(rows)


@lru_cache(None)
def scaled_schurs(weight):
    """Integral power-sum coefficients of weight! times each Schur state."""
    denominator = factorial(weight)
    basis = partitions(weight)
    characters = character_table(weight)
    return tuple(tuple((mu, denominator * characters[i][j] // partition_denominator(mu))
                       for j, mu in enumerate(basis) if characters[i][j])
                 for i in range(len(basis)))


def gram_bits(weight):
    """Homogeneous vacuum pairing in the Schur basis, reduced only at the end."""
    scale = factorial(weight) ** 2
    modulus = 2 * scale
    states = scaled_schurs(weight)
    result = []
    for left in states:
        row = 0
        for j, right in enumerate(states):
            value = sum(c * d * modular_operator(modulus, mu, (), nu)
                        for mu, c in left for nu, d in right) % modulus
            assert value % scale == 0
            row |= (value // scale) << j
        result.append(row)
    return result


def scaled_duals(weight, homogeneous=False):
    states = scaled_schurs(weight)
    if homogeneous:
        inverse = inverse_two(gram_bits(weight), len(states))
        result = []
        for row in inverse:
            out = defaultdict(int)
            for j, state in enumerate(states):
                if row >> j & 1:
                    for mu, c in state:
                        out[mu] += c
            result.append(tuple((mu, c) for mu, c in out.items() if c))
        return tuple(result)
    result = []
    for state in states:
        out = defaultdict(int)
        for mu, c in state:
            for nu, d in dual_polynomial(mu):
                out[nu] += c * d
        result.append(tuple((mu, c) for mu, c in out.items() if c))
    return tuple(result)


def check_duality(weight, homogeneous=False):
    scale = factorial(weight) ** 2
    modulus = 2 * scale
    states = scaled_schurs(weight)
    for i, left in enumerate(scaled_duals(weight, homogeneous)):
        for j, right in enumerate(states):
            value = sum(c * d * modular_operator(modulus, mu, (), nu)
                        for mu, c in left for nu, d in right) % modulus
            assert value % scale == 0 and value // scale == int(i == j)
    return len(states) ** 2


def piece_matrices(weight, homogeneous=False, only_shift=None, exchanged=True):
    """Source Q in either order; rows (x,y), columns (v,u)."""
    basis = partitions(weight)
    dimension = len(basis)
    denominator = factorial(weight)
    scale = denominator ** 3
    modulus = 2 * scale
    states = scaled_schurs(weight)
    duals = scaled_duals(weight, homogeneous)
    characters = character_table(weight)
    shifts = range(-2 * weight, 2 * weight + 1) if only_shift is None else [only_shift]
    pieces = {shift: [[0] * dimension ** 2 for _ in range(dimension ** 2)]
              for shift in shifts}
    for energy in range(4 * weight + 1):
        indices = []
        for index in range(2 * weight + 1):
            shift = energy - index if exchanged else 2 * weight - energy - index
            if shift in pieces and (not homogeneous or comb(2 * weight, index) % 2):
                indices.append((index, shift))
        if not indices:
            continue
        for intermediate in partitions(energy):
            right = []
            for state in states:
                for target in range(dimension):
                    if exchanged:
                        value = sum(c * characters[target][j] * modular_operator(
                            modulus, mu, rho, intermediate)
                            for mu, c in state for j, rho in enumerate(basis))
                    else:
                        value = sum(c * d * modular_operator(modulus, mu, intermediate, nu)
                                    for mu, c in state for nu, d in states[target])
                    right.append(value % modulus)
            nonzero_right = [(j, c) for j, c in enumerate(right) if c]
            if not nonzero_right:
                continue
            left_by_weight = []
            for dual in duals:
                for target in range(dimension):
                    values = [0] * (weight + 1)
                    for mu, c in dual:
                        if exchanged:
                            value = sum(d * modular_operator(modulus, mu, intermediate, nu)
                                        for nu, d in states[target])
                        else:
                            value = sum(characters[target][j] * modular_operator(
                                modulus, mu, rho, intermediate)
                                for j, rho in enumerate(basis))
                        values[sum(mu)] = (values[sum(mu)] + c * value) % modulus
                    left_by_weight.append(values)
            for index, shift in indices:
                matrix = pieces[shift]
                left = [sum(comb(weight + w, index) * c
                            for w, c in enumerate(values) if weight + w >= index) % modulus
                        for values in left_by_weight]
                for left_index, c in enumerate(left):
                    if not c:
                        continue
                    x, first = divmod(left_index, dimension)
                    for right_index, d in nonzero_right:
                        v, second = divmod(right_index, dimension)
                        y, u = (first, second) if exchanged else (second, first)
                        row, column = x * dimension + y, v * dimension + u
                        matrix[row][column] = (matrix[row][column]
                                               + (-1) ** index * c * d) % modulus
    for shift, matrix in pieces.items():
        for row in matrix:
            for j, value in enumerate(row):
                assert value % scale == 0, ("nonintegral complete piece", weight, shift)
                row[j] = value // scale
    if only_shift is None:
        total = [[sum(matrix[r][c] for matrix in pieces.values()) % 2
                  for c in range(dimension ** 2)] for r in range(dimension ** 2)]
        assert total == [[int(r == c) for c in range(dimension ** 2)]
                         for r in range(dimension ** 2)]
    return pieces


def check_integral_factors():
    """Direct rational arithmetic checks individual integral mode values."""
    states = tuple(mu for energy in range(5) for mu in partitions(energy))
    count = 0
    for weight in range(1, 4):
        denominator = factorial(weight)
        for dual in scaled_duals(weight):
            for output in states:
                out_basis = partitions(sum(output))
                output_index = out_basis.index(output)
                output_characters = character_table(sum(output))[output_index]
                for input_ in states:
                    input_terms = schur(input_)
                    by_weight = defaultdict(Fraction)
                    for mu, c in dual:
                        for nu, d in input_terms:
                            for k, rho in enumerate(out_basis):
                                by_weight[sum(mu)] += (
                                    Fraction(c, denominator) * d * output_characters[k]
                                    * operator(mu, rho, nu))
                    for index in range(2 * weight + 1):
                        value = sum(comb(weight + w, index) * c
                                    for w, c in by_weight.items() if weight + w >= index)
                        assert Fraction(value).denominator == 1
                        count += 1
    return count


def profile(pieces, dimension):
    return [[shift, rank_mod(matrix, 2), rank_mod(column_matrix(matrix, dimension), 2)]
            for shift, matrix in sorted(pieces.items())]


def verify(include_eight=False):
    reports = []
    duality_entries = exchange_entries = complete_entries = 0
    factor_entries = check_integral_factors()
    for homogeneous, weights in ((False, range(1, 5)), (True, (2, 4))):
        for weight in weights:
            duality_entries += check_duality(weight, homogeneous)
            pieces = piece_matrices(weight, homogeneous)
            dimension = len(partitions(weight))
            complete_entries += len(pieces) * dimension ** 4
            if weight <= 3 or homogeneous:
                other_order = piece_matrices(weight, homogeneous, exchanged=False)
                assert other_order == pieces
                exchange_entries += len(pieces) * dimension ** 4
            ranks = profile(pieces, dimension)
            reports.append({
                "weight": weight, "homogeneous": homogeneous,
                "dimension": dimension, "profile": ranks,
                "maximum_rank_product": max(r * c for _, r, c in ranks),
            })
            modular_operator.cache_clear()
    if include_eight:
        duality_entries += check_duality(8, True)
        pieces = piece_matrices(8, True, only_shift=0)
        ranks = profile(pieces, 22)
        assert ranks == [[0, 228, 228]]
        reports.append({"weight": 8, "homogeneous": True, "dimension": 22,
                        "central_only": True, "profile": ranks})
        complete_entries += 22 ** 4
        modular_operator.cache_clear()
    assert next(r for r in reports if r["weight"] == 4 and r["homogeneous"])[
        "maximum_rank_product"] == 144
    # Scaling the dual by zero destroys an actual diagonal pairing.
    assert gram_bits(2) != [0, 0]
    return {
        "status": "pass",
        "scope": "finite integral-coordinate and complete-piece checks; no improved exponent claimed",
        "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
        "individual_integral_angle_factors": factor_entries,
        "duality_entries": duality_entries,
        "complete_piece_integrality_entries": complete_entries,
        "both_order_complete_piece_entries": exchange_entries,
        "families": reports,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    parser.add_argument("--homogeneous-eight-central", action="store_true",
                        help="also compute the full 22-state homogeneous central piece at h=8")
    args = parser.parse_args()
    encoded = json.dumps(verify(args.homogeneous_eight_central), indent=2, sort_keys=True) + "\n"
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(encoded)
    print(encoded, end="")
