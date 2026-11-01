#!/usr/bin/env python3
"""Exact integral adjoints and F_2 Schur fusion/central-Gram checks.

No floating point or inversion of an even denominator is used.  The
all-weight theorem and its scope are in research/schur-fusion-gram.md;
the finite ranks are not an asymptotic rank certificate.
"""

from collections import defaultdict
from functools import lru_cache
from math import factorial
import argparse
import json
from pathlib import Path

from check_partition_channels import (
    dual_polynomial, operator, partition_denominator, partitions,
)
from check_homogeneous_duals import shift
from check_modular_crossing import inverse_two, modular_operator, rank_two
from check_integral_schur_crossing import (
    character_table, gram_bits, piece_matrices,
)


@lru_cache(None)
def rim_hooks(partition, size):
    """Murnaghan--Nakayama removals via beta-set bead moves."""
    descending = partition[::-1]
    length = len(descending)
    beta = tuple(value + length - 1 - i
                 for i, value in enumerate(descending))
    result = []
    for i, bead in enumerate(beta):
        moved = bead - size
        if moved < 0 or moved in beta:
            continue
        sign = (-1) ** sum(moved < other < bead for other in beta)
        changed = sorted(beta[:i] + beta[i + 1:] + (moved,), reverse=True)
        remainder = tuple(sorted(value - length + 1 + j
                                 for j, value in enumerate(changed)
                                 if value - length + 1 + j > 0))
        assert sum(remainder) == sum(partition) - size
        result.append((remainder, sign))
    return tuple(result)


@lru_cache(None)
def character(partition, cycles):
    if not cycles:
        return int(not partition)
    return sum(sign * character(remainder, cycles[:-1])
               for remainder, sign in rim_hooks(partition, cycles[-1]))


@lru_cache(None)
def characters(weight):
    basis = partitions(weight)
    return tuple(tuple(character(mu, nu) for nu in basis) for mu in basis)


def check_characters(weight):
    table = characters(weight)
    columns = tuple(zip(*table))
    cases = 0
    for i, mu in enumerate(partitions(weight)):
        for j in range(i + 1):
            value = sum(a * b for a, b in zip(columns[i], columns[j]))
            assert value == (partition_denominator(mu) if i == j else 0)
            cases += 1
    if weight <= 8:
        assert table == character_table(weight)
    return cases


def matrix_product_mod(left, right, modulus):
    columns = tuple(zip(*right))
    return [[sum(a * b for a, b in zip(row, col)) % modulus
             for col in columns] for row in left]


def transpose_bits(rows, columns):
    return [sum(((row >> j) & 1) << i for i, row in enumerate(rows))
            for j in range(columns)]


def multiply_bits(left, right):
    result = []
    for row in left:
        value = 0
        while row:
            bit = row & -row
            value ^= right[bit.bit_length() - 1]
            row ^= bit
        result.append(value)
    return result


def tensor_first(matrix, dimension):
    return [sum(((matrix[x] >> i) & 1) << (i * dimension + y)
                for i in range(dimension))
            for x in range(dimension) for y in range(dimension)]


def divided_shift(polynomial, direction, degree):
    result = dict(polynomial)
    for k in range(1, degree + 1):
        shifted = shift(result, direction)
        assert all(value % k == 0 for value in shifted.values())
        result = {mu: value // k for mu, value in shifted.items()}
    return result


@lru_cache(None)
def homogeneous_pairing_polynomial(mu):
    """B_h p_mu = G_h omega p_mu, in integer power-sum coordinates."""
    weight = sum(mu)
    result = defaultdict(int)
    for k in range(weight):
        lowered = divided_shift({mu: 1}, -1, k)
        raised = divided_shift(lowered, 1, k)
        for term, value in raised.items():
            result[term] += (-1) ** (weight - len(mu)) * value
    return tuple((mu, value) for mu, value in result.items() if value)


def check_integral_identities():
    states = tuple(mu for weight in range(6) for mu in partitions(weight))
    parameters = tuple(mu for weight in range(5) for mu in partitions(weight))
    adjoints = involutions = pairings = homogenizations = 0
    for mu in parameters:
        twice = defaultdict(int)
        for nu, c in dual_polynomial(mu):
            for rho, d in dual_polynomial(nu):
                twice[rho] += c * d
        assert {key: value for key, value in twice.items() if value} == {mu: 1}
        involutions += 1
        for output in states:
            for input_ in states:
                left = operator(mu, output, input_) * partition_denominator(output)
                right = sum(c * operator(nu, input_, output)
                            for nu, c in dual_polynomial(mu))
                right *= partition_denominator(input_)
                assert left == right, ("adjoint", mu, output, input_)
                adjoints += 1
    for h in range(1, 6):
        for mu in partitions(h):
            homogeneous = dict(homogeneous_pairing_polynomial(mu))
            for nu in partitions(h):
                assert operator(mu, (), nu) == (
                    homogeneous.get(nu, 0) * partition_denominator(nu))
                pairings += 1
            for energy in range(4):
                for input_ in partitions(energy):
                    for output in partitions(energy + h):
                        translated = sum(c * operator(nu, output, input_)
                                         for nu, c in dual_polynomial(mu))
                        homogenized = sum(c * operator(nu, output, input_)
                                          for nu, c in homogeneous.items())
                        assert translated == homogenized, ("homogenization", mu,
                                                           output, input_)
                        homogenizations += 1
    # Omitting J fails even over Q; the source's annihilation signs matter.
    assert operator((2,), (2,), ()) * partition_denominator((2,)) == 2
    assert operator((2,), (), (2,)) == -6
    return {"integral_adjoint_entries": adjoints,
            "integral_involutions": involutions,
            "homogeneous_pairing_entries": pairings,
            "integral_homogenization_entries": homogenizations}


def fusion_matrix(weight):
    basis, output = partitions(weight), partitions(2 * weight)
    dimension = len(basis)
    denominator = factorial(weight)
    scale, modulus = denominator ** 2, 2 * denominator ** 2
    scaled = []
    for row in characters(weight):
        assert all(denominator * value % partition_denominator(mu) == 0
                   for mu, value in zip(basis, row))
        scaled.append([denominator * value // partition_denominator(mu)
                       for mu, value in zip(basis, row)])
    transposed = tuple(zip(*scaled))
    raw = []
    for rho in output:
        source = [[modular_operator(modulus, a, rho, b) for b in basis]
                  for a in basis]
        value = matrix_product_mod(scaled, source, modulus)
        value = matrix_product_mod(value, transposed, modulus)
        raw.append([entry for row in value for entry in row])
    converted = matrix_product_mod(characters(2 * weight), raw, modulus)
    assert all(value % scale == 0 for row in converted for value in row)
    return [sum((value // scale) << j for j, value in enumerate(row))
            for row in converted]


def check_fusion(weight):
    dimension = len(partitions(weight))
    ambient = len(partitions(2 * weight))
    matrix = fusion_matrix(weight)
    fusion_rank = rank_two(matrix)
    gram = multiply_bits(transpose_bits(matrix, dimension ** 2), matrix)
    pairing = gram_bits(weight)
    inverse = inverse_two(pairing, dimension)
    central = multiply_bits(tensor_first(inverse, dimension), gram)
    central = multiply_bits(central, transpose_bits(
        tensor_first(pairing, dimension), dimension ** 2))
    vacuum = sum(1 << (x * dimension + x) for x in range(dimension))
    central = [row ^ (vacuum if i // dimension == i % dimension else 0)
               for i, row in enumerate(central)]
    assert transpose_bits(pairing, dimension) == pairing
    assert transpose_bits(inverse, dimension) == inverse
    vec_pairing = sum(row << (x * dimension) for x, row in enumerate(pairing))
    vec_inverse = sum(row << (x * dimension) for x, row in enumerate(inverse))
    column = [row ^ (vec_pairing if (vec_inverse >> i) & 1 else 0)
              for i, row in enumerate(gram)]
    reshuffled = [sum(((central[x * dimension + y] >> (v * dimension + u)) & 1)
                      << (v * dimension + y)
                      for v in range(dimension) for y in range(dimension))
                  for x in range(dimension) for u in range(dimension)]
    assert reshuffled == column
    transformed = multiply_bits(tensor_first(pairing, dimension), central)
    transformed = multiply_bits(transformed, transpose_bits(
        tensor_first(inverse, dimension), dimension ** 2))
    assert transformed == transpose_bits(column, dimension ** 2)
    direct_entries = 0
    if weight <= 4:
        direct = piece_matrices(weight, homogeneous=True, only_shift=0)[0]
        assert central == [sum(value << j for j, value in enumerate(row))
                           for row in direct]
        direct_entries = dimension ** 4
    central_rank = rank_two(central)
    assert central_rank == rank_two(column)
    gram_rank = rank_two(gram)
    lower = max(0, 2 * fusion_rank - ambient - 1)
    upper = min(dimension ** 2, fusion_rank + 1)
    assert lower <= central_rank <= upper
    assert max(0, 2 * fusion_rank - ambient) <= gram_rank <= fusion_rank
    expected = {1: (1, 0, 1), 2: (3, 1, 2), 4: (16, 12, 12),
                8: (229, 227, 228)}
    assert (fusion_rank, gram_rank, central_rank) == expected[weight]
    if weight == 2:
        # Output order: 1111,211,31,22,4.  The product s_11^2 has bits 11.
        assert transpose_bits(matrix, dimension ** 2) == [5, 8, 8, 18]
        annihilator, ordinary_product = 5, 11
        assert all((annihilator & column).bit_count() % 2 == 0
                   for column in transpose_bits(matrix, dimension ** 2))
        assert (annihilator & ordinary_product).bit_count() % 2 == 1
        assert rank_two(transpose_bits(matrix, dimension ** 2)
                        + [ordinary_product]) == fusion_rank + 1
    return {"h": weight, "dimension": dimension, "ambient_dimension": ambient,
            "fusion_rank": fusion_rank, "fusion_gram_rank": gram_rank,
            "central_rank": central_rank, "proved_lower_bound": lower,
            "proved_upper_bound": upper,
            "integral_fusion_entries": ambient * dimension ** 2,
            "row_column_identity_entries": dimension ** 4,
            "independent_direct_central_entries": direct_entries}


def verify(include_eight=False):
    result = check_integral_identities()
    weights = (1, 2, 4, 8) if include_eight else (1, 2, 4)
    character_weights = sorted(set(weights) | {2 * weight for weight in weights})
    result["integer_character_orthogonality_entries"] = sum(
        check_characters(weight) for weight in character_weights)
    records = []
    for weight in weights:
        records.append(check_fusion(weight))
        modular_operator.cache_clear()
    result["weights"] = records
    result["negative_controls"] = ["omitted adjoint involution",
                                   "ordinary multiplication image containment"]
    result["scope"] = ("exact finite checks of the integral adjoint and central "
                       "fusion-Gram theorem; asymptotic fusion rank remains open")
    result["new_factoring_exponent_claimed"] = False
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight-eight", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    evidence = json.dumps(verify(args.weight_eight), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(evidence)
    print(evidence, end="")
