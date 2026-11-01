#!/usr/bin/env python3
"""Independent finite checks of mixed interpolation/addition/carry blocks."""

from itertools import product
import json
from random import Random

from check_logarithmic_wire_network import (
    CharacterAlgebra, carry_prefix_query, interpolation_matrix, linear_solve,
)


def segments(indices):
    result = []
    for index in indices:
        if result and index == result[-1][-1] + 1:
            result[-1].append(index)
        else:
            result.append([index])
    return result


def coefficient_map(inputs, size, modulus, inverse):
    coefficients = [
        [sum(inputs[power * size + point] * inverse[degree][point]
             for point in range(size)) % modulus for degree in range(size)]
        for power in range(2 * size)
    ]
    result = []
    for index in range(size ** 2):
        quotient, residue = divmod(index, size)
        value = coefficients[residue][quotient]
        if quotient:
            value += coefficients[size + residue][quotient - 1]
        result.append(value % modulus)
    return tuple(result)


def carry_weight(coefficients, selected, digits, endpoints, base, algebra):
    """Sum original carry settings for fixed coefficients, without Fourier DP."""
    result = algebra.phase(0)
    for component, (incoming_frequency, outgoing_frequency) in zip(segments(selected), endpoints):
        total = algebra.zero
        initial = component[0] == 0
        final = component[-1] == len(digits) - 1
        for incoming in range(1 if initial else algebra.modulus):
            carry = incoming
            for index in component:
                value = coefficients[index] + carry
                if value % base != digits[index]:
                    break
                carry = value // base
                assert 0 <= carry < algebra.modulus
            else:
                if final and carry:
                    continue
                total = algebra.add(total, algebra.phase(
                    incoming_frequency * incoming + outgoing_frequency * carry))
        result = algebra.multiply(result, total)
    return result


def compiled_query(matrix, beta, selected, digits, endpoints, base, algebra):
    source_count = len(beta)
    rank, _ = linear_solve(matrix, [0] * len(matrix), algebra.modulus)
    assert rank == len(selected)
    transpose = [tuple(row) for row in zip(*matrix)]
    _, gamma = linear_solve(transpose, beta, algebra.modulus)
    if gamma is None:
        return algebra.zero
    frequency = dict(zip(selected, gamma))
    result = algebra.phase(0)
    for component, (incoming_frequency, outgoing_frequency) in zip(segments(selected), endpoints):
        query = carry_prefix_query(
            [digits[index] for index in component],
            [frequency[index] for index in component],
            incoming_frequency, outgoing_frequency, base,
            component[0] == 0, component[-1] == len(digits) - 1, algebra,
        )
        result = algebra.multiply(result, query)
    return algebra.scale(result, pow(
        algebra.modulus, source_count - len(selected), algebra.prime))


def verify():
    rng = Random(20261030)
    size, modulus = 2, 3
    count = 2 * size ** 2
    inverse = interpolation_matrix(size, modulus)
    columns = []
    for source in range(count):
        basis = [int(index == source) for index in range(count)]
        columns.append(coefficient_map(basis, size, modulus, inverse))
    full_matrix = [tuple(row) for row in zip(*columns)]
    # This finite cache is only the exhaustive reference input list, not the
    # production block evaluator or its storage strategy.
    reference = [(word, coefficient_map(word, size, modulus, inverse))
                 for word in product(range(modulus), repeat=count)]

    surjective_cases = query_cases = outside_cases = disconnected_cases = 0
    nonzero_cases = direct_settings = 0
    fiber_scalars = set()
    for coefficient_prime in (2, 5, 7):
        algebra = CharacterAlgebra(modulus, coefficient_prime)
        for mask in range(1, 1 << (size ** 2)):
            selected = tuple(index for index in range(size ** 2) if mask >> index & 1)
            matrix = [full_matrix[index] for index in selected]
            rank, _ = linear_solve(matrix, [0] * len(selected), modulus)
            assert rank == len(selected)
            surjective_cases += 1
            transpose = [tuple(row) for row in zip(*matrix)]
            outside = next(tuple(int(index == source) for index in range(count))
                           for source in range(count)
                           if linear_solve(transpose, [int(index == source) for index in range(count)],
                                           modulus)[1] is None)
            for trial in range(3):
                base = 2
                digits = [0] * (size ** 2) if trial == 0 else [
                    rng.randrange(base) for _ in range(size ** 2)]
                endpoints = [(0, 0) if trial == 0 else
                             (rng.randrange(modulus), rng.randrange(modulus))
                             for _ in segments(selected)]
                gamma = [0] * len(selected) if trial == 0 else [
                    rng.randrange(modulus) for _ in selected]
                beta = tuple(sum(gamma[row] * matrix[row][column]
                                 for row in range(len(selected))) % modulus
                             for column in range(count))
                if trial == 2:
                    beta = tuple((a + b) % modulus for a, b in zip(beta, outside))
                    assert linear_solve(transpose, beta, modulus)[1] is None
                    outside_cases += 1
                expected = compiled_query(
                    matrix, beta, selected, digits, endpoints, base, algebra)
                weights = {coefficients: carry_weight(
                    coefficients, selected, digits, endpoints, base, algebra)
                    for coefficients in product(range(modulus), repeat=size ** 2)}
                direct = algebra.zero
                for word, coefficients in reference:
                    weight = weights[coefficients]
                    if weight != algebra.zero:
                        phase = sum(frequency * value for frequency, value in zip(beta, word))
                        direct = algebra.add(direct, algebra.multiply(weight, algebra.phase(phase)))
                    direct_settings += 1
                assert direct == expected
                nonzero_cases += expected != algebra.zero
                disconnected_cases += len(segments(selected)) > 1
                query_cases += 1
                fiber_scalars.add((coefficient_prime,
                                  pow(modulus, count - len(selected), coefficient_prime)))

    # The root algebra can have zero divisors: Phi_3=(X-2)(X-4) over F_7.
    algebra = CharacterAlgebra(3, 7)
    assert algebra.multiply((5, 1), (3, 1)) == algebra.zero
    assert (5, 1) != algebra.zero and (3, 1) != algebra.zero

    # Surjectivity is essential. L(w)=(w,0) never reaches a point with c_1=1.
    # The unrestricted output sum with the false M^(r-t) factor is nonzero.
    nonsurjective = [[1], [0]]
    rank, _ = linear_solve(nonsurjective, [0, 0], 3)
    assert rank == 1 < 2
    actual_count = sum(1 for word in range(3) if (word, 0)[1] == 1)
    false_count_mod_five = pow(3, -1, 5) * 3 % 5
    assert actual_count == 0 and false_count_mod_five == 1

    return {
        "scope": "exact mixed linear/carry boundary queries; no full-network contraction claim",
        "surjective_coefficient_maps": surjective_cases,
        "compiled_vs_exhaustive_block_queries": query_cases,
        "row_space_negative_cases": outside_cases,
        "disconnected_carry_block_queries": disconnected_cases,
        "nonzero_block_values": nonzero_cases,
        "original_interpolation_input_settings_checked": direct_settings,
        "fiber_scalars": sorted(fiber_scalars),
        "negative_controls": {
            "coefficient_zero_divisors_allowed": True,
            "surjectivity_cannot_be_dropped": True,
        },
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
