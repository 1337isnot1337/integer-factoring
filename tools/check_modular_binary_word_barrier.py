#!/usr/bin/env python3
"""Exact finite certificates for the characteristic-two binary word barrier."""

from __future__ import annotations

import argparse
from collections import defaultdict
from functools import reduce
import hashlib
from itertools import combinations_with_replacement, permutations, product
import json
from math import comb, factorial, prod
from operator import xor
from pathlib import Path


def bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def insert(basis, vector):
    while vector:
        pivot = vector.bit_length() - 1
        if pivot in basis:
            vector ^= basis[pivot]
        else:
            basis[pivot] = vector
            return True
    return False


def multiply(left, right):
    answer = set()
    for a, b in product(left, right):
        monomial = tuple(sorted(a + b))
        if monomial in answer:
            answer.remove(monomial)
        else:
            answer.add(monomial)
    return answer


def add(*polynomials):
    answer = set()
    for polynomial in polynomials:
        answer.symmetric_difference_update(polynomial)
    return answer


def minor(rows, columns):
    return {
        tuple(sorted(4 * r + c for r, c in zip(rows, permutation)))
        for permutation in permutations(columns)
    }


def adjoint_coefficient(i, j, a, b):
    """A[i,a] adj(A)[b,j] in characteristic two."""
    cofactor = minor(
        [r for r in range(4) if r != j],
        [c for c in range(4) if c != b],
    )
    return multiply({(4 * i + a,)}, cofactor)


def simple_adjoint_coefficients():
    inputs = [[(a, b)] for a, b in permutations(range(4), 2)]
    inputs += [[(0, 0), (3, 3)], [(1, 1), (3, 3)]]
    outputs = [[(i, j)] for i, j in permutations(range(4), 2)]
    outputs += [[(0, 0), (2, 2)], [(1, 1), (2, 2)]]
    return [
        add(*(adjoint_coefficient(i, j, a, b)
              for a, b in input_entries for i, j in output_entries))
        for input_entries in inputs for output_entries in outputs
    ]


def reshuffle(polynomial):
    answer = set()
    for monomial in polynomial:
        new = []
        for variable in monomial:
            x, y, u, v = ((variable >> i) & 1 for i in (3, 2, 1, 0))
            new.append(8 * x + 4 * u + 2 * y + v)
        answer.add(tuple(sorted(new)))
    return answer


def evaluate_crossing(polynomial):
    return sum(all(v in (0, 6, 9, 15) for v in monomial)
               for monomial in polynomial) & 1


def intersection(left, right):
    basis = {}
    answer = {}
    for vector, provenance in [(p, 0) for p in left] + [(p, p) for p in right]:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = vector, provenance
                break
            other, other_provenance = basis[pivot]
            vector ^= other
            provenance ^= other_provenance
        if not vector and provenance:
            insert(answer, provenance)
    return answer


def polynomial_checks():
    determinant = minor(range(4), range(4))
    trace = add(*(adjoint_coefficient(i, i, a, a)
                  for i in (0, 1) for a in (1, 3)))
    witness = add(determinant, trace)
    factored = add(
        multiply(minor((0, 1), (0, 2)), minor((2, 3), (1, 3))),
        multiply(minor((0, 1), (1, 3)), minor((2, 3), (0, 2))),
    )
    assert witness == factored
    assert len(witness) == 8
    assert reshuffle(witness) == witness
    assert all(reshuffle({monomial}) == {monomial} for monomial in witness)
    assert evaluate_crossing(witness) == 1

    coefficients = simple_adjoint_coefficients()
    # X=diag(0,1,0,1), ell(Y)=Y00+Y11. The chosen quotient
    # coordinates express ell as the sum of its last two output forms.
    assert trace == add(coefficients[14 * 13 + 12],
                        coefficients[14 * 13 + 13])

    monomials = list(combinations_with_replacement(range(16), 4))
    index = {m: i for i, m in enumerate(monomials)}

    def encode(polynomial):
        return sum(1 << index[m] for m in polynomial)

    coefficient_bits = [encode(p) for p in coefficients]
    basis = {}
    for vector in coefficient_bits:
        insert(basis, vector)
    assert len(basis) == 196
    adjoint_common = intersection(
        coefficient_bits, [encode(reshuffle(p)) for p in coefficients]
    )
    assert len(adjoint_common) == 17
    evaluation_mask = sum(
        1 << i for i, monomial in enumerate(monomials)
        if all(v in (0, 6, 9, 15) for v in monomial)
    )
    assert all((p & evaluation_mask).bit_count() % 2 == 0
               for p in adjoint_common.values())
    assert insert(basis, encode(determinant))
    assert len(basis) == 197
    all_polynomials = coefficients + [determinant]
    common = intersection(
        [encode(p) for p in all_polynomials],
        [encode(reshuffle(p)) for p in all_polynomials],
    )
    assert len(common) == 18
    common_basis = common.copy()
    assert not insert(common_basis, encode(witness))

    # Frobenius powers are checked as exact monomial multisets.
    power = witness
    for j in range(1, 7):
        power = multiply(power, power)
        expected = {tuple(v for v in m for _ in range(1 << j)) for m in witness}
        assert power == expected
        assert reshuffle(power) == power
        assert evaluate_crossing(power) == 1

    return {
        "witness_monomials": [list(m) for m in sorted(witness)],
        "simple_adjoint_coefficient_dimension": len(coefficients),
        "simple_adjoint_common_dimension": len(adjoint_common),
        "combined_common_dimension": len(common),
        "frobenius_powers_checked": 7,
    }


def dimension_vandermonde(shape):
    shape = tuple(shape) + (0,) * (4 - len(shape))
    total = sum(shape)
    numerator = factorial(total) * prod(
        shape[i] - shape[j] + j - i
        for i in range(4) for j in range(i + 1, 4)
    )
    denominator = prod(factorial(shape[i] + 3 - i) for i in range(4))
    assert numerator % denominator == 0
    return numerator // denominator


def dimension_hooks(shape):
    hooks = [
        shape[r] - c + sum(other > c for other in shape[r + 1:])
        for r in range(len(shape)) for c in range(shape[r])
    ]
    denominator = prod(hooks)
    numerator = factorial(sum(shape))
    assert numerator % denominator == 0
    return numerator // denominator


def tableaux_words(shape):
    counts = [0] * len(shape)
    word = []

    def visit():
        if len(word) == sum(shape):
            yield tuple(word)
            return
        for row in range(len(shape)):
            if counts[row] < shape[row] and (row == 0 or counts[row] < counts[row - 1]):
                counts[row] += 1
                word.append(row)
                yield from visit()
                word.pop()
                counts[row] -= 1

    yield from visit()


def encode_word(word):
    answer = 0
    for letter in word:
        answer = 4 * answer + letter
    return answer


def decode_word(value, length):
    return tuple((value >> (2 * (length - 1 - i))) & 3 for i in range(length))


def tableau_vector(word, shape):
    rows = [[i for i, r in enumerate(word) if r == row]
            for row in range(len(shape))]
    columns = [
        [rows[row][column] for row in range(len(shape)) if column < shape[row]]
        for column in range(max(shape))
    ]
    words = {word}
    for positions in columns:
        assert positions == sorted(positions)
        new = set()
        for current in words:
            for labels in permutations(range(len(positions))):
                changed = list(current)
                for position, label in zip(positions, labels):
                    changed[position] = label
                new.add(tuple(changed))
        words = new
    vector = sum(1 << encode_word(w) for w in words)
    assert min(words) == word
    return vector


def upper_distribution(vector, length, source, target):
    """Coefficients in s for e_source -> e_source+s e_target."""
    coefficients = defaultdict(int)
    for index in bits(vector):
        word = decode_word(index, length)
        positions = [i for i, letter in enumerate(word) if letter == source]
        for mask in range(1 << len(positions)):
            changed = list(word)
            for j, position in enumerate(positions):
                if mask >> j & 1:
                    changed[position] = target
            coefficients[mask.bit_count()] ^= 1 << encode_word(changed)
    return coefficients


def tableaux_checks():
    vector_count = distribution_count = 0
    for k in (1, 2):
        for shape in ((2 * k, k, k, 0), (k, k, k, k)):
            basis = {}
            words = list(tableaux_words(shape))
            assert len(words) == dimension_hooks(shape)
            for word in words:
                vector = tableau_vector(word, shape)
                assert insert(basis, vector)
                for source in range(1, 4):
                    for target in range(source):
                        coefficients = upper_distribution(vector, len(word), source, target)
                        assert coefficients[0] == vector
                        assert all(not value for degree, value in coefficients.items() if degree)
                        distribution_count += 1
                vector_count += 1
            assert len(basis) == dimension_vandermonde(shape)
    return vector_count, distribution_count


def dimension_checks():
    cases = 0
    for k in range(1, 129):
        multinomial = comb(4 * k, 2 * k) * comb(2 * k, k)
        assert multinomial * (4 * k + 1) * (2 * k + 1) >= 1 << (6 * k)
        for shape in ((2 * k, k, k, 0), (k, k, k, k)):
            dimension = dimension_vandermonde(shape)
            assert dimension == dimension_hooks(shape)
            assert dimension * (4 * k + 4)**6 >= multinomial
            assert dimension * (4 * k + 4)**8 >= 1 << (6 * k)
            cases += 1
    return cases


def position_class_checks():
    cases = 0
    for length in range(4, 16385):
        blocks = [4 * (1 << j) for j in bits(length // 4)]
        retained = sum(blocks)
        assert retained == 4 * (length // 4)
        assert length - 3 <= retained <= length
        assert all(block % 4 == 0 and (block // 4).bit_count() == 1 for block in blocks)
        denominator = prod((block + 4)**8 for block in blocks)
        assert denominator <= (length + 4)**(8 * len(blocks))
        assert len(blocks) <= length.bit_length()
        assert sum(3 * block // 2 for block in blocks) == 3 * retained // 2
        cases += 1
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    polynomial = polynomial_checks()
    tableaux, distributions = tableaux_checks()
    hook_cases = dimension_checks()
    position_cases = position_class_checks()
    # rho=3/4>1/sqrt(2), or equivalently log2(4/3)<1/2.
    assert 9 > 8 and 16 < 18
    result = {
        "status": "passed",
        "scope": "Exact finite certificates supporting the written modular all-length rank proof; no new factoring exponent.",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "field_characteristic": 2,
        **polynomial,
        "independent_standard_tableau_vectors": tableaux,
        "exact_upper_unipotent_distribution_checks": distributions,
        "independent_hook_formula_and_integer_bounds": hook_cases,
        "binary_position_class_cases": position_cases,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
