#!/usr/bin/env python3
"""Exact controls for punctured sparse recovery and Frobenius digit boxes.

The finite checks supplement the ordinary proofs in the research note.
All fields are constructed with explicit binary defining polynomials;
the tensor-product checks also include an actual non-field algebra.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, product
import json
from math import gcd
from pathlib import Path
import random

from check_carry_fourier_elimination import (
    BinaryField, construct_field_and_root, irreducible,
)


def construct_field(degree: int) -> BinaryField:
    for lower in range(1, 1 << degree, 2):
        candidate = (1 << degree) | lower
        if irreducible(candidate):
            return BinaryField(candidate)
    raise AssertionError("no irreducible binary polynomial found")


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return True


def degree_cutoff(n: int, sparsity: int) -> int:
    h = n * (sparsity - 1) + 1
    a = h.bit_length()
    return 64 * (a + 1) * (1 << ((a + 1) // 2))


def degree_family(n: int, sparsity: int) -> list[int]:
    threshold = 1 << (n * (sparsity - 1))
    low = sparsity.bit_length()
    bound = degree_cutoff(n, sparsity)
    answer, accumulated = [], 1
    for degree in range(max(3, low + 1) | 1, bound + 1, 2):
        if is_prime(degree):
            order = (1 << degree) - 1
            assert order > sparsity
            assert gcd(order, accumulated) == 1
            answer.append(degree)
            accumulated *= order
            if accumulated > threshold:
                return answer
    raise AssertionError("the proved prime-degree cutoff was exceeded")


class TensorAlgebra:
    """F[Z]/g(Z), with g the binary defining polynomial of E.

    Tuples contain coefficients in F. This directly represents F tensor E
    and intentionally makes no irreducibility assumption over F.
    """

    def __init__(self, root_field: BinaryField, weight_field: BinaryField):
        self.root_field = root_field
        self.weight_field = weight_field
        self.degree = weight_field.degree
        self.zero = (0,) * self.degree
        self.one = (1,) + (0,) * (self.degree - 1)

    def add(self, left: tuple[int, ...], right: tuple[int, ...]):
        return tuple(a ^ b for a, b in zip(left, right))

    def mul(self, left: tuple[int, ...], right: tuple[int, ...]):
        degree, field = self.degree, self.root_field
        result = [0] * (2 * degree - 1)
        for i, a in enumerate(left):
            if a:
                for j, b in enumerate(right):
                    if b:
                        result[i + j] ^= field.mul(a, b)
        polynomial = self.weight_field.polynomial
        for i in range(2 * degree - 2, degree - 1, -1):
            value = result[i]
            if value:
                for j in range(degree):
                    if (polynomial >> j) & 1:
                        result[i - degree + j] ^= value
        return tuple(result[:degree])

    def power(self, value: tuple[int, ...], exponent: int):
        result = self.one
        while exponent:
            if exponent & 1:
                result = self.mul(result, value)
            value = self.mul(value, value)
            exponent >>= 1
        return result

    def root_scalar(self, value: int):
        return (value,) + (0,) * (self.degree - 1)

    def weight_scalar(self, value: int):
        return tuple((value >> i) & 1 for i in range(self.degree))

    def coefficient_frobenius(self, value: tuple[int, ...], width: int):
        return tuple(self.root_field.power(a, 1 << width) for a in value)

    def geometric(self, value: tuple[int, ...], length: int):
        total, power = self.zero, self.one
        step_total, step_power = self.one, value
        while length:
            if length & 1:
                total = self.add(total, self.mul(power, step_total))
                power = self.mul(power, step_power)
            step_total = self.add(step_total, self.mul(step_power, step_total))
            step_power = self.mul(step_power, step_power)
            length >>= 1
        return total


def binary_coefficient_evaluate(field: BinaryField, polynomial, x: int, y: int):
    """Two coordinates model F4 or F2[epsilon]/epsilon^2 coefficients.

    The scalar extension is free in either coefficient basis. Polynomial
    evaluation only adds scalar multiples of coefficients, so this exact
    representation tests nonbinary and zero-divisor coefficients alike.
    """
    answer = [0, 0]
    for coefficient, u, v in polynomial:
        monomial = field.mul(field.power(x, u), field.power(y, v))
        for bit in range(2):
            if (coefficient >> bit) & 1:
                answer[bit] ^= monomial
    return tuple(answer)


def reduced_polynomial(polynomial, order: int):
    coefficients = {}
    for coefficient, u, v in polynomial:
        exponent = (u % order, v % order)
        coefficients[exponent] = coefficients.get(exponent, 0) ^ coefficient
    return {exponent: value for exponent, value in coefficients.items() if value}


def anchor_product(polynomial):
    _, u0, v0 = polynomial[0]
    answer = 1
    for _, u, v in polynomial[1:]:
        difference = abs(u - u0) if u != u0 else abs(v - v0)
        assert difference > 0
        answer *= difference
    return answer


def check_sparse_recovery(counts: Counter):
    rng = random.Random(16012026)
    fields = {degree: construct_field(degree) for degree in (3, 5, 7)}
    models = []
    support = [(0, 0), (0, 7), (7, 0), (7, 7), (1, 2), (2, 1)]
    for pair in combinations(support, 2):
        for coefficients in product((1, 2, 3), repeat=2):
            models.append([(coefficient, *exponent)
                           for coefficient, exponent in zip(coefficients, pair)])
    for _ in range(90):
        support3 = rng.sample([(u, v) for u in range(8) for v in range(8)], 3)
        models.append([(rng.randrange(1, 4), *exponent) for exponent in support3])
    models.extend([
        [(2, 0, 0), (2, 7, 7)],
        [(3, 0, 0), (1, 7, 0), (2, 0, 7)],
        [(1, 0, 0)], [(2, 3, 4)], [(3, 7, 7)],
    ])
    for polynomial in models:
        sparsity = len(polynomial)
        family = degree_family(3, sparsity)
        hit = False
        difference_product = anchor_product(polynomial)
        assert difference_product <= 1 << (3 * (sparsity - 1))
        bad_order_product = 1
        for degree in family:
            field = fields[degree]
            order = field.size - 1
            reduced = reduced_polynomial(polynomial, order)
            full_grid_zero = True
            for x, y in product(range(2, field.size), repeat=2):
                value = binary_coefficient_evaluate(field, polynomial, x, y)
                counts["punctured_grid_values"] += 1
                if value != (0, 0):
                    full_grid_zero = False
                    if degree != 3:
                        break
            assert full_grid_zero == (not reduced)
            counts["punctured_interpolation_equivalences"] += 1
            if full_grid_zero:
                assert difference_product % order == 0
                bad_order_product *= order
                assert difference_product % bad_order_product == 0
                counts["anchor_collision_divisibilities"] += 1
            else:
                hit = True
        assert hit
        counts["nonbinary_sparse_polynomials"] += 1

    # One entire small field is deliberately bad; the next degree separates.
    cancellation = [(2, 0, 0), (2, 7, 7)]
    assert not reduced_polynomial(cancellation, 7)
    assert all(binary_coefficient_evaluate(fields[3], cancellation, x, y) == (0, 0)
               for x, y in product(range(2, 8), repeat=2))
    assert any(binary_coefficient_evaluate(fields[5], cancellation, x, y) != (0, 0)
               for x, y in product(range(2, 32), repeat=2))
    counts["deliberate_whole_field_cancellations"] += 1

    # Sparsity equal to the group order invalidates the punctured conclusion.
    boundary = [(2, exponent, 0) for exponent in range(7)]
    assert all(binary_coefficient_evaluate(fields[3], boundary, x, 2) == (0, 0)
               for x in range(2, 8))
    assert binary_coefficient_evaluate(fields[3], boundary, 1, 2) == (0, 1)
    assert reduced_polynomial(boundary, 7)
    counts["strict_sparsity_boundary_controls"] += 1

    # Each degree may change every coefficient, provided the exponent
    # support stays fixed and all coefficients remain nonzero.
    for _ in range(50):
        support = rng.sample([(u, v) for u in range(8) for v in range(8)], 3)
        anchor = anchor_product([(1, *exponent) for exponent in support])
        hit = False
        for degree in degree_family(3, 3):
            field = fields[degree]
            coefficients = [rng.randrange(1, field.size) for _ in support]
            reduced = {}
            for coefficient, (u, v) in zip(coefficients, support):
                key = (u % (field.size - 1), v % (field.size - 1))
                reduced[key] = reduced.get(key, 0) ^ coefficient
            reduced = {key: coefficient for key, coefficient in reduced.items() if coefficient}
            found = False
            for x, y in product(range(2, field.size), repeat=2):
                value = 0
                for coefficient, (u, v) in zip(coefficients, support):
                    value ^= field.product(coefficient, field.power(x, u), field.power(y, v))
                counts["degree_dependent_grid_values"] += 1
                if value:
                    found = True
                    break
            assert found == bool(reduced)
            if not found:
                assert anchor % (field.size - 1) == 0
            hit |= found
        assert hit
        counts["degree_dependent_coefficient_supports"] += 1

    samples = []
    for n in (1, 2, 3, 8, 16, 32, 64):
        for sparsity in (1, 2, 3, 4, 7, 8, 15, 16):
            family = degree_family(n, sparsity)
            accumulated = 1
            for degree in family:
                assert is_prime(degree) and degree % 2
                order = (1 << degree) - 1
                assert gcd(order, accumulated) == 1
                accumulated *= order
            assert accumulated > 1 << (n * (sparsity - 1))
            assert family[-1] <= degree_cutoff(n, sparsity)
            counts["prime_degree_product_certificates"] += 1
            if n in (1, 8, 64) and sparsity in (1, 3, 16):
                samples.append({"n": n, "K": sparsity, "degrees": family,
                                "proved_cutoff": degree_cutoff(n, sparsity),
                                "product_bit_length": accumulated.bit_length()})
    return samples


def binary_polynomial_value(field: BinaryField, polynomial: int, value: int):
    answer = 0
    for i in range(polynomial.bit_length() - 1, -1, -1):
        answer = field.mul(answer, value) ^ ((polynomial >> i) & 1)
    return answer


def check_prefix_recovery(counts: Counter):
    fields = {degree: construct_field(degree) for degree in (3, 5, 7)}
    records = []
    for initial in (2, 3, 4, 5, 6, 7, 9, 11, 13, 17, 25, 49):
        value, recovered = initial, []
        while value > 1:
            n = value.bit_length()
            divisors = [d for d in range(1, value + 1) if value % d == 0]
            family = degree_family(n, len(divisors))
            assert set(family) <= set(fields)

            def prefix_nonzero(k, prefix):
                terms = [(1, d, value // d) for d in divisors
                         if d != 1 and d >> k == prefix]
                counts["bivariate_prefix_queries"] += 1
                for degree in family:
                    field = fields[degree]
                    for x, y in product(range(2, field.size), repeat=2):
                        answer = binary_coefficient_evaluate(field, terms, x, y)
                        counts["prefix_field_evaluations"] += 1
                        if answer != (0, 0):
                            assert terms
                            return True
                assert not terms
                return False

            prefix = 0
            for k in range(n - 1, -1, -1):
                prefix *= 2
                if not prefix_nonzero(k, prefix):
                    prefix += 1
            assert prefix == min(d for d in divisors if d > 1)
            assert is_prime(prefix)
            recovered.append(prefix)
            value //= prefix
        check = 1
        for prime in recovered:
            check *= prime
        assert check == initial
        counts["complete_prefix_factorizations"] += 1
        records.append({"N": initial, "recovered_primes": recovered})
    return records


def exhibit_zero_divisors(algebra: TensorAlgebra, counts: Counter):
    field, degree = algebra.root_field, algebra.degree
    assert field.degree % degree == 0
    order = (1 << degree) - 1
    for candidate in range(2, field.size):
        generator = field.power(candidate, (field.size - 1) // order)
        if generator != 1:
            break
    else:
        raise AssertionError("missing subfield generator")
    root = None
    for exponent in range(1, order):
        value = field.power(generator, exponent)
        if binary_polynomial_value(field, algebra.weight_field.polynomial, value) == 0:
            root = value
            break
    assert root is not None
    linear = (root, 1) + (0,) * (degree - 2)
    # Synthetic division of the binary polynomial by Z+root.
    quotient = [0] * degree
    quotient[-1] = 1
    polynomial = algebra.weight_field.polynomial
    for i in range(degree - 1, 0, -1):
        quotient[i - 1] = ((polynomial >> i) & 1) ^ field.mul(root, quotient[i])
    assert (polynomial & 1) ^ field.mul(root, quotient[0]) == 0
    quotient = tuple(quotient)
    assert linear != algebra.zero and quotient != algebra.zero
    assert algebra.mul(linear, quotient) == algebra.zero
    counts["explicit_tensor_zero_divisors"] += 1
    return {"linear_factor": list(linear), "nonzero_annihilator": list(quotient)}


def check_tensor_boxes(counts: Counter):
    rng = random.Random(31012026)
    records = []
    for degree, width, m in ((3, 3, 8), (3, 3, 12), (5, 5, 12)):
        field, root, construction = construct_field_and_root(m)
        weights = construct_field(degree)
        algebra = TensorAlgebra(field, weights)
        modulus, base = (1 << m) + 1, 1 << width
        assert width % degree == 0
        assert gcd(modulus, weights.size - 1) == 1
        root_element = algebra.root_scalar(root)
        assert algebra.power(root_element, modulus) == algebra.one
        zero_divisors = None
        if field.degree % degree == 0:
            zero_divisors = exhibit_zero_divisors(algebra, counts)
        for _ in range(8):
            value = tuple(rng.randrange(field.size) for _ in range(degree))
            assert algebra.power(value, base) == algebra.coefficient_frobenius(value, width)
            counts["tensor_frobenius_actions"] += 1
        frequencies = sorted({0, 1, 2, 3, modulus // 2, modulus - 2, modulus - 1,
                              *(rng.randrange(modulus) for _ in range(7))})
        allowed_weights = list(range(2, weights.size))
        if len(allowed_weights) > 10:
            allowed_weights = [2, 3, weights.size - 1] + rng.sample(range(4, weights.size - 1), 7)
        for alpha in allowed_weights:
            scalar = algebra.weight_scalar(alpha)
            assert algebra.power(scalar, base) == scalar
            for position in range(4):
                assert weights.power(alpha, base ** position) == alpha
                counts["position_independent_digit_weights"] += 1
            denominator = 1 ^ weights.power(alpha, modulus)
            assert denominator != 0
            scalar_inverse = algebra.weight_scalar(weights.inverse(denominator))
            for frequency in frequencies:
                value = algebra.mul(scalar, algebra.power(root_element, frequency))
                unit = algebra.add(algebra.one, value)
                inverse = algebra.mul(algebra.geometric(value, modulus), scalar_inverse)
                assert algebra.mul(unit, inverse) == algebra.one
                counts["explicit_mixed_character_inverses"] += 1
                box = algebra.geometric(value, base)
                literal, power = algebra.zero, algebra.one
                for _ in range(base):
                    literal = algebra.add(literal, power)
                    power = algebra.mul(power, value)
                assert box == literal
                numerator = algebra.power(unit, base)
                assert numerator == algebra.coefficient_frobenius(unit, width)
                assert box == algebra.mul(numerator, inverse)
                numerator_inverse = algebra.power(inverse, base)
                box_inverse = algebra.mul(unit, numerator_inverse)
                assert algebra.mul(box, box_inverse) == algebra.one
                counts["full_box_coboundaries"] += 1
                for partial_width in range(width + 1):
                    length = 1 << partial_width
                    offset = (3 * frequency + alpha) % base
                    phase = algebra.power(value, offset)
                    phase_inverse = algebra.power(value, modulus * (weights.size - 1) - offset)
                    partial = algebra.mul(phase, algebra.geometric(value, length))
                    claimed = algebra.mul(phase, algebra.mul(algebra.power(unit, length), inverse))
                    assert partial == claimed
                    partial_inverse = algebra.mul(
                        phase_inverse, algebra.mul(unit, algebra.power(inverse, length)))
                    assert algebra.mul(partial, partial_inverse) == algebra.one
                    counts["dyadic_interval_unit_formulas"] += 1
        # Restoring alpha=1 produces the genuine zero-frequency failure.
        assert algebra.geometric(algebra.one, base) == algebra.zero
        counts["excluded_unit_weight_negative_controls"] += 1
        records.append({**construction, "weight_degree": degree,
                        "weight_polynomial_binary": bin(weights.polynomial),
                        "digit_width": width, "tested_frequencies": frequencies,
                        "tested_nontrivial_weights": allowed_weights,
                        "tensor_binary_dimension": 2 * m * degree,
                        "zero_divisor_witness": zero_divisors})
    for m in range(1, 65):
        for degree in (3, 5, 7, 11, 13, 17, 19, 23):
            assert gcd((1 << m) + 1, (1 << degree) - 1) == 1
            counts["coprime_character_and_weight_orders"] += 1
    return records


def exact_pretest_width(n: int):
    """Use beta=3/4 and the optional 100n/log2(n) term, with integers."""
    q = n.bit_length() - 1
    low, high = 0, n
    target = n ** 4
    while low + 1 < high:
        middle = (low + high) // 2
        if middle ** 4 * q ** 3 >= target:
            high = middle
        else:
            low = middle
    assert high ** 4 * q ** 3 >= target
    assert (high - 1) ** 4 * q ** 3 < target
    return max(high, (100 * n + q - 1) // q)


def check_integer_parameters(counts: Counter):
    records = []
    for q in (128, 256, 512, 1024, 2048, 4096):
        for increment in (0, 1, 17):
            n = (1 << q) + increment
            u = exact_pretest_width(n)
            sparsity = 1 << ((n - 1) // u)
            assert sparsity ** 100 <= n
            cutoff = degree_cutoff(n, sparsity)
            p = 1 << (q // 40)
            h = p ** 9
            d = 2 * h * h
            assert p ** 40 <= n
            assert d < 8 * p ** 18
            w0 = (n + d - 1) // d
            usable = d * cutoff * 16 <= n
            if q >= 512:
                assert usable
            for ell in (3, 5, 7, 11, 127, cutoff - 1, cutoff):
                rounded = ell * ((w0 + ell - 1) // ell)
                assert w0 <= rounded < w0 + ell
                assert rounded % ell == 0
                assert d * (rounded - w0) < d * ell <= d * cutoff
                t = d.bit_length() - 1
                m = d * ((2 * rounded + t + d - 1) // d)
                assert m % d == 0
                assert 2 * rounded + t <= m < 2 * rounded + t + d
                # Exponent comparisons certify the coefficient and root bounds
                # without ever constructing 2**rounded or 2**m.
                assert (2 * rounded + t) <= m
                counts["exact_width_and_modulus_roundings"] += 1
            counts["roughness_and_capacity_parameters"] += 1
            if increment == 0:
                records.append({"log2_n": q, "pretest_width_bit_length": u.bit_length(),
                                "K": sparsity, "degree_cutoff_bit_length": cutoff.bit_length(),
                                "log2_P": p.bit_length() - 1,
                                "log2_D": d.bit_length() - 1,
                                "DL_bit_length": (d * cutoff).bit_length(),
                                "capacity_check_DL_at_most_n_over_16": usable})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts = Counter()
    families = check_sparse_recovery(counts)
    factorizations = check_prefix_recovery(counts)
    tensors = check_tensor_boxes(counts)
    parameters = check_integer_parameters(counts)
    result = {"status": "all exact controls passed",
              "claim": "Punctured sparse recovery and everywhere-unit weighted Frobenius boxes",
              "score_improvement_claimed": False,
              "summary": dict(sorted(counts.items())),
              "prime_degree_family_samples": families,
              "prefix_factorization_samples": factorizations,
              "tensor_algebras": tensors,
              "integer_parameter_samples": parameters}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "summary": result["summary"]}, indent=2))


if __name__ == "__main__":
    main()
