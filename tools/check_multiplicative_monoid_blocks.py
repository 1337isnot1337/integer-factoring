#!/usr/bin/env python3
"""Exact modular-character blocks and Fourier--Mellin closure controls.

The defining coefficient fields, unit generators, characters, polynomial
blocks, and interpolation algorithms are all explicit.  Characteristic-two
unit-group radicals are retained rather than replaced by scalar characters.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from math import gcd, lcm, prod
from pathlib import Path
import random

from check_carry_fourier_elimination import BinaryField, field_rank, prime_factors
from check_frobenius_sparse_recovery import construct_field


def factor_powers(value: int) -> list[tuple[int, int]]:
    result = []
    for prime in prime_factors(value):
        exponent = 0
        while value % prime == 0:
            value //= prime
            exponent += 1
        result.append((prime, exponent))
    assert value == 1
    return result


def valuation(value: int, prime: int, ceiling: int) -> int:
    if value == 0:
        return ceiling
    result = 0
    while value % prime == 0 and result < ceiling:
        value //= prime
        result += 1
    return result


def odd_part(value: int) -> int:
    return value >> valuation(value, 2, value.bit_length())


def totient_prime_power(prime: int, exponent: int) -> int:
    return 1 if exponent == 0 else (prime - 1) * prime ** (exponent - 1)


def primitive_unit(prime: int, exponent: int) -> int:
    modulus = prime ** exponent
    order = totient_prime_power(prime, exponent)
    for candidate in range(1, modulus):
        if gcd(candidate, modulus) == 1 and all(
                pow(candidate, order // q, modulus) != 1
                for q in prime_factors(order)):
            return candidate
    raise AssertionError("odd prime-power unit group has no generator")


def field_element_of_order(field: BinaryField, order: int) -> int:
    assert (field.size - 1) % order == 0
    for candidate in range(1, field.size):
        value = field.power(candidate, (field.size - 1) // order)
        if all(field.power(value, order // p) != 1 for p in prime_factors(order)):
            assert field.power(value, order) == 1
            return value
    raise AssertionError("missing prescribed field root")


def coefficient_field(modulus: int, additive_root: bool = False):
    odd_order = 1
    for prime, exponent in factor_powers(modulus):
        odd_order = lcm(odd_order, odd_part(totient_prime_power(prime, exponent)))
    required_order = lcm(odd_order, modulus) if additive_root else odd_order
    degree = 1
    while pow(2, degree, required_order) != 1 % required_order:
        degree += 1
    least_degree = (2 * modulus).bit_length()
    degree *= (least_degree + degree - 1) // degree
    field = construct_field(degree)
    root = field_element_of_order(field, required_order)
    odd_root = field.power(root, required_order // odd_order)
    additive = field.power(root, required_order // modulus) if additive_root else None
    return field, odd_order, odd_root, additive


def binomial_bit(top: int, bottom: int) -> int:
    return int(bottom & top == bottom)


def poly_mul(field: BinaryField, left, right, length: int):
    result = [0] * length
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right[:max(0, length - i)]):
                if b:
                    result[i + j] ^= field.mul(a, b)
    return tuple(result)


def poly_power(field: BinaryField, value, exponent: int, length: int):
    result = (1,) + (0,) * (length - 1)
    while exponent:
        if exponent & 1:
            result = poly_mul(field, result, value, length)
        value = poly_mul(field, value, value, length)
        exponent >>= 1
    return result


def block_indices(shape):
    return tuple(product(*(range(length) for length in shape)))


def block_mul(field: BinaryField, shape, left, right):
    indices = block_indices(shape)
    positions = {index: i for i, index in enumerate(indices)}
    result = [0] * len(indices)
    for i, a in enumerate(left):
        if not a:
            continue
        for j, b in enumerate(right):
            if not b:
                continue
            index = tuple(x + y for x, y in zip(indices[i], indices[j]))
            if all(x < length for x, length in zip(index, shape)):
                result[positions[index]] ^= field.mul(a, b)
    return tuple(result)


class PrimePowerBlocks:
    def __init__(self, prime, exponent, field, common_order, common_root):
        self.prime, self.exponent = prime, exponent
        self.modulus, self.field = prime ** exponent, field
        self.generator = primitive_unit(prime, exponent)
        order = totient_prime_power(prime, exponent)
        self.h = odd_part(prime - 1)
        self.d = (prime - 1) // self.h
        self.odd_order = order // self.d
        self.root = field.power(common_root, common_order // self.odd_order)
        self.logs = {pow(self.generator, k, self.modulus): k for k in range(order)}
        self.lengths = [1 + min(valuation(j, prime, exponent - 1), exponent - 1)
                        for j in range(self.odd_order)]
        self.shapes = [()] + [(self.d, length) for length in self.lengths]
        assert sum(prod(shape) for shape in self.shapes) == self.modulus

    def image(self, residue):
        result = [(1,)]
        if residue == 0:
            return result + [(0,) * prod(shape) for shape in self.shapes[1:]]
        level = valuation(residue, self.prime, self.exponent)
        unit = residue // self.prime ** level
        log = self.logs[unit]
        even = log * pow(self.odd_order, -1, self.d) % self.d
        odd = (log * pow(self.d, -1, self.odd_order) % self.odd_order
               if self.odd_order > 1 else 0)
        for j, length in enumerate(self.lengths):
            scalar = self.field.power(self.root, j * odd)
            result.append(tuple(scalar if y == level and binomial_bit(even, x) else 0
                                for x, y in block_indices((self.d, length))))
        return result


class MonoidBlocks:
    def __init__(self, modulus):
        self.modulus = modulus
        self.field, odd_order, root, _ = coefficient_field(modulus)
        self.local = [PrimePowerBlocks(p, e, self.field, odd_order, root)
                      for p, e in factor_powers(modulus)]
        self.block_choices = tuple(product(*(range(len(part.shapes))
                                             for part in self.local)))
        self.shapes = [sum((part.shapes[index] for part, index in zip(self.local, choice)), ())
                       for choice in self.block_choices]
        assert sum(prod(shape) for shape in self.shapes) == modulus
        self.images = [self.image(a) for a in range(modulus)]

    def image(self, residue):
        local_images = [part.image(residue % part.modulus) for part in self.local]
        result = []
        for choice in self.block_choices:
            coefficient_lists = [image[index] for image, index in zip(local_images, choice)]
            result.append(tuple(self.field.product(*coefficients)
                                for coefficients in product(*coefficient_lists)))
        return result

    def multiply(self, left, right):
        return [block_mul(self.field, shape, a, b)
                for shape, a, b in zip(self.shapes, left, right)]

    def transform(self, values):
        result = [[0] * prod(shape) for shape in self.shapes]
        for value, image in zip(values, self.images):
            if value:
                for out, block in zip(result, image):
                    for j, entry in enumerate(block):
                        out[j] ^= self.field.mul(value, entry)
        return [tuple(block) for block in result]


def evaluate_polynomial(field, coefficients, point):
    value = 0
    for coefficient in reversed(coefficients):
        value = field.mul(value, point) ^ coefficient
    return value


def interpolation_weights(field, length):
    """Coefficient columns of the Lagrange basis at the nodes 0,...,length-1."""
    assert length <= field.size
    columns = []
    for point in range(length):
        numerator, denominator = [1], 1
        for other in range(length):
            if point == other:
                continue
            denominator = field.mul(denominator, point ^ other)
            result = [0] * (len(numerator) + 1)
            for i, value in enumerate(numerator):
                result[i] ^= field.mul(value, other)
                result[i + 1] ^= value
            numerator = result
        inverse = field.inverse(denominator)
        columns.append(tuple(field.mul(value, inverse) for value in numerator))
    return columns


def grid_product(field, shape, left, right):
    """Literal bilinear interpolation of the untruncated product, then truncate."""
    indices = block_indices(shape)
    grids = tuple(product(*(range(2 * length - 1) for length in shape)))
    weights = [interpolation_weights(field, 2 * length - 1) for length in shape]
    output = [0] * len(indices)
    for point in grids:
        def evaluate(values):
            total = 0
            for coefficient, index in zip(values, indices):
                total ^= field.product(coefficient, *(
                    field.power(x, power) for x, power in zip(point, index)))
            return total
        product_value = field.mul(evaluate(left), evaluate(right))
        for i, index in enumerate(indices):
            output[i] ^= field.product(product_value, *(
                weight[x][power] for weight, x, power in zip(weights, point, index)))
    return tuple(output), len(grids)


def local_product_count(prime, exponent):
    h = odd_part(prime - 1)
    d = (prime - 1) // h
    return 1 + h * (2 * d - 1) * sum(
        totient_prime_power(prime, r) * (2 * (exponent - r) - 1)
        for r in range(exponent))


class PrimeMoments:
    def __init__(self, prime):
        self.prime = prime
        self.field, order, root, self.zeta = coefficient_field(prime, additive_root=True)
        self.local = PrimePowerBlocks(prime, 1, self.field, order, root)
        self.d, self.h = self.local.d, self.local.odd_order
        self.u = (1, 1) + (0,) * (self.d - 2)
        self.images = [self.local.image(a)[1:] for a in range(prime)]
        self.gauss = self.moments([self.field.power(self.zeta, a) for a in range(prime)])[1]

    def moments(self, values):
        blocks = [[0] * self.d for _ in range(self.h)]
        for a in range(1, self.prime):
            for block, image in zip(blocks, self.images[a]):
                for k, value in enumerate(image):
                    block[k] ^= self.field.mul(values[a], value)
        return values[0], tuple(tuple(block) for block in blocks)

    def iota(self, values):
        # U^{-1}=U^{d-1}, and X maps to U^{-1}+1.
        transformed_x = list(poly_power(self.field, self.u, self.d - 1, self.d))
        transformed_x[0] ^= 1
        answer, power = [0] * self.d, (1,) + (0,) * (self.d - 1)
        for value in values:
            for k, coefficient in enumerate(power):
                answer[k] ^= self.field.mul(value, coefficient)
            power = poly_mul(self.field, power, transformed_x, self.d)
        return tuple(answer)

    def fourier_formula(self, moments):
        zero, blocks = moments
        output = []
        for character in range(self.h):
            value = list(poly_mul(self.field, self.gauss[character],
                                  self.iota(blocks[-character % self.h]), self.d))
            if character == 0:
                value[-1] ^= zero
            output.append(tuple(value))
        return zero ^ blocks[0][0], tuple(output)

    def pointwise_formula(self, left, right):
        left_zero, left_blocks = left
        right_zero, right_blocks = right
        result = [[0] * self.d for _ in range(self.h)]
        for a, b in product(range(self.h), repeat=2):
            output = result[(a + b) % self.h]
            for i, j in product(range(self.d), repeat=2):
                output[i & j] ^= self.field.mul(left_blocks[a][i], right_blocks[b][j])
        return self.field.mul(left_zero, right_zero), tuple(tuple(block) for block in result)

    def product_formula(self, left, right):
        a0, a = left
        b0, b = right
        zero = self.field.product(a0, b0) ^ self.field.product(a0, b[0][0])
        zero ^= self.field.product(b0, a[0][0])
        blocks = tuple(poly_mul(self.field, x, y, self.d) for x, y in zip(a, b))
        return zero, blocks

    def fourier(self, values):
        output = []
        for s in range(self.prime):
            total = 0
            for a, value in enumerate(values):
                total ^= self.field.mul(value, self.field.power(self.zeta, a * s))
            output.append(total)
        return output

    def inverse_fourier_formula(self, moments):
        zero, blocks = self.fourier_formula(moments)
        minus = poly_power(self.field, self.u, self.d // 2, self.d)
        return zero, tuple(poly_mul(self.field, minus, block, self.d) for block in blocks)


class PrimePowerMoments:
    """Retain every valuation layer and its characteristic-two unit moments."""

    def __init__(self, prime, exponent):
        self.prime, self.exponent = prime, exponent
        self.modulus = prime ** exponent
        self.field, order, root, self.zeta = coefficient_field(self.modulus, additive_root=True)
        self.local = PrimePowerBlocks(prime, exponent, self.field, order, root)
        self.d, self.h = self.local.d, self.local.odd_order
        self.u = (1, 1) + (0,) * (self.d - 2)
        self.images = [self.local.image(a)[1:] for a in range(self.modulus)]
        _, fourier_delta = self.moments([self.field.power(self.zeta, a)
                                         for a in range(self.modulus)])
        self.gauss = tuple(tuple(block[x * length + length - 1] for x in range(self.d))
                           for block, length in zip(fourier_delta, self.local.lengths))

    def moments(self, values):
        blocks = [[0] * (self.d * length) for length in self.local.lengths]
        for a in range(1, self.modulus):
            for block, image in zip(blocks, self.images[a]):
                for k, value in enumerate(image):
                    block[k] ^= self.field.mul(values[a], value)
        return values[0], tuple(tuple(block) for block in blocks)

    def iota(self, values):
        transformed_x = list(poly_power(self.field, self.u, self.d - 1, self.d))
        transformed_x[0] ^= 1
        answer, power = [0] * self.d, (1,) + (0,) * (self.d - 1)
        for value in values:
            for k, coefficient in enumerate(power):
                answer[k] ^= self.field.mul(value, coefficient)
            power = poly_mul(self.field, power, transformed_x, self.d)
        return tuple(answer)

    def fourier_formula(self, moments):
        zero, blocks = moments
        output = []
        for character, length in enumerate(self.local.lengths):
            incoming = blocks[-character % self.h]
            value = [0] * (self.d * length)
            for w in range(length):
                v = length - 1 - w
                column = tuple(incoming[x * length + v] for x in range(self.d))
                column = poly_mul(self.field, self.gauss[character], self.iota(column), self.d)
                for x in range(self.d):
                    value[x * length + w] = column[x]
                if character == 0:
                    correction = zero
                    for higher in range(v + 1, self.exponent):
                        correction ^= blocks[0][higher]
                    value[(self.d - 1) * length + w] ^= correction
            output.append(tuple(value))
        new_zero = zero
        for v in range(self.exponent):
            new_zero ^= blocks[0][v]
        return new_zero, tuple(output)

    def fourier(self, values):
        output = []
        for s in range(self.modulus):
            total = 0
            for a, value in enumerate(values):
                total ^= self.field.mul(value, self.field.power(self.zeta, a * s))
            output.append(total)
        return output

    def inverse_fourier_formula(self, moments):
        zero, blocks = self.fourier_formula(moments)
        minus = poly_power(self.field, self.u, self.d // 2, self.d)
        output = []
        for block, length in zip(blocks, self.local.lengths):
            result = [0] * len(block)
            for v in range(length):
                column = poly_mul(self.field, minus,
                                  tuple(block[x * length + v] for x in range(self.d)), self.d)
                for x in range(self.d):
                    result[x * length + v] = column[x]
            output.append(tuple(result))
        return zero, tuple(output)

    def pointwise_formula(self, left, right):
        a0, a = left
        b0, b = right
        result = [[0] * (self.d * length) for length in self.local.lengths]
        for v in range(self.exponent):
            active = [chi for chi, length in enumerate(self.local.lengths) if v < length]
            for chi, psi in product(active, repeat=2):
                out = (chi + psi) % self.h
                for i, j in product(range(self.d), repeat=2):
                    result[out][(i & j) * self.local.lengths[out] + v] ^= self.field.mul(
                        a[chi][i * self.local.lengths[chi] + v],
                        b[psi][j * self.local.lengths[psi] + v])
        return self.field.mul(a0, b0), tuple(tuple(block) for block in result)

    def product_formula(self, left, right):
        a0, a = left
        b0, b = right
        blocks = tuple(block_mul(self.field, (self.d, length), x, y)
                       for length, x, y in zip(self.local.lengths, a, b))
        for v in range(self.exponent):
            a0 ^= a[0][v]
            b0 ^= b[0][v]
        zero = self.field.mul(a0, b0)
        for v in range(self.exponent):
            zero ^= blocks[0][v]
        return zero, blocks


def check_monoids(counts, models):
    rng = random.Random(3110202601)
    for modulus in (3, 5, 7, 9, 15, 17, 25, 27, 65, 257):
        model = MonoidBlocks(modulus)
        pairs = product(range(modulus), repeat=2) if modulus <= 27 else (
            (rng.randrange(modulus), rng.randrange(modulus)) for _ in range(80))
        for a, b in pairs:
            assert model.multiply(model.images[a], model.images[b]) == model.images[a * b % modulus]
            counts["basis_product_identities"] += 1
        if modulus <= 65:
            matrix = [list(sum((tuple(block) for block in image), ())) for image in model.images]
            assert field_rank(model.field, matrix) == modulus
            counts["full_change_of_basis_ranks"] += 1
        literal_count = sum(prod(2 * length - 1 for length in shape) for shape in model.shapes)
        formula_count = prod(local_product_count(p, e) for p, e in factor_powers(modulus))
        assert literal_count == formula_count
        assert formula_count <= modulus * 4 ** len(factor_powers(modulus))
        for _ in range(3 if modulus <= 27 else 1):
            left = [rng.randrange(model.field.size) for _ in range(modulus)]
            right = [rng.randrange(model.field.size) for _ in range(modulus)]
            direct = [0] * modulus
            for a, b in product(range(modulus), repeat=2):
                direct[a * b % modulus] ^= model.field.mul(left[a], right[b])
            assert model.transform(direct) == model.multiply(model.transform(left), model.transform(right))
            counts["full_monoid_convolutions"] += 1
        models.append({
            "modulus": modulus, "field_degree": model.field.degree,
            "binary_defining_polynomial": bin(model.field.polynomial),
            "block_shapes": [list(shape) for shape in model.shapes],
            "total_dimension": sum(prod(shape) for shape in model.shapes),
            "bilinear_product_count": formula_count,
        })


def check_grid_products(counts):
    field, rng = construct_field(6), random.Random(3110202602)
    for shape in ((), (2,), (4,), (8,), (2, 2), (2, 3), (4, 2), (2, 4), (2, 2, 2)):
        for _ in range(5):
            left = tuple(rng.randrange(field.size) for _ in block_indices(shape))
            right = tuple(rng.randrange(field.size) for _ in block_indices(shape))
            actual, rank = grid_product(field, shape, left, right)
            assert actual == block_mul(field, shape, left, right)
            assert rank == prod(2 * length - 1 for length in shape)
            counts["literal_multivariate_interpolations"] += 1
            counts["charged_bilinear_scalar_products"] += rank


def check_prime_closure(counts, models):
    rng = random.Random(3110202603)
    for prime in (3, 5, 7, 13, 17, 31):
        model = PrimeMoments(prime)
        field = model.field
        assert all(block[0] for block in model.gauss)
        for character in range(model.h):
            value = poly_mul(field, model.gauss[character],
                             model.iota(model.gauss[-character % model.h]), model.d)
            expected = list(poly_power(field, model.u, model.d // 2, model.d))
            if character == 0:
                expected[-1] ^= 1
            assert value == tuple(expected)
            counts["Gauss_unit_identities"] += 1
        for _ in range(18):
            left = [rng.randrange(field.size) for _ in range(prime)]
            right = [rng.randrange(field.size) for _ in range(prime)]
            ml, mr = model.moments(left), model.moments(right)
            assert model.fourier_formula(ml) == model.moments(model.fourier(left))
            assert model.inverse_fourier_formula(model.fourier_formula(ml)) == ml
            counts["complete_Fourier_Mellin_transforms"] += 1
            actual_pointwise = model.moments([field.mul(a, b) for a, b in zip(left, right)])
            assert model.pointwise_formula(ml, mr) == actual_pointwise
            counts["character_AND_product_closures"] += 1
            actual_product, actual_sum = [0] * prime, [0] * prime
            for a, b in product(range(prime), repeat=2):
                value = field.mul(left[a], right[b])
                actual_product[a * b % prime] ^= value
                actual_sum[(a + b) % prime] ^= value
            assert model.product_formula(ml, mr) == model.moments(actual_product)
            sum_moments = model.inverse_fourier_formula(model.pointwise_formula(
                model.fourier_formula(ml), model.fourier_formula(mr)))
            assert sum_moments == model.moments(actual_sum)
            counts["coupled_addition_and_multiplication_gates"] += 1
        # Ordinary scalar characters miss the entire 2-primary radical.
        unit = pow(model.local.generator, model.h, prime)
        radical = [0] * prime
        radical[1] ^= 1
        radical[unit] ^= 1
        _, radical_moments = model.moments(radical)
        assert all(block[0] == 0 for block in radical_moments)
        assert any(any(block) for block in radical_moments)
        counts["scalar_character_loss_controls"] += 1
        models.append({"prime": prime, "even_unit_order": model.d,
                       "odd_unit_order": model.h, "field_degree": field.degree})


def check_source_arithmetic(counts, models):
    for m in range(2, 26, 2):
        modulus = (1 << m) + 1
        factors = factor_powers(modulus)
        minimum_order = 1 << (valuation(m, 2, m.bit_length()) + 1)
        for prime, exponent in factors:
            order, power = 1, 2 % prime
            while power != 1:
                power = 2 * power % prime
                order += 1
            assert order & -order == minimum_order
            assert (prime - 1) % minimum_order == 0
            h = odd_part(prime - 1)
            d = (prime - 1) // h
            dimension = 1 + d * h * sum(
                totient_prime_power(prime, r) * (exponent - r) for r in range(exponent))
            assert dimension == prime ** exponent
            count = local_product_count(prime, exponent)
            closed = 1 + h * (2 * d - 1) * (
                2 * (prime ** exponent - 1) // (prime - 1) - prime ** (exponent - 1))
            assert count == closed
            assert count * prime <= 2 * (prime + 1) * prime ** exponent
            counts["source_prime_power_dimension_and_count_controls"] += 1
        rank = prod(local_product_count(p, e) for p, e in factors)
        assert rank <= modulus * 4 ** len(factors)
        assert (minimum_order + 1) ** len(factors) <= modulus
        counts["source_moduli"] += 1
        models.append({"m": m, "M": modulus, "factorization": factors,
                       "minimum_even_unit_order": minimum_order,
                       "bilinear_products_over_the_splitting_algebra": rank})


def check_prime_power_closure(counts, models):
    rng = random.Random(3110202604)
    for prime, exponent in ((3, 2), (5, 2), (3, 3), (7, 2)):
        model = PrimePowerMoments(prime, exponent)
        field, modulus = model.field, model.modulus
        assert all(block[0] for block in model.gauss)
        for character in range(model.h):
            actual = poly_mul(field, model.gauss[character],
                              model.iota(model.gauss[-character % model.h]), model.d)
            expected = list(poly_power(field, model.u, model.d // 2, model.d))
            if character == 0:
                expected[-1] ^= 1
            assert actual == tuple(expected)
            counts["prime_power_Gauss_unit_identities"] += 1
        for _ in range(7):
            left = [rng.randrange(field.size) for _ in range(modulus)]
            right = [rng.randrange(field.size) for _ in range(modulus)]
            ml, mr = model.moments(left), model.moments(right)
            assert model.fourier_formula(ml) == model.moments(model.fourier(left))
            assert model.inverse_fourier_formula(model.fourier_formula(ml)) == ml
            counts["prime_power_Fourier_Mellin_transforms"] += 1
            actual_pointwise = model.moments([field.mul(a, b) for a, b in zip(left, right)])
            assert model.pointwise_formula(ml, mr) == actual_pointwise
            counts["valuation_layer_AND_closures"] += 1
            actual_product, actual_sum = [0] * modulus, [0] * modulus
            for a, b in product(range(modulus), repeat=2):
                value = field.mul(left[a], right[b])
                actual_product[a * b % modulus] ^= value
                actual_sum[(a + b) % modulus] ^= value
            assert model.product_formula(ml, mr) == model.moments(actual_product)
            addition = model.inverse_fourier_formula(model.pointwise_formula(
                model.fourier_formula(ml), model.fourier_formula(mr)))
            assert addition == model.moments(actual_sum)
            counts["prime_power_coupled_arithmetic_gates"] += 1
        # The actual two-output copy tensor is an OR rule.  It is not
        # replaced by a pair of independent marginal distributions.
        labels = [(chi, x, v) for chi, length in enumerate(model.local.lengths)
                  for x in range(model.d) for v in range(length)]
        selected_pairs = product(labels, repeat=2) if modulus <= 9 else (
            (rng.choice(labels), rng.choice(labels)) for _ in range(90))
        for (chi, i, v), (psi, j, w) in selected_pairs:
            for a in range(modulus):
                left = model.images[a][chi][i * model.local.lengths[chi] + v]
                right = model.images[a][psi][j * model.local.lengths[psi] + w]
                character = (chi + psi) % model.h
                expected = (model.images[a][character][(i | j) * model.local.lengths[character] + v]
                            if v == w else 0)
                assert field.mul(left, right) == expected
                counts["joint_copy_OR_coefficients"] += 1
        models.append({"prime": prime, "exponent": exponent, "modulus": modulus,
                       "even_unit_order": model.d, "odd_unit_order": model.h,
                       "field_degree": field.degree,
                       "valuation_lengths": model.local.lengths})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts, monoids, closure, prime_power_closure, source = Counter(), [], [], [], []
    check_monoids(counts, monoids)
    check_grid_products(counts)
    check_prime_closure(counts, closure)
    check_prime_power_closure(counts, prime_power_closure)
    check_source_arithmetic(counts, source)
    evidence = {
        "status": "all exact checks passed",
        "scope": "constructive monoid blocks, every-prime-power Fourier--Mellin closure, and joint copy tensors; no factoring exponent improvement",
        "counts": dict(sorted(counts.items())),
        "monoid_models": monoids,
        "prime_closure_models": closure,
        "prime_power_closure_models": prime_power_closure,
        "source_arithmetic_models": source,
    }
    payload = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(json.dumps({"status": evidence["status"], "counts": evidence["counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
