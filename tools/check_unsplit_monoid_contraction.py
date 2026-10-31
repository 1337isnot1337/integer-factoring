#!/usr/bin/env python3
"""Exact descent, additive-subspace multiplication, and bundled MID controls.

No odd-character splitting field is used.  Every monoid basis change has
binary coefficients, including examples whose odd characters do not lie in
the coefficient field.  Fast multiplication includes evaluation, interpolation,
and polynomial reduction rather than counting pointwise products alone.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from math import prod
from pathlib import Path
import random

from check_carry_fourier_elimination import polynomial_remainder
from check_frobenius_sparse_recovery import construct_field
from check_multiplicative_monoid_blocks import (
    binomial_bit, block_indices, evaluate_polynomial, factor_powers,
    field_element_of_order, interpolation_weights, odd_part, primitive_unit,
    totient_prime_power, valuation,
)


def next_power_two(value):
    return 1 << (value - 1).bit_length()


def xor_all(values):
    answer = 0
    for value in values:
        answer ^= value
    return answer


def root_only_field(modulus):
    maximum_axis = 1
    for p, e in factor_powers(modulus):
        h, d = odd_part(p - 1), (p - 1) // odd_part(p - 1)
        maximum_axis = max(maximum_axis, d, e,
                           *(h * totient_prime_power(p, r) for r in range(e)))
    required_degree = (next_power_two(2 * maximum_axis - 1)).bit_length() - 1
    order = 1
    while pow(2, order, modulus) != 1:
        order += 1
    degree = order * max(1, (required_degree + order - 1) // order)
    field = construct_field(degree)
    return field, field_element_of_order(field, modulus)


def psi_polynomial(prime, h, r):
    if r == 0:
        return (1 << h) | 1
    step = h * prime ** (r - 1)
    return sum(1 << (j * step) for j in range(prime))


def binary_zeta(values):
    result = list(values)
    step = 1
    while step < len(result):
        for i in range(len(result)):
            if not i & step:
                result[i] ^= result[i | step]
        step <<= 1
    return result


def inverse_even_variable(values):
    translated = binary_zeta(values)
    return binary_zeta([translated[-i % len(values)] for i in range(len(values))])


class Axis:
    def __init__(self, length, prime=None, h=None, r=None):
        self.length, self.prime, self.h, self.r = length, prime, h, r
        self.polynomial = (psi_polynomial(prime, h, r) if prime is not None
                           else 1 << length)
        self.period = h * prime ** r if prime is not None else None
        self.remainders = [tuple(j for j in range(length)
                                 if polynomial_remainder(1 << i, self.polynomial) >> j & 1)
                           for i in range(2 * length - 1)]

    def reduce(self, values):
        if self.prime is None:
            return list(values[:self.length])
        folded = [0] * self.period
        for i, value in enumerate(values):
            folded[i % self.period] ^= value
        if self.r == 0:
            return folded
        step = self.period // self.prime
        for offset in range(step):
            value = folded[self.length + offset]
            for j in range(self.prime - 1):
                folded[j * step + offset] ^= value
        return folded[:self.length]


class UnsplitLocal:
    def __init__(self, prime, exponent):
        self.prime, self.exponent, self.modulus = prime, exponent, prime ** exponent
        self.h = odd_part(prime - 1)
        self.d = (prime - 1) // self.h
        self.odd_order = self.h * prime ** (exponent - 1)
        self.generator = primitive_unit(prime, exponent)
        self.logs = {pow(self.generator, k, self.modulus): k
                     for k in range(totient_prime_power(prime, exponent))}
        self.axes = [()] + [(
            Axis(self.h * totient_prime_power(prime, r), prime, self.h, r),
            Axis(self.d), Axis(exponent - r)) for r in range(exponent)]

    def unit_coordinates(self, residue):
        layer = valuation(residue, self.prime, self.exponent)
        unit = residue // self.prime ** layer
        log = self.logs[unit]
        even = log * pow(self.odd_order, -1, self.d) % self.d
        odd = (log * pow(self.d, -1, self.odd_order) % self.odd_order
               if self.odd_order > 1 else 0)
        return layer, even, odd

    def image(self, residue):
        result = [(1,)]
        for r, axes in enumerate(self.axes[1:]):
            c, d, length = (axis.length for axis in axes)
            if residue == 0:
                result.append((0,) * (c * d * length))
                continue
            layer, even, odd = self.unit_coordinates(residue)
            z_monomial = polynomial_remainder(1 << odd, axes[0].polynomial)
            result.append(tuple(int((z_monomial >> z & 1) and binomial_bit(even, x) and v == layer)
                                for z, x, v in block_indices((c, d, length))))
        return result

    def forward_flat(self, values):
        blocks = [[0] * prod(axis.length for axis in axes) for axes in self.axes]
        blocks[0][0] = xor_all(values)
        layers = [[[0] * self.d for _ in range(self.odd_order // self.prime ** v)]
                  for v in range(self.exponent)]
        for residue in range(1, self.modulus):
            v, even, odd = self.unit_coordinates(residue)
            layers[v][odd % len(layers[v])][even] ^= values[residue]
        for v, layer in enumerate(layers):
            moments = [binary_zeta(column) for column in layer]
            for r in range(self.exponent - v):
                axes = self.axes[r + 1]
                c, d, length = (axis.length for axis in axes)
                for x in range(d):
                    reduced = axes[0].reduce([column[x] for column in moments])
                    for z, value in enumerate(reduced):
                        blocks[r + 1][(z * d + x) * length + v] = value
        return sum((tuple(block) for block in blocks), ())

    def inverse_flat(self, values):
        blocks, offset = [], 0
        for axes in self.axes:
            size = prod(axis.length for axis in axes)
            blocks.append(values[offset:offset + size])
            offset += size
        layers = []
        for v in range(self.exponent):
            odd_length = self.odd_order // self.prime ** v
            moments = [[0] * self.d for _ in range(odd_length)]
            for r in range(self.exponent - v):
                axes = self.axes[r + 1]
                c, d, length = (axis.length for axis in axes)
                period = self.h * self.prime ** r
                for x in range(d):
                    coefficient = [blocks[r + 1][(z * d + x) * length + v] for z in range(c)]
                    if r > 0:
                        smaller = period // self.prime
                        folded = [0] * smaller
                        for z, value in enumerate(coefficient):
                            folded[z % smaller] ^= value
                    for z in range(odd_length):
                        local = z % period
                        if local < c:
                            moments[z][x] ^= coefficient[local]
                        if r > 0:
                            moments[z][x] ^= folded[z % smaller]
            layers.append([binary_zeta(column) for column in moments])
        result = [0] * self.modulus
        for residue in range(1, self.modulus):
            v, even, odd = self.unit_coordinates(residue)
            result[residue] = layers[v][odd % len(layers[v])][even]
        result[0] = blocks[0][0] ^ xor_all(result[1:])
        return result


def multiply_block(field, axes, left, right):
    shape = tuple(axis.length for axis in axes)
    indices = block_indices(shape)
    positions = {index: i for i, index in enumerate(indices)}
    result = [0] * len(indices)
    for i, a in enumerate(left):
        if not a:
            continue
        for j, b in enumerate(right):
            if not b:
                continue
            choices = [axis.remainders[x + y]
                       for axis, x, y in zip(axes, indices[i], indices[j])]
            if any(not choice for choice in choices):
                continue
            value = field.mul(a, b)
            for index in product(*choices):
                result[positions[index]] ^= value
    return tuple(result)


class UnsplitMonoid:
    def __init__(self, modulus):
        self.modulus = modulus
        self.field, self.root = root_only_field(modulus)
        self.local = [UnsplitLocal(p, e) for p, e in factor_powers(modulus)]
        choices = tuple(product(*(range(len(part.axes)) for part in self.local)))
        self.axes = [sum((part.axes[index] for part, index in zip(self.local, choice)), ())
                     for choice in choices]
        self.shapes = [tuple(axis.length for axis in axes) for axes in self.axes]
        self.images = []
        for residue in range(modulus):
            local_images = [part.image(residue % part.modulus) for part in self.local]
            image = []
            for choice in choices:
                image.append(tuple(prod(values) for values in product(*(
                    local[index] for local, index in zip(local_images, choice)))))
            self.images.append(image)
        self.crt_shape = tuple(part.modulus for part in self.local)
        tensor_positions = {index: i for i, index in enumerate(block_indices(self.crt_shape))}
        self.crt_positions = [tensor_positions[tuple(a % part.modulus for part in self.local)]
                              for a in range(modulus)]
        local_offsets = []
        for part in self.local:
            offsets, offset = [], 0
            for axes in part.axes:
                offsets.append(offset)
                offset += prod(axis.length for axis in axes)
            local_offsets.append(offsets)
        self.tensor_block_positions = []
        for choice in choices:
            ranges = [range(offsets[index], offsets[index] + prod(axis.length for axis in part.axes[index]))
                      for part, offsets, index in zip(self.local, local_offsets, choice)]
            self.tensor_block_positions.append(tuple(tensor_positions[index] for index in product(*ranges)))
        assert sum(prod(shape) for shape in self.shapes) == modulus
        self.offsets, offset = [], 0
        for shape in self.shapes:
            self.offsets.append(offset)
            offset += prod(shape)
        rows = [sum((sum(self.images[a], ())[i] & 1) << a for a in range(modulus))
                | (1 << (modulus + i)) for i in range(modulus)]
        for column in range(modulus):
            pivot = next(i for i in range(column, modulus) if rows[i] >> column & 1)
            rows[column], rows[pivot] = rows[pivot], rows[column]
            for i in range(modulus):
                if i != column and rows[i] >> column & 1:
                    rows[i] ^= rows[column]
        assert all(row & ((1 << modulus) - 1) == 1 << i for i, row in enumerate(rows))
        self.inverse_rows = [row >> modulus for row in rows]

    def transform(self, values):
        answer = [[0] * prod(shape) for shape in self.shapes]
        for value, image in zip(values, self.images):
            if value:
                for out, block in zip(answer, image):
                    for i, bit in enumerate(block):
                        if bit:
                            out[i] ^= value
        return [tuple(block) for block in answer]

    def inverse_transform(self, blocks):
        flat = sum((tuple(block) for block in blocks), ())
        return [xor_all(value for i, value in enumerate(flat) if row >> i & 1)
                for row in self.inverse_rows]

    def multiply(self, left, right):
        return [multiply_block(self.field, axes, a, b)
                for axes, a, b in zip(self.axes, left, right)]

    def transform_fast(self, values):
        tensor = [0] * self.modulus
        for position, value in zip(self.crt_positions, values):
            tensor[position] = value
        for axis, part in enumerate(self.local):
            tensor, _ = axis_map(tensor, self.crt_shape, axis, part.forward_flat)
        return [tuple(tensor[position] for position in positions)
                for positions in self.tensor_block_positions]

    def inverse_transform_fast(self, blocks):
        tensor = [0] * self.modulus
        for block, positions in zip(blocks, self.tensor_block_positions):
            for value, position in zip(block, positions):
                tensor[position] = value
        for axis, part in enumerate(self.local):
            tensor, _ = axis_map(tensor, self.crt_shape, axis, part.inverse_flat)
        return [tensor[position] for position in self.crt_positions]


class AdditiveFFT:
    """Sparse subspace-polynomial evaluation and interpolation, including cosets."""

    def __init__(self, field, depth):
        assert depth <= field.degree
        self.field, self.depth = field, depth
        self.subspace = [{1: 1}]
        self.constants, self.inverse_constants = [], []
        for k in range(depth):
            previous = self.subspace[-1]
            c = xor_all(field.mul(value, field.power(1 << k, power))
                        for power, value in previous.items())
            assert c
            self.constants.append(c)
            self.inverse_constants.append(field.inverse(c))
            following = {}
            for power, value in previous.items():
                following[2 * power] = following.get(2 * power, 0) ^ field.mul(value, value)
                following[power] = following.get(power, 0) ^ field.mul(c, value)
            self.subspace.append({power: value for power, value in following.items() if value})

    def shifted_divisor(self, depth, coset):
        result = dict(self.subspace[depth])
        constant = xor_all(self.field.mul(value, self.field.power(coset, power))
                           for power, value in result.items())
        if constant:
            result[0] = constant
        return result

    def evaluate(self, coefficients, coset=0):
        length = len(coefficients)
        assert length > 0 and length & (length - 1) == 0
        depth = length.bit_length() - 1
        if depth == 0:
            return list(coefficients)
        half, remainder = length // 2, list(coefficients)
        divisor = self.shifted_divisor(depth - 1, coset)
        quotient = [0] * half
        for i in range(length - 1, half - 1, -1):
            value = remainder[i]
            quotient[i - half] = value
            if value:
                for power, coefficient in divisor.items():
                    remainder[i - half + power] ^= self.field.mul(value, coefficient)
        assert not any(remainder[half:])
        first = remainder[:half]
        second = [a ^ self.field.mul(self.constants[depth - 1], b)
                  for a, b in zip(first, quotient)]
        return self.evaluate(first, coset) + self.evaluate(second, coset ^ half)

    def interpolate(self, values, coset=0):
        length = len(values)
        assert length > 0 and length & (length - 1) == 0
        depth = length.bit_length() - 1
        if depth == 0:
            return list(values)
        half = length // 2
        first = self.interpolate(values[:half], coset)
        second = self.interpolate(values[half:], coset ^ half)
        quotient = [self.field.mul(a ^ b, self.inverse_constants[depth - 1])
                    for a, b in zip(first, second)]
        result = first + [0] * half
        for i, value in enumerate(quotient):
            for power, coefficient in self.shifted_divisor(depth - 1, coset).items():
                result[i + power] ^= self.field.mul(value, coefficient)
        return result


def axis_map(values, shape, axis, operation, output_length=None):
    indices = block_indices(shape)
    input_positions = {index: i for i, index in enumerate(indices)}
    output_shape = list(shape)
    if output_length is not None:
        output_shape[axis] = output_length
    output_shape = tuple(output_shape)
    output_positions = {index: i for i, index in enumerate(block_indices(output_shape))}
    answer = [0] * prod(output_shape)
    other_shape = shape[:axis] + shape[axis + 1:]
    for other in block_indices(other_shape):
        line = [values[input_positions[other[:axis] + (j,) + other[axis:]]]
                for j in range(shape[axis])]
        image = operation(line)
        assert len(image) == output_shape[axis]
        for j, value in enumerate(image):
            answer[output_positions[other[:axis] + (j,) + other[axis:]]] = value
    return tuple(answer), output_shape


def reduce_tensor(values, shape, axes):
    for i, axis in enumerate(axes):
        values, shape = axis_map(values, shape, i, axis.reduce, axis.length)
    return tuple(values)


def fast_block_product(field, axes, left, right):
    shape = tuple(axis.length for axis in axes)
    grid = tuple(next_power_two(2 * length - 1) for length in shape)
    indices = block_indices(shape)
    positions = {index: i for i, index in enumerate(block_indices(grid))}
    fft = AdditiveFFT(field, max((length.bit_length() - 1 for length in grid), default=0))
    def forward(values):
        result = [0] * prod(grid)
        for index, value in zip(indices, values):
            result[positions[index]] = value
        for axis in range(len(grid)):
            result, _ = axis_map(result, grid, axis, fft.evaluate)
        return result
    a, b = forward(left), forward(right)
    result = [field.mul(x, y) for x, y in zip(a, b)]
    for axis in range(len(grid)):
        result, _ = axis_map(result, grid, axis, fft.interpolate)
    return reduce_tensor(result, grid, axes), prod(grid)


def inverse_z_x(field, axes, values):
    z_axis, x_axis = axes
    c, d = z_axis.length, x_axis.length
    result = [0] * (c * d)
    for x in range(d):
        reversed_polynomial = [0] * z_axis.period
        for z in range(c):
            reversed_polynomial[-z % z_axis.period] ^= values[z * d + x]
        reduced = z_axis.reduce(reversed_polynomial)
        for z, value in enumerate(reduced):
            result[z * d + x] = value
    for z in range(c):
        result[z * d:(z + 1) * d] = inverse_even_variable(result[z * d:(z + 1) * d])
    return tuple(result)


def descended_fourier(model, values):
    assert len(model.local) == 1
    part = model.local[0]
    incoming = model.transform_fast(values)
    exponential = model.transform_fast([model.field.power(model.root, a) for a in range(model.modulus)])
    zero = values[0]
    outgoing = [(values[0],)]
    for r, axes in enumerate(model.axes[1:]):
        c, d, length = (axis.length for axis in axes)
        gauss = tuple(exponential[r + 1][(z * d + x) * length + length - 1]
                      for z, x in block_indices((c, d)))
        result = [0] * (c * d * length)
        for w in range(length):
            v = length - 1 - w
            column = tuple(incoming[r + 1][(z * d + x) * length + v]
                           for z, x in block_indices((c, d)))
            transformed = inverse_z_x(model.field, axes[:2], column)
            image, _ = fast_block_product(model.field, axes[:2], gauss, transformed)
            for i, value in enumerate(image):
                result[i * length + w] = value
            if r == 0:
                correction = values[0]
                for higher in range(v + 1, part.exponent):
                    correction ^= xor_all(incoming[1][(z * d) * length + higher]
                                          for z in range(c))
                for z in range(c):
                    result[(z * d + d - 1) * length + w] ^= correction
        outgoing.append(tuple(result))
    # The first block is augmentation, equal to input f(0) after Fourier.
    # The actual output zero coefficient is checked separately below.
    return outgoing


def explicit_cp(model):
    """Literal optimal-size tensor grids, with every output factor descended."""
    field, factors = model.field, []
    for block_index, (axes, shape) in enumerate(zip(model.axes, model.shapes)):
        grid = tuple(2 * length - 1 for length in shape)
        weights = [interpolation_weights(field, length) for length in grid]
        for point in block_indices(grid):
            input_factor = []
            for image in model.images:
                input_factor.append(xor_all(field.product(value, *(
                    field.power(node, degree) for node, degree in zip(point, index)))
                    for value, index in zip(image[block_index], block_indices(shape))))
            output_polynomial = tuple(field.product(*(
                weight[node][degree] for weight, node, degree in zip(weights, point, index)))
                for index in block_indices(grid))
            output_block = reduce_tensor(output_polynomial, grid, axes)
            blocks = [(0,) * prod(other) for other in model.shapes]
            blocks[block_index] = output_block
            factors.append((tuple(input_factor), tuple(model.inverse_transform(blocks))))
    return factors


def check_basis_and_products(counts, records):
    rng = random.Random(3110202605)
    for modulus in (3, 5, 7, 9, 15, 23, 25, 27, 29, 65, 257):
        model = UnsplitMonoid(modulus)
        pairs = product(range(modulus), repeat=2) if modulus <= 29 else (
            (rng.randrange(modulus), rng.randrange(modulus)) for _ in range(60))
        for a, b in pairs:
            assert model.multiply(model.images[a], model.images[b]) == model.images[a * b % modulus]
            counts["binary_basis_products"] += 1
        rank_count = sum(prod(2 * length - 1 for length in shape) for shape in model.shapes)
        grid_count = sum(prod(next_power_two(2 * length - 1) for length in shape)
                         for shape in model.shapes)
        assert rank_count <= modulus * 8 ** len(model.local)
        assert grid_count <= modulus * 64 ** len(model.local)
        for _ in range(3):
            values = [rng.randrange(model.field.size) for _ in range(modulus)]
            assert model.inverse_transform(model.transform(values)) == values
            assert model.transform_fast(values) == model.transform(values)
            assert model.inverse_transform_fast(model.transform_fast(values)) == values
            counts["full_binary_change_of_basis_inversions"] += 1
        for _ in range(2 if modulus <= 29 else 1):
            left = [rng.randrange(model.field.size) for _ in range(modulus)]
            right = [rng.randrange(model.field.size) for _ in range(modulus)]
            direct = [0] * modulus
            for a, b in product(range(modulus), repeat=2):
                direct[a * b % modulus] ^= model.field.mul(left[a], right[b])
            tl, tr = model.transform_fast(left), model.transform_fast(right)
            fast = []
            for axes, a, b in zip(model.axes, tl, tr):
                result, products = fast_block_product(model.field, axes, a, b)
                assert result == multiply_block(model.field, axes, a, b)
                fast.append(result)
                counts["fully_charged_subspace_products"] += products
            assert model.inverse_transform_fast(fast) == direct
            counts["complete_fast_monoid_convolutions"] += 1
        nonsplit = any((model.field.size - 1) % part.odd_order for part in model.local)
        counts["genuinely_unsplit_models"] += int(nonsplit)
        records.append({"modulus": modulus, "field_degree": model.field.degree,
                        "binary_defining_polynomial": bin(model.field.polynomial),
                        "odd_characters_not_all_in_the_field": nonsplit,
                        "block_shapes": [list(shape) for shape in model.shapes],
                        "optimal_grid_bilinear_products": rank_count,
                        "power_of_two_grid_products": grid_count})


def check_subspace_fft(counts):
    rng, field = random.Random(3110202606), construct_field(7)
    fft = AdditiveFFT(field, 6)
    for depth in range(7):
        length = 1 << depth
        for _ in range(5):
            polynomial = [rng.randrange(field.size) for _ in range(length)]
            coset = rng.randrange(field.size)
            values = fft.evaluate(polynomial, coset)
            assert values == [evaluate_polynomial(field, polynomial, coset ^ x) for x in range(length)]
            assert fft.interpolate(values, coset) == polynomial
            counts["affine_subspace_evaluation_and_interpolation"] += 1


def check_fourier_descent(counts):
    rng = random.Random(3110202607)
    for modulus in (3, 5, 7, 9, 23, 25, 27, 29):
        model = UnsplitMonoid(modulus)
        for _ in range(5):
            values = [rng.randrange(model.field.size) for _ in range(modulus)]
            direct = [xor_all(model.field.mul(value, model.field.power(model.root, a * s))
                              for a, value in enumerate(values)) for s in range(modulus)]
            transformed = descended_fourier(model, values)
            assert transformed == model.transform(direct)
            assert model.inverse_transform(transformed) == direct
            assert direct[0] == xor_all(values)
            counts["unsplit_complete_Fourier_transforms"] += 1


def check_literal_cp_and_mid(counts, records):
    models = {}
    for modulus in (3, 5, 7, 9):
        model = UnsplitMonoid(modulus)
        factors = explicit_cp(model)
        assert len(factors) == sum(prod(2 * length - 1 for length in shape) for shape in model.shapes)
        for a, b, c in product(range(modulus), repeat=3):
            actual = xor_all(model.field.product(left[a], left[b], right[c])
                             for left, right in factors)
            assert actual == int(a * b % modulus == c)
            counts["literal_scalar_CP_coefficients"] += 1
        models[modulus] = model, factors
    modulus, a_count, omega = 5, 2, 2
    model, factors = models[modulus]
    vectors = tuple(product(range(modulus), repeat=a_count))
    omega_inverse = pow(omega, -1, modulus)
    forward = {a: tuple(sum(a[r] * pow(omega, a_count * r * k, modulus)
                            for r in range(a_count)) % modulus for k in range(a_count)) for a in vectors}
    assert len(set(forward.values())) == len(vectors)
    for k1 in range(a_count):
        output_map = {z: tuple(pow(omega_inverse, k1 * j, modulus) * sum(
            z[k] * pow(omega_inverse, a_count * k * j, modulus) for k in range(a_count)) % modulus
                              for j in range(a_count)) for z in vectors}
        assert len(set(output_map.values())) == len(vectors)
        output_inverse = {c: z for z, c in output_map.items()}
        cp_labels = tuple(product(range(len(factors)), repeat=a_count))
        left_tables = [[model.field.product(*(factors[label[k]][0][forward[a][k]]
                                             for k in range(a_count))) for a in vectors]
                       for label in cp_labels]
        right_tables = [[model.field.product(*(factors[label[k]][1][output_inverse[c][k]]
                                              for k in range(a_count))) for c in vectors]
                        for label in cp_labels]
        for ai, bi, ci in product(range(len(vectors)), repeat=3):
            actual = xor_all(model.field.product(left[ai], left[bi], right[ci])
                             for left, right in zip(left_tables, right_tables))
            products = tuple(forward[vectors[ai]][k] * forward[vectors[bi]][k] % modulus
                             for k in range(a_count))
            assert actual == int(output_map[products] == vectors[ci])
            counts["complete_MID_bundle_CP_coefficients"] += 1
        records.append({"M": modulus, "A": a_count, "k1": k1,
                        "bundle_dimension": modulus ** a_count,
                        "CP_terms": len(cp_labels),
                        "input_setting_count": modulus ** (2 * a_count)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts, models, mid = Counter(), [], []
    check_subspace_fft(counts)
    check_basis_and_products(counts, models)
    check_fourier_descent(counts)
    check_literal_cp_and_mid(counts, mid)
    evidence = {
        "status": "all exact checks passed",
        "scope": "original-field monoid blocks, complete subspace multiplication, descended Fourier, and three-bundle MID decomposition; no stronger factoring exponent",
        "counts": dict(sorted(counts.items())), "models": models, "MID_models": mid,
    }
    if args.output:
        args.output.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": evidence["status"], "counts": evidence["counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
