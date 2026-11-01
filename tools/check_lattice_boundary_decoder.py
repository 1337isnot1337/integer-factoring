#!/usr/bin/env python3
"""Exact affine-lattice, nearest-plane, and grouped boundary decoding controls."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import random

from check_adaptive_transform_constructor import TransformModel, ceil_div
from check_operand_length_transform import LengthTransform
from check_observable_boundary_constructor import BoundaryNetwork, fa_message, lagrange_coefficients
from check_factoring_primitive_recipients import determinant, gram_schmidt, lll_reduce, nearest


def enhanced_model(h, width):
    a, d = 2 * h, 2 * h * h
    m = d * ceil_div(2 * width + (d.bit_length() - 1) + a, d)
    return TransformModel(h, width, m)


class EnhancedLengthTransform(LengthTransform):
    def __init__(self, h, n, ell):
        super().__init__(h, n, ell)
        width = self.digit_bits
        d, a = self.digits, self.a
        m = d * ceil_div(2 * width + self.log_digits + a, d)
        TransformModel.__init__(self, h, width, m)
        self.capacity = self.length * width


def orthogonal_vectors(basis):
    vectors = []
    for row in basis:
        current = list(map(Fraction, row))
        for previous in vectors:
            norm = sum(value * value for value in previous)
            multiple = sum(x * y for x, y in zip(row, previous)) / norm
            current = [x - multiple * y for x, y in zip(current, previous)]
        vectors.append(current)
    return vectors


class BoundaryLattice:
    def __init__(self, model, height, recipients, counts):
        self.model, self.height, self.recipients = model, height, tuple(recipients)
        self.nodes = tuple(model.powers[(model.a * j) % model.length] for j in recipients)
        self.count = min(height, len(recipients))
        self.counts = counts
        assert height < model.a
        if height <= len(recipients):
            self.basis = None
            return
        size, modulus = len(recipients), model.modulus
        high = []
        for column in range(size, height):
            values = [(-pow(node, column, modulus)) % modulus for node in self.nodes]
            high.append(lagrange_coefficients(self.nodes, values, modulus))
        basis = [[modulus * int(i == j) for j in range(height)] for i in range(size)]
        basis += [list(values) + [int(i == j) for j in range(height - size)]
                  for i, values in enumerate(high)]
        assert abs(determinant(basis)) == modulus ** size
        assert all(sum(value * pow(node, i, modulus) for i, value in enumerate(row)) % modulus == 0
                   for row in basis for node in self.nodes)
        self.basis, stats = lll_reduce(basis)
        assert abs(determinant(self.basis)) == modulus ** size
        self.orthogonal = orthogonal_vectors(self.basis)
        self.norms = [sum(value * value for value in row) for row in self.orthogonal]
        assert gram_schmidt(self.basis)[1] == self.norms
        assert min(self.norms) > model.base ** 2 * height
        assert all(sum(value * pow(node, i, modulus) for i, value in enumerate(row)) % modulus == 0
                   for row in self.basis for node in self.nodes)
        counts["exact_LLL_bases"] += 1
        counts["LLL_swaps"] += stats["swaps"]
        counts["LLL_size_reductions"] += stats["size_reductions"]

    def decode(self, residue, prescribed):
        model, height = self.model, self.height
        if len(prescribed) != len(self.recipients) or any(not 0 <= x < model.modulus for x in prescribed):
            return None
        values = [value * pow(model.powers[(residue * recipient) % model.length], -1,
                              model.modulus) % model.modulus
                  for recipient, value in zip(self.recipients, prescribed)]
        if height == 0:
            candidate = ()
        elif self.basis is None:
            candidate = lagrange_coefficients(self.nodes[:height], values[:height], model.modulus)
            self.counts["direct_Vandermonde_decodings"] += 1
        else:
            initial = lagrange_coefficients(self.nodes, values, model.modulus)
            offset = list(initial) + [0] * (height - self.count)
            residual = [Fraction(model.base - 1, 2) - x for x in offset]
            lattice = [0] * height
            for row, orthogonal, norm in reversed(list(zip(self.basis, self.orthogonal, self.norms))):
                coefficient = nearest(sum(x * y for x, y in zip(residual, orthogonal)) / norm)
                residual = [x - coefficient * y for x, y in zip(residual, row)]
                lattice = [x + coefficient * y for x, y in zip(lattice, row)]
            candidate = tuple(x + y for x, y in zip(offset, lattice))
            self.counts["nearest_plane_decodings"] += 1
        if any(not 0 <= value < model.base for value in candidate):
            return None
        if any(fa_message(model, residue, candidate, j) != value
               for j, value in zip(self.recipients, prescribed)):
            return None
        return candidate


def check_complete_input_maps(counts, records):
    rng = random.Random(713952)
    models = [enhanced_model(2, 1), enhanced_model(2, 2), enhanced_model(4, 1)]
    for model in models:
        a, b = model.a, model.base
        hidden_sets = [tuple(range(start, start + size))
                       for size in range(1, a // 2 + 1)
                       for start in range(a - size + 1)] if a == 4 else [
                           tuple(range(a // 2)), tuple(range(1, a // 2 + 1)),
                           tuple(range(a // 2, a)), (2, 3)]
        for hidden in hidden_sets:
            recipients = tuple(j for j in range(a) if j not in hidden)
            for height in range(a):
                decoder = BoundaryLattice(model, height, recipients, counts)
                residue = (height + hidden[0]) % a
                images = {}
                for digits in product(range(b), repeat=height):
                    values = tuple(fa_message(model, residue, digits, j) for j in recipients)
                    assert values not in images
                    images[values] = digits
                    assert decoder.decode(residue, values) == digits
                    counts["complete_bounded_input_recoveries"] += 1
                for _ in range(8):
                    values = tuple(rng.randrange(model.modulus) for _ in recipients)
                    assert decoder.decode(residue, values) == images.get(values)
                    counts["arbitrary_observation_comparisons"] += 1
                values = [0] * len(recipients)
                values[0] = model.modulus
                assert decoder.decode(residue, values) is None
                counts["invalid_observation_rejections"] += 1
                counts["complete_bounded_input_maps"] += 1
        records.append({"A": a, "w": model.digit_bits, "m": model.modulus_bits,
                        "M": str(model.modulus), "kind": "complete digit boxes"})


def check_large_digit_boxes(counts, records):
    rng = random.Random(194278)
    for model in [enhanced_model(2, 35), enhanced_model(4, 14), enhanced_model(8, 3)]:
        a = model.a
        for hidden in [tuple(range(a // 2)), tuple(range(a // 4, 3 * a // 4)),
                       tuple(range(a // 2, a))]:
            recipients = tuple(j for j in range(a) if j not in hidden)
            for height in sorted({a // 2 + 1, a - 1}):
                decoder = BoundaryLattice(model, height, recipients, counts)
                assert decoder.basis is not None
                inputs = [(0,) * height, (model.base - 1,) * height]
                inputs += [tuple((model.base - 1) * (i % 2) for i in range(height))]
                inputs += [tuple(rng.randrange(model.base) for _ in range(height)) for _ in range(15)]
                for index, digits in enumerate(inputs):
                    residue = index % a
                    values = tuple(fa_message(model, residue, digits, j) for j in recipients)
                    assert decoder.decode(residue, values) == digits
                    counts["large_digit_box_recoveries"] += 1
        records.append({"A": a, "w": model.digit_bits, "m": model.modulus_bits,
                        "M": str(model.modulus), "kind": "corners and deterministic samples"})


def lattice_entry(network, group, boundary, cache, counts, padding=False):
    if not network.boundary_valid(group, boundary, padding):
        return 0
    recipients = tuple(j for j in range(network.model.a) if j not in network.blocks[group])
    inputs = {}
    for family, residue in product(range(2), network.blocks[group]):
        height = network.heights[family, residue]
        key = (height, recipients)
        if key not in cache:
            cache[key] = BoundaryLattice(network.model, height, recipients, counts)
        prescribed = tuple(boundary["fa", family, residue, j] for j in recipients)
        digits = cache[key].decode(residue, prescribed)
        if digits is None or not network.input_valid(family, residue, digits):
            return 0
        inputs[family, residue] = digits
    return int(network.accepts(group, boundary, inputs))


def check_complete_entries(counts, records):
    rng = random.Random(9823)
    examples = [
        (enhanced_model(2, 1), 3, 5),
        (EnhancedLengthTransform(2, 9, 1), 1, 269),
        (EnhancedLengthTransform(2, 10, 10), 617, 1),
        (EnhancedLengthTransform(2, 19, 1), 1, 267013),
        (EnhancedLengthTransform(2, 19, 19), 267013, 1),
        (EnhancedLengthTransform(2, 20, 10), 719, 1451),
    ]
    for model, left, right in examples:
        cache = {}
        target = left * right
        for cut in sorted({0, target.bit_length() // 2, target.bit_length()}):
            network = BoundaryNetwork(model, 2, target, cut, left >> cut)
            trace = network.trace(left, right)
            for group in range(2):
                boundary = network.boundary(group, trace)
                answer = lattice_entry(network, group, boundary, cache, counts)
                reference, trials = network.entry(group, boundary, False)
                elementary, _ = network.entry(group, boundary, True)
                assert answer == reference == elementary == 1
                counts["complete_local_entry_comparisons"] += 1
                counts["reference_local_assignments"] += trials
                for _ in range(2):
                    changed = dict(boundary)
                    edge = rng.choice(tuple(boundary))
                    bound = model.carry_bound + 1 if edge[0] == "carry" else model.modulus
                    changed[edge] = (changed[edge] + rng.randrange(1, bound)) % bound
                    answer = lattice_entry(network, group, changed, cache, counts)
                    reference, trials = network.entry(group, changed, False)
                    assert answer == reference
                    counts["mutated_local_entry_comparisons"] += 1
                    counts["reference_local_assignments"] += trials
                assert lattice_entry(network, group, boundary, cache, counts, padding=True) == 0
                counts["invalid_padding_rejections"] += 1
        records.append({"N": target, "left": left, "right": right,
                        "w": model.digit_bits, "m": model.modulus_bits})


def check_parameter_bounds(counts, records):
    upper_f = Fraction(276407, 100000) + Fraction(1, 1024)
    upper_c, epsilon = Fraction(100205, 100000), Fraction(1, 512)
    exact_bound = 6 * (1 + epsilon) ** 2 * upper_f + 2 * epsilon * (1 + epsilon) * upper_f
    exact_bound += 4 * epsilon * (1 + epsilon) * upper_c
    assert exact_bound < 17
    for p in list(range(2, 130)) + [255, 256, 257, 511, 512, 1023, 1024, 4095, 4096]:
        h = 1 << ((p ** 9 - 1).bit_length())
        a, d = 2 * h, 2 * h * h
        t, k, n0 = d.bit_length() - 1, ceil_div(a, p), p ** 40
        denominator = a * a - a - 1
        preliminary0 = ceil_div(n0 + 1, denominator)
        candidates = {n0, n0 + 1, n0 + d - 1, n0 + d, 2 * n0, 8 * n0}
        candidates.update(denominator * (preliminary0 + shift) + delta
                          for shift in range(4) for delta in (-2, -1, 0, 1))
        for n in sorted(value for value in candidates if value >= n0):
            for even in (False, True):
                w = ceil_div(n + 1, denominator)
                if even:
                    w += w % 2
                m = d * ceil_div(2 * w + t + a, d)
                f, c = m + 1, w + t
                width = 6 * k * k * f + 2 * k * f + 4 * k * c
                assert 1024 * a ** 3 <= n
                assert m >= 2 * w + t + a
                assert m + 1 <= 2 * w + t + a + d
                assert p * p * width < 17 * n
                assert (1 << t) * 2 > a - 1
                counts["exact_parameter_cases"] += 1
                if n == n0 and p in (2, 3, 4, 8, 16, 128) and even:
                    records.append({"P": p, "n": str(n), "A": a, "D": d,
                                    "w": str(w), "m": str(m), "lambda": str(width),
                                    "B_tab": 1, "scaled_volume": [str(p * p * width), str(n)]})
    return [exact_bound.numerator, exact_bound.denominator]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts, boxes, entries, parameters = Counter(), [], [], []
    check_complete_input_maps(counts, boxes)
    check_large_digit_boxes(counts, boxes)
    check_complete_entries(counts, entries)
    bound = check_parameter_bounds(counts, parameters)
    evidence = {"status": "all exact checks passed", "counts": dict(sorted(counts.items())),
                "scope": "polynomial entry decoder and full written top constructor; factoring exponent unchanged",
                "input_models": boxes, "entry_models": entries, "parameter_models": parameters,
                "uniform_volume_upper_bound": bound}
    if args.output:
        args.output.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": evidence["status"], "counts": evidence["counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
