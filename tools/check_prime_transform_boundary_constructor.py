#!/usr/bin/env python3
"""Exact narrow-prime transforms and bounded closest-vector decoding controls.

The finite CVP control enumerates a rational sphere.  The proved asymptotic
constructor invokes the published deterministic single-exponential CVP
algorithm; this checker does not claim that runtime for its enumerator.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import gcd, isqrt
from pathlib import Path
import random

from check_adaptive_transform_constructor import TransformModel, ceil_div
from check_operand_length_transform import LengthTransform
from check_observable_boundary_constructor import BoundaryNetwork, fa_message, lagrange_coefficients
from check_factoring_primitive_recipients import determinant, gram_schmidt, lll_reduce
from check_lattice_boundary_decoder import orthogonal_vectors


def prime(number):
    if number < 2:
        return False
    if number % 2 == 0:
        return number == 2
    return all(number % divisor for divisor in range(3, isqrt(number) + 1, 2))


def install_prime(model):
    threshold = 1 << (2 * model.digit_bits + model.log_digits)
    candidate = threshold + 1
    trials = 0
    while not prime(candidate):
        candidate += model.length
        trials += 1
    model.prime_trials = trials + 1
    model.prime_threshold = threshold
    model.modulus = candidate
    model.modulus_bits = candidate.bit_length()
    nonresidue = next(value for value in range(2, candidate)
                     if pow(value, (candidate - 1) // 2, candidate) == candidate - 1)
    model.root = pow(nonresidue, (candidate - 1) // model.length, candidate)
    model.inverse_root = pow(model.root, -1, candidate)
    model.inverse_length = pow(model.length, -1, candidate)
    model.powers = [pow(model.root, j, candidate) for j in range(model.length)]
    model.inverse_powers = [pow(model.inverse_root, j, candidate) for j in range(model.length)]
    assert pow(model.root, model.length // 2, candidate) == candidate - 1
    assert pow(model.root, model.length, candidate) == 1
    assert candidate > model.digits * model.base ** 2
    model.exact_recovery = model.digits * (model.base - 1) ** 2 < candidate


class PrimeTransform(TransformModel):
    def __init__(self, h, width):
        super().__init__(h, width)
        install_prime(self)


class PrimeLengthTransform(LengthTransform):
    def __init__(self, h, n, ell):
        super().__init__(h, n, ell)
        install_prime(self)


class ExactSphereDecoder:
    def __init__(self, model, height, recipients, counts):
        self.model, self.height, self.recipients = model, height, tuple(recipients)
        self.nodes = tuple(model.powers[model.a * j % model.length] for j in recipients)
        self.counts = counts
        self.basis = None
        assert height <= model.a - 1
        if height <= len(recipients):
            return
        size, modulus = len(recipients), model.modulus
        basis = [[modulus * int(i == j) for j in range(height)] for i in range(size)]
        for index, column in enumerate(range(size, height)):
            values = [(-pow(node, column, modulus)) % modulus for node in self.nodes]
            low = lagrange_coefficients(self.nodes, values, modulus)
            basis.append(list(low) + [int(index == j) for j in range(height - size)])
        assert abs(determinant(basis)) == modulus ** size
        self.basis, _ = lll_reduce(basis)
        assert abs(determinant(self.basis)) == modulus ** size
        assert all(sum(x * pow(node, i, modulus) for i, x in enumerate(row)) % modulus == 0
                   for row in self.basis for node in self.nodes)
        self.mu, self.norms = gram_schmidt(self.basis)
        self.orthogonal = orthogonal_vectors(self.basis)
        self.radius_squared = Fraction(height * (model.base - 1) ** 2, 4)
        assert modulus > 2 * height * (model.base - 1) ** 2
        counts["exact_affine_lattice_bases"] += 1

    def sphere(self, target):
        size = self.height
        coordinates = [sum(x * y for x, y in zip(target, vector)) / norm
                       for vector, norm in zip(self.orthogonal, self.norms)]
        coefficients, output = [0] * size, []

        def visit(index, remaining):
            if index < 0:
                point = tuple(sum(coefficients[i] * self.basis[i][j] for i in range(size))
                              for j in range(size))
                assert sum((x - y) ** 2 for x, y in zip(point, target)) <= self.radius_squared
                output.append(point)
                return
            center = coordinates[index] - sum(coefficients[j] * self.mu[j][index]
                                              for j in range(index + 1, size))
            square = remaining / self.norms[index]
            radius = isqrt(ceil_div(square.numerator, square.denominator)) + 1
            floor_center = center.numerator // center.denominator
            for value in range(floor_center - radius, floor_center + radius + 2):
                consumed = (value - center) ** 2 * self.norms[index]
                self.counts["finite_sphere_branch_checks"] += 1
                if consumed <= remaining:
                    coefficients[index] = value
                    visit(index - 1, remaining - consumed)

        visit(size - 1, self.radius_squared)
        assert len(output) <= 1
        self.counts["finite_sphere_decodings"] += 1
        return output

    def decode(self, residue, prescribed):
        model, height = self.model, self.height
        if len(prescribed) != len(self.recipients) or any(not 0 <= x < model.modulus for x in prescribed):
            return None
        values = [value * pow(model.powers[residue * j % model.length], -1, model.modulus)
                  % model.modulus for j, value in zip(self.recipients, prescribed)]
        if height == 0:
            candidate = ()
        elif self.basis is None:
            candidate = lagrange_coefficients(self.nodes[:height], values[:height], model.modulus)
        else:
            low = lagrange_coefficients(self.nodes, values, model.modulus)
            offset = list(low) + [0] * (height - len(self.recipients))
            target = [Fraction(model.base - 1, 2) - x for x in offset]
            points = self.sphere(target)
            if not points:
                return None
            candidate = tuple(x + y for x, y in zip(offset, points[0]))
        if any(not 0 <= value < model.base for value in candidate):
            self.counts["postdecoded_box_rejections"] += 1
            return None
        if any(fa_message(model, residue, candidate, j) != value
               for j, value in zip(self.recipients, prescribed)):
            return None
        return candidate


def check_prime_transforms(counts, records):
    rng = random.Random(671209)
    for model in [PrimeTransform(2, 1), PrimeTransform(2, 2), PrimeTransform(4, 1),
                  PrimeTransform(8, 2), PrimeTransform(2, 10)]:
        assert prime(model.modulus)
        assert (model.modulus - 1) % model.length == 0
        assert all(gcd(model.powers[i] - model.powers[j], model.modulus) == 1
                   for i in range(model.length) for j in range(i))
        for _ in range(16):
            left, right = [rng.randrange(1 << model.capacity) for _ in range(2)]
            values = model.coefficients(left, right)
            assert values == model.ordinary_convolution(left, right)
            recovered, _, valid = model.carry(values)
            assert valid and recovered == left * right
            counts["prime_transform_products"] += 1
            counts["prime_transform_coefficients"] += model.length
        records.append({"A": model.a, "w": model.digit_bits, "X": str(model.prime_threshold),
                        "Q": str(model.modulus), "root": str(model.root),
                        "prime_candidates_tested": model.prime_trials,
                        "dyadic_interval": model.modulus < 2 * model.prime_threshold})
    # This example lies outside the analytic cutoff. It records why a
    # prime in the requested dyadic interval cannot be assumed for all scales.
    small = PrimeTransform(2, 1)
    assert small.modulus >= 2 * small.prime_threshold
    counts["below_cutoff_dyadic_prime_failures"] += 1
    for n in range(2, 8):
        for ell in range(1, n + 1):
            model = PrimeLengthTransform(2, n, ell)
            for left in range(1 << (ell - 1), 1 << ell):
                for right in range(1 << model.right_bits):
                    model.verify_pair(left, right, counts)


def check_input_maps(counts, records):
    rng = random.Random(193123)
    for model in [PrimeTransform(2, 1), PrimeTransform(2, 2), PrimeTransform(4, 1)]:
        a, b = model.a, model.base
        hidden_sets = [tuple(range(start, start + size))
                       for size in range(1, a // 2 + 1)
                       for start in range(a - size + 1)] if a == 4 else [
                           tuple(range(a // 2)), tuple(range(1, a // 2 + 1)),
                           tuple(range(a // 2, a)), (2, 3)]
        for hidden in hidden_sets:
            recipients = tuple(j for j in range(a) if j not in hidden)
            for height in range(a):
                decoder = ExactSphereDecoder(model, height, recipients, counts)
                residue = (height + hidden[0]) % a
                images = {}
                for digits in product(range(b), repeat=height):
                    values = tuple(fa_message(model, residue, digits, j) for j in recipients)
                    assert values not in images
                    images[values] = digits
                    assert decoder.decode(residue, values) == digits
                    counts["complete_bounded_input_recoveries"] += 1
                for _ in range(10):
                    values = tuple(rng.randrange(model.modulus) for _ in recipients)
                    assert decoder.decode(residue, values) == images.get(values)
                    counts["arbitrary_observation_comparisons"] += 1
                invalid = (model.modulus,) + (0,) * (len(recipients) - 1)
                assert decoder.decode(residue, invalid) is None
                counts["invalid_observation_rejections"] += 1
                counts["complete_bounded_input_maps"] += 1
        records.append({"A": a, "w": model.digit_bits, "Q": model.modulus})


def cvp_entry(network, group, boundary, cache, counts, padding=False):
    if not network.boundary_valid(group, boundary, padding):
        return 0
    recipients = tuple(j for j in range(network.model.a) if j not in network.blocks[group])
    inputs = {}
    for family, residue in product(range(2), network.blocks[group]):
        height = network.heights[family, residue]
        key = (height, recipients)
        if key not in cache:
            cache[key] = ExactSphereDecoder(network.model, height, recipients, counts)
        prescribed = tuple(boundary["fa", family, residue, j] for j in recipients)
        digits = cache[key].decode(residue, prescribed)
        if digits is None or not network.input_valid(family, residue, digits):
            return 0
        inputs[family, residue] = digits
    return int(network.accepts(group, boundary, inputs))


def check_local_entries(counts, records):
    rng = random.Random(20843)
    examples = [(PrimeTransform(2, 1), 3, 5),
                (PrimeLengthTransform(2, 9, 1), 1, 269),
                (PrimeLengthTransform(2, 10, 10), 617, 1),
                (PrimeLengthTransform(2, 19, 1), 1, 267013),
                (PrimeLengthTransform(2, 19, 19), 267013, 1),
                (PrimeLengthTransform(2, 20, 10), 719, 1451)]
    for model, left, right in examples:
        target, cache = left * right, {}
        for cut in sorted({0, target.bit_length() // 2, target.bit_length()}):
            network = BoundaryNetwork(model, 2, target, cut, left >> cut)
            trace = network.trace(left, right)
            for group in range(2):
                boundary = network.boundary(group, trace)
                actual = cvp_entry(network, group, boundary, cache, counts)
                reference, trials = network.entry(group, boundary, False)
                assert actual == reference == 1
                counts["complete_local_entry_comparisons"] += 1
                counts["reference_local_assignments"] += trials
                for _ in range(2):
                    changed = dict(boundary)
                    edge = rng.choice(tuple(changed))
                    bound = model.carry_bound + 1 if edge[0] == "carry" else model.modulus
                    changed[edge] = (changed[edge] + rng.randrange(1, bound)) % bound
                    reference, trials = network.entry(group, changed, False)
                    assert cvp_entry(network, group, changed, cache, counts) == reference
                    counts["mutated_local_entry_comparisons"] += 1
                    counts["reference_local_assignments"] += trials
                assert cvp_entry(network, group, boundary, cache, counts, padding=True) == 0
                counts["invalid_padding_rejections"] += 1
        records.append({"N": target, "left": left, "right": right,
                        "Q": model.modulus, "w": model.digit_bits})


def parameters(n):
    u = (n + 1).bit_length()
    p = isqrt(n) // (16 * u * u)
    a = 1 << ((p * u - 1).bit_length())
    d, denominator = a * a // 2, a * a - a - 1
    w = 2 * ceil_div(n + 1, 2 * denominator)
    t = d.bit_length() - 1
    f, c, k = 2 * w + t + 1, w + t, ceil_div(a, p)
    width = 6 * k * k * f + 2 * k * f + 4 * k * c
    return u, p, a, d, w, t, f, c, width


def check_large_parameters(counts, records):
    width_bound = (1 + Fraction(1, 2 ** 40)) / (1 - Fraction(1, 1024) - Fraction(1, 1024 ** 2))
    width_bound += Fraction(2, 64 * 41 ** 2)
    assert width_bound < Fraction(1001, 1000)
    volume_bound = 6 * Fraction(42, 41) ** 2 * Fraction(2003, 1000)
    volume_bound += 2 * Fraction(42, 41 ** 2) * Fraction(2003, 1000)
    volume_bound += 4 * Fraction(42, 41 ** 2) * Fraction(1002, 1000)
    assert volume_bound < 13
    ns = set()
    for exponent in list(range(40, 257)) + [511, 512, 1023, 1024, 2047, 2048, 4095, 4096]:
        anchor = 1 << exponent
        ns.update(anchor + offset for offset in (-3, -2, -1, 0, 1, 2, 3, 127))
        u = (anchor + 1).bit_length()
        initial_p = isqrt(anchor) // (16 * u * u)
        for p in range(initial_p, initial_p + 3):
            square_boundary = (16 * u * u * p) ** 2
            ns.update(square_boundary + offset for offset in (-2, -1, 0, 1, 2))
    for n in sorted(value for value in ns if value >= 2 ** 40):
        u, p, a, d, w, t, f, c, width = parameters(n)
        assert u >= 41 and p >= 2 and a >= 1024 and a >= 2 * p
        assert 64 * u * u * a * a < n
        assert w > 64 * u * u
        assert 2 * w + t > 128 * u * u
        assert a * a >= 10 ** 5
        assert 1000 * a * a * w < 1001 * n
        assert 1000 * a * a * f < 2003 * n
        assert 1000 * a * a * c < 1002 * n
        assert p * p * width < 13 * n
        assert a <= n and w <= n
        counts["exact_large_parameter_cases"] += 1
        if n in (1 << 40, 1 << 64, 1 << 128, 1 << 256, 1 << 1024, 1 << 4096):
            records.append({"n": str(n), "u": u, "P": str(p), "A": str(a),
                            "D": str(d), "w": w, "s": 2 * w + t, "f": f,
                            "lambda": str(width), "B_tab": 1,
                            "scaled_volume": [str(p * p * width), str(n)]})
    return {"digit_width": [width_bound.numerator, width_bound.denominator],
            "volume": [volume_bound.numerator, volume_bound.denominator]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts, primes, boxes, entries, schedules = Counter(), [], [], [], []
    check_prime_transforms(counts, primes)
    check_input_maps(counts, boxes)
    check_local_entries(counts, entries)
    rational_bounds = check_large_parameters(counts, schedules)
    evidence = {"status": "all exact checks passed", "counts": dict(sorted(counts.items())),
                "scope": "narrow-prime constructor and exact bounded CVP reduction; asymptotic CVP theorem is externally cited; factoring exponent unchanged",
                "prime_models": primes, "input_models": boxes, "entry_models": entries,
                "parameter_models": schedules, "rational_bounds": rational_bounds}
    if args.output:
        args.output.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": evidence["status"], "counts": evidence["counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
