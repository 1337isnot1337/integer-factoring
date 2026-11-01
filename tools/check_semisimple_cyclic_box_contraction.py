#!/usr/bin/env python3
"""Exact zero-safe box controls in tensor/cyclic recovery algebras."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import random

from check_carry_fourier_elimination import BinaryField
from check_frobenius_rational_trace import field_of_degree
from check_frobenius_sparse_recovery import TensorAlgebra
from frobenius_trace_infinity_correction import zero_safe_binomial_trace


COUNTS = Counter()


class CyclicTensor:
    def __init__(self, d, t, r):
        field = lambda degree: BinaryField(0b11) if degree == 1 else field_of_degree(degree)
        self.base = TensorAlgebra(field(d), field(t))
        self.d, self.t, self.r = d, t, r
        self.period = d * t * (r - 1)
        self.zero = (self.base.zero,) * r
        self.one = (self.base.one,) + (self.base.zero,) * (r - 1)

    def add(self, left, right):
        return tuple(self.base.add(a, b) for a, b in zip(left, right))

    def mul(self, left, right):
        result = [self.base.zero] * self.r
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                index = (i + j) % self.r
                result[index] = self.base.add(result[index], self.base.mul(a, b))
        return tuple(result)

    def power(self, value, exponent):
        result = self.one
        while exponent:
            if exponent & 1:
                result = self.mul(result, value)
            value = self.mul(value, value)
            exponent >>= 1
        return result

    def scalar(self, value):
        return (self.base.weight_scalar(value),) + (self.base.zero,) * (self.r - 1)

    def random(self, rng):
        return tuple(tuple(rng.randrange(1 << self.d) for _ in range(self.t))
                     for _ in range(self.r))

    def geometric(self, value, width):
        total, power = self.zero, self.one
        for _ in range(1 << width):
            total = self.add(total, power)
            power = self.mul(power, value)
        return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(0x5E1151)
    models = []
    for d, t, r in [(1, 1, 3), (2, 1, 3), (2, 2, 3), (2, 3, 3), (2, 2, 5), (3, 2, 3)]:
        ring = CyclicTensor(d, t, r)
        inverse_exponent = (1 << (ring.period + 1)) - 3
        cyclic_idempotent = (ring.base.one,) * r
        generator = (ring.base.zero, ring.base.one) + (ring.base.zero,) * (r - 2)
        annihilator = ring.add(ring.one, generator)
        assert cyclic_idempotent != ring.zero and annihilator != ring.zero
        assert ring.mul(cyclic_idempotent, annihilator) == ring.zero
        assert ring.mul(cyclic_idempotent, cyclic_idempotent) == cyclic_idempotent
        COUNTS["genuine_zero_divisor_controls"] += 1
        values = [ring.zero, ring.one, cyclic_idempotent, ring.add(ring.one, cyclic_idempotent)]
        values += [ring.random(rng) for _ in range(5)]
        for value in values:
            assert ring.power(value, 1 << ring.period) == value
            pseudo = ring.power(value, inverse_exponent)
            assert ring.mul(value, ring.mul(pseudo, value)) == value
            assert ring.mul(pseudo, ring.mul(value, pseudo)) == pseudo
            COUNTS["exact_Frobenius_and_reflexive_inverse_controls"] += 1
            for width in [1, 2, 5]:
                numerator = ring.add(ring.one, ring.power(value, 1 << width))
                denominator = ring.add(ring.one, value)
                actual = ring.mul(numerator, ring.power(denominator, inverse_exponent))
                assert actual == ring.geometric(value, width)
                COUNTS["exact_all_point_geometric_controls"] += 1
        for _ in range(3):
            scalar = ring.random(rng)
            numerators = [(ring.random(rng), ring.random(rng), rng.randrange(t))]
            denominators = [(ring.random(rng), ring.random(rng), rng.randrange(t))]
            actual = zero_safe_binomial_trace(t, ring.period, scalar, numerators, denominators,
                                              multiply=ring.mul, add=ring.add,
                                              zero=ring.zero, one=ring.one)
            full, punctured = ring.zero, ring.zero
            for field_value in range(1 << t):
                point = ring.scalar(field_value)
                value = scalar
                for a, b, exponent in numerators:
                    value = ring.mul(value, ring.add(a, ring.mul(b, ring.power(point, 1 << exponent))))
                for a, b, exponent in denominators:
                    denominator = ring.add(a, ring.mul(b, ring.power(point, 1 << exponent)))
                    value = ring.mul(value, ring.power(denominator, inverse_exponent))
                full = ring.add(full, value)
                if field_value:
                    punctured = ring.add(punctured, value)
                COUNTS["reference_scalar_field_evaluations"] += 1
            assert actual == (full, punctured)
            COUNTS["exact_zero_safe_algebra_trace_controls"] += 1
        models.append({"root_degree": d, "scalar_degree": t, "cyclic_prime": r,
                       "binary_dimension": d * t * r, "Frobenius_period": ring.period})
    bad = CyclicTensor(1, 1, 2)
    generator = (bad.base.zero, bad.base.one)
    epsilon = bad.add(bad.one, generator)
    assert epsilon != bad.zero and bad.mul(epsilon, epsilon) == bad.zero
    assert bad.geometric(generator, 1) == epsilon
    assert bad.add(bad.one, bad.power(generator, 2)) == bad.zero
    assert bad.power(generator, 1 << bad.period) != generator
    COUNTS["even_cyclic_nilpotent_countercontrols"] += 1
    result = {"status": "passed", "claim": "zero-safe box contractions in odd-prime cyclic recovery algebras",
              "counts": dict(sorted(COUNTS.items())), "models": models,
              "scope": "The Frobenius-period hypothesis is essential; arbitrary nilpotent algebras are excluded."}
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
