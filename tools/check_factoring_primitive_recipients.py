#!/usr/bin/env python3
"""Exact lattice preimages for primitive Fourier recipient sets.

Source fan-out: Rudich--Rousseau c91cc3f. Lattice reduction: classical
LLL with delta=3/4, implemented with exact rational arithmetic.
The finite checks do not establish a new factoring exponent.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import gcd
from pathlib import Path
import random

from check_factoring_network_structure import Fanout


def nearest(value: Fraction) -> int:
    """The exact floor(x+1/2) rule, including negative inputs and ties."""
    return (2 * value.numerator + value.denominator) // (2 * value.denominator)


def gram_schmidt(basis: list[list[int]]) -> tuple[list[list[Fraction]], list[Fraction]]:
    size = len(basis)
    orthogonal: list[list[Fraction]] = []
    norms: list[Fraction] = []
    mu = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for i in range(size):
        vector = [Fraction(x) for x in basis[i]]
        for j in range(i):
            mu[i][j] = sum(basis[i][a] * orthogonal[j][a] for a in range(size)) / norms[j]
            vector = [x - mu[i][j] * y for x, y in zip(vector, orthogonal[j])]
        orthogonal.append(vector)
        norms.append(sum(x * x for x in vector))
        if not norms[-1]:
            raise ValueError("singular lattice basis")
    return mu, norms


def lll_reduce(basis: list[list[int]]) -> tuple[list[list[int]], dict[str, int]]:
    """Integral row operations with exact incremental Gram--Schmidt updates."""
    basis = [row[:] for row in basis]
    size = len(basis)
    mu, norms = gram_schmidt(basis)
    index, swaps, reductions = 1, 0, 0
    while index < size:
        for previous in range(index - 1, -1, -1):
            multiple = nearest(mu[index][previous])
            if multiple:
                basis[index] = [x - multiple * y
                                for x, y in zip(basis[index], basis[previous])]
                for j in range(previous):
                    mu[index][j] -= multiple * mu[previous][j]
                mu[index][previous] -= multiple
                reductions += 1
        old_mu = mu[index][index - 1]
        if norms[index] >= (Fraction(3, 4) - old_mu * old_mu) * norms[index - 1]:
            index += 1
            continue
        basis[index], basis[index - 1] = basis[index - 1], basis[index]
        new_norm = norms[index] + old_mu * old_mu * norms[index - 1]
        new_mu = old_mu * norms[index - 1] / new_norm
        norms[index] = norms[index - 1] * norms[index] / new_norm
        norms[index - 1] = new_norm
        mu[index][index - 1] = new_mu
        for j in range(index - 1):
            mu[index][j], mu[index - 1][j] = mu[index - 1][j], mu[index][j]
        for i in range(index + 1, size):
            old = mu[i][index]
            mu[i][index] = mu[i][index - 1] - old_mu * old
            mu[i][index - 1] = old + new_mu * mu[i][index]
        index = max(1, index - 1)
        swaps += 1

    # Fresh orthogonalization independently checks all incremental updates.
    fresh_mu, fresh_norms = gram_schmidt(basis)
    assert (fresh_mu, fresh_norms) == (mu, norms)
    assert all(abs(mu[i][j]) <= Fraction(1, 2) for i in range(size) for j in range(i))
    assert all(norms[i] >= (Fraction(3, 4) - mu[i][i - 1] ** 2) * norms[i - 1]
               for i in range(1, size))
    return basis, {"swaps": swaps, "size_reductions": reductions}


def inverse_mod(matrix: list[list[int]], modulus: int) -> list[list[int]]:
    size = len(matrix)
    rows = [[x % modulus for x in row] + [int(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for column in range(size):
        # The input is Vandermonde on points with pairwise unit differences.
        # Its elimination has unit pivots; this is not a generic composite-
        # ring inverse routine that assumes every nonzero element is a unit.
        pivot = next((i for i in range(column, size)
                      if gcd(rows[i][column], modulus) == 1), None)
        if pivot is None:
            raise ValueError("unit Vandermonde pivot missing")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = pow(rows[column][column], -1, modulus)
        rows[column] = [x * scale % modulus for x in rows[column]]
        for i in range(size):
            if i != column:
                multiple = rows[i][column]
                rows[i] = [(x - multiple * y) % modulus
                           for x, y in zip(rows[i], rows[column])]
    assert all(rows[i][j] == int(i == j) for i in range(size) for j in range(size))
    return [row[size:] for row in rows]


def determinant(matrix: list[list[int]]) -> int:
    """Independent fraction-free determinant with exact-divisibility checks."""
    matrix = [row[:] for row in matrix]
    size, previous, sign = len(matrix), 1, 1
    if size == 0:
        return 1
    for column in range(size - 1):
        pivot = next((i for i in range(column, size) if matrix[i][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            sign = -sign
        value = matrix[column][column]
        for i in range(column + 1, size):
            for j in range(column + 1, size):
                numerator = value * matrix[i][j] - matrix[i][column] * matrix[column][j]
                quotient, remainder = divmod(numerator, previous)
                assert not remainder
                matrix[i][j] = quotient
            matrix[i][column] = 0
        previous = value
    return sign * matrix[-1][-1]


def inverse_rational(matrix: list[list[int]]) -> list[list[Fraction]]:
    size = len(matrix)
    rows = [[Fraction(x) for x in row] + [Fraction(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(i for i in range(column, size) if rows[i][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [x / scale for x in rows[column]]
        for i in range(size):
            if i != column and rows[i][column]:
                multiple = rows[i][column]
                rows[i] = [x - multiple * y for x, y in zip(rows[i], rows[column])]
    return [row[size:] for row in rows]


def negacyclic_orbit(vector: list[int]) -> list[list[int]]:
    columns = []
    for _ in vector:
        columns.append(vector)
        vector = [-vector[-1], *vector[:-1]]
    return [list(row) for row in zip(*columns)]


class RecipientLift:
    """Construct a source preimage for one common-valuation recipient set.

    A successful instance carries an explicit norm certificate. At settings
    outside the proved theorem, failure of that certificate is not a rank
    lower bound or a proof of nonsurjectivity.
    """

    def __init__(self, fan: Fanout, recipients: tuple[int, ...]):
        if not recipients or len(set(recipients)) != len(recipients):
            raise ValueError("a nonempty set of distinct recipients is required")
        if not all(0 < j < fan.a for j in recipients):
            raise ValueError("recipient outside the nonzero source range")
        q = recipients[0] & -recipients[0]
        if any(j & -j != q for j in recipients):
            raise ValueError("this lift requires one common two-adic valuation")
        self.fan, self.recipients, self.q = fan, recipients, q
        self.dimension = fan.h // q
        self.indices = tuple(j // q for j in recipients)
        self.root = pow(fan.root, q, fan.modulus)
        self.height = q * (fan.base - 1) + 1
        self.offset = (q // 2) * (fan.base - 1)
        dimension, count, modulus = self.dimension, len(recipients), fan.modulus
        if count > dimension:
            raise ValueError("too many recipients in this cyclotomic class")
        self.evaluation = [[pow(self.root, i * j, modulus) for i in range(dimension)]
                           for j in self.indices]
        self.pivot_inverse = inverse_mod([row[:count] for row in self.evaluation], modulus)
        basis = []
        for i in range(dimension):
            row = [0] * dimension
            if i < count:
                row[i] = modulus
            else:
                row[i] = 1
                for j in range(count):
                    row[j] = (-sum(self.pivot_inverse[j][a] * self.evaluation[a][i]
                                   for a in range(count))) % modulus
            basis.append(row)
        assert abs(determinant(basis)) == modulus ** count
        reduced, self.lll_stats = lll_reduce(basis)
        assert abs(determinant(reduced)) == modulus ** count
        self.vector = reduced[0]
        self.l1 = sum(abs(x) for x in self.vector)
        self.norm_squared = sum(x * x for x in self.vector)
        assert self.l1 > 0
        assert all(sum(x * y for x, y in zip(self.vector, row)) % modulus == 0
                   for row in self.evaluation)
        # LLL's Euclidean determinant guarantee, expressed with integers.
        assert self.norm_squared ** (2 * dimension) <= (
            (1 << (dimension * (dimension - 1))) * modulus ** (4 * count))
        self.orbit = negacyclic_orbit(self.vector)
        self.orbit_determinant = determinant(self.orbit)
        assert self.orbit_determinant
        assert all(sum(abs(x) for x in row) == self.l1 for row in self.orbit)
        for column in zip(*self.orbit):
            assert all(sum(x * y for x, y in zip(column, row)) % modulus == 0
                       for row in self.evaluation)
        self.certified = self.l1 <= self.height - 1
        self.orbit_inverse = inverse_rational(self.orbit) if self.certified else None

    def lift(self, messages: tuple[int, ...]) -> tuple[int, ...]:
        if not self.certified:
            raise ValueError("short-vector norm does not certify the required digit box")
        if len(messages) != len(self.recipients) or any(not 0 <= x < self.fan.modulus for x in messages):
            raise ValueError("invalid residue messages")
        modulus = self.fan.modulus
        target = [value * pow(self.fan.phase(j), -1, modulus) % modulus
                  for j, value in zip(self.recipients, messages)]
        # Translate the signed folded interval to [0,H-1].
        target = [(value + self.offset * sum(row)) % modulus
                  for value, row in zip(target, self.evaluation)]
        count = len(target)
        initial = [sum(self.pivot_inverse[i][j] * target[j] for j in range(count)) % modulus
                   for i in range(count)] + [0] * (self.dimension - count)
        center = Fraction(self.height - 1, 2)
        displacement = [Fraction(x) - center for x in initial]
        coordinates = [sum(a * x for a, x in zip(row, displacement))
                       for row in self.orbit_inverse]
        rounded = [nearest(x) for x in coordinates]
        folded = [initial[i] - sum(a * x for a, x in zip(row, rounded))
                  for i, row in enumerate(self.orbit)]
        assert all(0 <= x < self.height for x in folded)
        assert all(sum(x * y for x, y in zip(folded, row)) % modulus == value
                   for row, value in zip(self.evaluation, target))
        source = [0] * self.fan.h
        for i, remaining in enumerate(folded):
            for block in range(self.q):
                value = min(self.fan.base - 1, remaining)
                remaining -= value
                source[i + self.dimension * block] = (value if block % 2 == 0
                                                       else self.fan.base - 1 - value)
            assert remaining == 0
        assert all(0 <= x < self.fan.base for x in source)
        assert tuple(self.fan.message_at(tuple(source), j) for j in self.recipients) == messages
        return tuple(source)

    def record(self) -> dict:
        dimension, count = self.dimension, len(self.recipients)
        guaranteed = ((self.height - 1) ** (4 * dimension)
                      >= dimension ** (2 * dimension)
                      * (1 << (dimension * (dimension - 1)))
                      * self.fan.modulus ** (4 * count))
        return {
            "h": self.fan.h, "digit_bits": self.fan.digit_bits,
            "gate_index": self.fan.gate_index, "recipients": self.recipients,
            "common_two_part": self.q, "lattice_dimension": dimension,
            "norm_certified": self.certified, "generic_lll_bound_certifies": guaranteed,
            "norm_l1": str(self.l1), "digit_height": str(self.height),
            "short_vector": [str(x) for x in self.vector],
            "orbit_determinant_bits": abs(self.orbit_determinant).bit_length(),
            **self.lll_stats,
        }


def run_checks() -> dict:
    rng = random.Random(20261101)
    cases = []
    controls = {"certified_recipient_sets": 0, "preimages": 0,
                "sharp_boundary_controls": 0, "same_valuation_sets": 0,
                "rounding_cases": 0, "all_input_sets_at_h8": 0}
    for denominator in range(1, 8):
        for numerator in range(-20, 21):
            value = Fraction(numerator, denominator)
            rounded = nearest(value)
            assert abs(value - rounded) <= Fraction(1, 2)
            controls["rounding_cases"] += 1
    schedules = [(8, (j,)) for j in range(1, 16, 2)]
    schedules += [(16, indices) for indices in
                  ((1, 3), (3, 15), (1, 17), (1, 31),
                   (1, 3, 7), (1, 9, 17), (3, 11, 27), (5, 19, 29))]
    # Common-valuation examples include a boundary certified by the actual
    # short vector even when the generic LLL norm estimate is too coarse.
    schedules += [(8, (2,)), (8, (4,)), (16, (2,)), (16, (4,)), (16, (4, 12))]
    for h, recipients in schedules:
        fan = Fanout(h, 2 * h * h, rng.randrange(2 * h))
        lift = RecipientLift(fan, recipients)
        record = lift.record()
        if all(j % 2 for j in recipients) and len(recipients) < h // 4:
            assert record["generic_lll_bound_certifies"] and lift.certified
        if lift.certified:
            targets = [tuple(0 for _ in recipients),
                       tuple(fan.modulus - 1 for _ in recipients),
                       tuple(1 for _ in recipients)]
            targets += [tuple(rng.randrange(fan.modulus) for _ in recipients) for _ in range(5)]
            for target in targets:
                source = lift.lift(target)
                assert tuple(fan.message_at(source, j) for j in recipients) == target
                controls["preimages"] += 1
            controls["certified_recipient_sets"] += 1
            if lift.q > 1:
                controls["same_valuation_sets"] += 1
        else:
            # These selected failures are also independently excluded by
            # the folded-domain cardinality, rather than by LLL failure.
            assert lift.height ** lift.dimension < fan.modulus ** len(recipients)
        if h == 8 and len(recipients) == 1 and recipients[0] % 2:
            controls["all_input_sets_at_h8"] += 1
        cases.append(record)
    for indices in ((1, 3), (1, 7), (1, 9)):
        fan = Fanout(8, 128)
        lift = RecipientLift(fan, indices)
        assert fan.base ** fan.h < fan.modulus ** len(indices)
        assert not lift.certified
        cases.append(lift.record())
        controls["sharp_boundary_controls"] += 1
    assert controls["all_input_sets_at_h8"] == 8
    return {
        "claim": "arbitrary primitive-recipient full-range threshold; no factoring score increase",
        "source_revision": "c91cc3f4947caca8977e81a5d14244f73b5ae2b1",
        "algorithm": "exact-rational delta=3/4 LLL and negacyclic-orbit rounding",
        "controls": controls, "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run_checks()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["controls"], sort_keys=True))


if __name__ == "__main__":
    main()
