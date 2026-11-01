#!/usr/bin/env python3
"""Exact key counts and weighted dominance controls for the cofactor quota."""

from __future__ import annotations

import argparse
import json
from math import gcd, lcm, prod
from pathlib import Path
import random

from cofactor_collision_relations import (RationalRowBasis, partial_leading_word,
                                         verify_partial_identity)
from prime_form_factor_base import build_prime_form_factor_base
from relation_two_saturation import first_useful_relation_top
from weighted_norm_relations import least_power_above


def word_norm(word: dict[int, int]) -> int:
    return prod(prime ** abs(exponent) for prime, exponent in word.items())


def bounded_word(rng: random.Random, primes: tuple[int, ...], bound: int) -> dict[int, int]:
    word: dict[int, int] = {}
    for _ in range(30):
        prime = rng.choice(primes)
        candidate = dict(word)
        candidate[prime] = candidate.get(prime, 0) + rng.choice((-1, 1))
        if word_norm(candidate) <= bound:
            word = candidate
    return word


def rational_span_controls() -> int:
    """Compare retained-row extraction with fully enumerated saturation images."""
    rng = random.Random(401501)
    for _ in range(300):
        dimension = rng.randrange(1, 6)
        rank = rng.randrange(1, dimension + 1)
        moduli = (rng.randrange(2, 10), rng.randrange(2, 10))
        identity = (0, 0)

        def multiply(left, right):
            return tuple((x + y) % modulus for x, y, modulus in zip(left, right, moduli))

        def inverse(value):
            return tuple(-x % modulus for x, modulus in zip(value, moduli))

        generators = [tuple(rng.randrange(modulus) for modulus in moduli)
                      for _ in range(dimension)]

        def evaluate(row):
            return tuple(sum(coefficient * value[j] for coefficient, value in zip(row, generators))
                         % moduli[j] for j in range(2))

        # First rows of a unimodular matrix form a known saturated lattice.
        saturated_basis = [[int(i == j) for j in range(dimension)] for i in range(dimension)]
        if dimension > 1:
            for _ in range(12):
                first, second = rng.sample(range(dimension), 2)
                scalar = rng.randrange(-2, 3)
                saturated_basis[first] = [a + scalar * b for a, b in
                                          zip(saturated_basis[first], saturated_basis[second])]
        saturated_basis = saturated_basis[:rank]
        images = [evaluate(row) for row in saturated_basis]
        subgroup, frontier = {identity}, [identity]
        while frontier:
            point = frontier.pop()
            for value in images:
                candidate = multiply(point, value)
                if candidate not in subgroup:
                    subgroup.add(candidate)
                    frontier.append(candidate)
        tops = [point for point in subgroup if point != identity and multiply(point, point) == identity]
        rejected = tops[0] if tops and rng.randrange(2) else None
        useful = {point for point in tops if point != rejected}
        relations = []
        for row, value in zip(saturated_basis, images):
            order = lcm(*(modulus // gcd(coordinate, modulus)
                          for coordinate, modulus in zip(value, moduli)))
            relations.append([order * coefficient for coefficient in row])
        candidates = []
        for row in relations:
            candidates.append(row)
            scalar = 1 << rng.randrange(1, 7)
            candidates.append([scalar * coefficient for coefficient in row])
        for _ in range(12):
            coefficients = [rng.randrange(-5, 6) for _ in range(rank)]
            candidates.append([sum(coefficients[j] * relations[j][i] for j in range(rank))
                               for i in range(dimension)])
        rng.shuffle(candidates)
        retained = RationalRowBasis(dimension)
        for row in candidates:
            assert evaluate(row) == identity
            retained.add(row)
        assert len(retained.rows) == rank

        def inspect(top):
            return top if top != rejected else None

        top, certificate, _ = first_useful_relation_top(
            generators, retained.rows, identity, multiply, inverse, inspect)
        assert (top is not None) == bool(useful)
        if top is not None:
            assert top in useful and certificate == top
    return 300


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    key_controls = 0
    for n in (45, 77, 105, 221, 959):
        delta = -n if n % 4 == 3 else -4 * n
        for bound in range(1, 41):
            keys = {None}
            for residual in range(3, bound + 1, 2):
                if gcd(residual, 2 * delta) != 1:
                    continue
                for middle in range(delta % 2, 2 * residual, 2):
                    if (middle * middle - delta) % (4 * residual) == 0:
                        canonical = min(middle, (-middle) % (2 * residual))
                        keys.add((residual, canonical))
            assert len(keys) <= bound * bound.bit_length()
            key_controls += 1
    rng = random.Random(15911)
    primes = (2, 3, 5, 7, 11)
    dominance_controls = 0
    for _ in range(1000):
        height = rng.randrange(1, 100001)
        prime = rng.choice(primes)
        first = least_power_above(prime, height) + rng.randrange(10)
        gap = least_power_above(prime, height * height)
        assert prime ** gap > height * height
        if gap:
            assert prime ** (gap - 1) <= height * height
        left, right = bounded_word(rng, primes, height), bounded_word(rng, primes, height)
        sign = rng.choice((-1, 1))
        opposite = bool(rng.randrange(2))
        other_sign = -sign if opposite else sign
        second = first + (rng.randrange(1, 10) if opposite else gap)
        row = {p: sign * left.get(p, 0) - other_sign * right.get(p, 0) for p in primes}
        row[prime] += other_sign * second - sign * first
        assert prime ** abs(row[prime]) > prod(p ** abs(row[p]) for p in primes if p != prime)
        dominance_controls += 1
    # Equality p^gap=H^2 cannot replace the strict threshold: the row can vanish.
    equality = {"prime": 2, "height": 16, "gap": 8, "first_exponent": 5,
                "left_word": {2: -4}, "right_word": {2: 4}, "cancelled_row": {2: 0}}
    assert 2 ** equality["gap"] == equality["height"] ** 2
    assert (equality["first_exponent"] + equality["gap"] - 4
            - (equality["first_exponent"] + 4)) == 0
    n = 959
    records = {item.prime: (item.root, item.positive)
               for item in build_prime_form_factor_base(-n, 2)["records"]}
    fixtures = []
    for form in ((15, 1, 16), (15, -1, 16), (15, 11, 18)):
        split = partial_leading_word(n, form, sorted(records), records)
        assert verify_partial_identity(n, form, split, records)
        fixtures.append({"form": form, "key": split["key"], "orientation": split["orientation"]})
    assert fixtures[0]["key"] == fixtures[1]["key"] != fixtures[2]["key"]
    assert (fixtures[0]["orientation"], fixtures[1]["orientation"]) == (1, -1)
    result = {"claim": "A finite per-generator admissible quota gives a full-rank cofactor relation family",
              "literal_key_count_controls": key_controls,
              "exact_weighted_dominance_controls": dominance_controls,
              "finite_abelian_saturation_image_controls": rational_span_controls(),
              "strict_threshold_control": equality,
              "composite_key_fixture": {"integer": n, "records": fixtures},
              "scope": "The time to accumulate the per-generator admissible quota remains unproved."}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
