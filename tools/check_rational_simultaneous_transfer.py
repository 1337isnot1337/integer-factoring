#!/usr/bin/env python3
"""Finite algebra and exact cutoff checks for the varying-ring transfer.

The small tensor fixtures test the reduction and type-channel identities.
They are not claimed to supply the asymptotic contracting family assumed
by the conditional theorem.
"""

import argparse
from fractions import Fraction
from itertools import combinations_with_replacement, product
import json
from math import gcd, isqrt
from pathlib import Path
from random import Random


def safe_modulus(denominators):
    denominator_product = 1
    for denominator in denominators:
        if denominator == 0:
            raise ValueError("zero denominator")
        denominator_product *= abs(denominator)
    return 2 * denominator_product + 1, denominator_product


def reduce_fraction(value, modulus):
    assert gcd(value.denominator, modulus) == 1
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def modulus_checks():
    cases = homomorphisms = vanished = composites = 0
    rng = Random(202610300755)
    for denominator_product in range(1, 1501):
        modulus, actual_product = safe_modulus((denominator_product,))
        assert actual_product == denominator_product
        assert modulus == 2 * denominator_product + 1 and modulus % 2 == 1
        assert gcd(modulus, denominator_product) == 1
        composites += any(modulus % divisor == 0 for divisor in range(2, isqrt(modulus) + 1))
        cases += 1
    for _ in range(1200):
        denominators = [rng.randrange(1, 8) for _ in range(rng.randrange(1, 5))]
        modulus, denominator_product = safe_modulus(denominators)
        values = [Fraction(rng.randrange(-40, 41), denominator) for denominator in denominators]
        left, right = values[0], values[-1]
        a, b = reduce_fraction(left, modulus), reduce_fraction(right, modulus)
        assert reduce_fraction(left + right, modulus) == (a + b) % modulus
        assert reduce_fraction(left * right, modulus) == a * b % modulus
        assert all(denominator_product % denominator == 0 for denominator in denominators)
        homomorphisms += 2
        zero_image = Fraction(modulus, denominators[0])
        assert zero_image != 0 and reduce_fraction(zero_image, modulus) == 0
        vanished += 1
    assert vanished > 0
    assert composites > 0
    return {"safe_modulus_cases": cases, "composite_moduli": composites,
            "homomorphism_identities": homomorphisms,
            "nonzero_rationals_reducing_to_zero": vanished}


def local_factor(piece, first_side, pair):
    x, y, coefficient = piece
    if pair != (x, y):
        return Fraction(0)
    return coefficient if first_side else Fraction(1)


def local_entry(piece, coords):
    x, y, coefficient = piece
    a, b, u, v = coords
    return coefficient if (a, b, u, v) == (x, y, y, x) else Fraction(0)


def digits(value, blocks):
    return tuple((value >> (blocks - 1 - position)) & 1 for position in range(blocks))


def type_checks():
    scalar_factors = crossing_entries = grouped_entries = channels = zero_pieces = 0
    rings = []
    fixture_modulus = 2 * (2 ** 16) + 1
    for coefficients in ((Fraction(fixture_modulus, 2), Fraction(2 - fixture_modulus, 2)),
                         (Fraction(1, 2), Fraction(1, 2))):
        pieces = tuple((x, y, coefficient) for x in range(2) for y in range(2)
                       for coefficient in coefficients)
        denominators = []
        for piece in pieces:
            for cut, first_side, a, b in product(range(2), (False, True), range(2), range(2)):
                value = local_factor(piece, first_side, (a, b))
                denominators.append(value.denominator)
                scalar_factors += 1
        modulus, denominator_product = safe_modulus(denominators)
        rings.append({"modulus": modulus, "denominator_product_bits": denominator_product.bit_length()})
        assert modulus == fixture_modulus and modulus % 3 == 0
        for coords in product(range(2), repeat=4):
            a, b, u, v = coords
            total = Fraction(0)
            for piece in pieces:
                expected = local_entry(piece, coords)
                for left_pair, right_pair in (((a, b), (v, u)), ((a, u), (v, b))):
                    left = local_factor(piece, True, left_pair)
                    right = local_factor(piece, False, right_pair)
                    assert left * right == expected
                    assert reduce_fraction(left, modulus) * reduce_fraction(right, modulus) % modulus == (
                        reduce_fraction(expected, modulus))
                total += expected
            assert total == int(a == v and b == u)
            crossing_entries += 1
        zero_pieces += sum(reduce_fraction(piece[2], modulus) == 0 for piece in pieces)

        blocks = 2
        dimension = 1 << blocks
        total = [0] * (dimension ** 4)
        for ordered_type in combinations_with_replacement(range(len(pieces)), blocks):
            counts = tuple(ordered_type.count(i) for i in range(len(pieces)))
            words = tuple(word for word in product(range(len(pieces)), repeat=blocks)
                          if tuple(word.count(i) for i in range(len(pieces))) == counts)
            channels += len(words)
            for index, coords in enumerate(product(range(dimension), repeat=4)):
                split = tuple(digits(value, blocks) for value in coords)
                expected = Fraction(0)
                for word in words:
                    term = Fraction(1)
                    for position, piece in enumerate(word):
                        term *= local_entry(pieces[piece], tuple(value[position] for value in split))
                    expected += term
                reduced = reduce_fraction(expected, modulus)
                total[index] = (total[index] + reduced) % modulus
                for cut in range(2):
                    actual = 0
                    for word in words:
                        left = right = 1
                        for position, piece in enumerate(word):
                            a, b, u, v = (value[position] for value in split)
                            lp, rp = ((a, b), (v, u)) if cut == 0 else ((a, u), (v, b))
                            left = left * reduce_fraction(local_factor(pieces[piece], True, lp), modulus) % modulus
                            right = right * reduce_fraction(local_factor(pieces[piece], False, rp), modulus) % modulus
                        actual = (actual + left * right) % modulus
                    assert actual == reduced
                    grouped_entries += 1
        for index, (a, b, u, v) in enumerate(product(range(dimension), repeat=4)):
            assert total[index] == int(a == v and b == u)
            crossing_entries += 1
    assert zero_pieces == 4
    return {"enumerated_scalar_factors": scalar_factors,
            "type_factorization_entries": grouped_entries,
            "type_channels": channels, "crossing_entries": crossing_entries,
            "vanishing_pieces_preserving_crossing": zero_pieces, "fixture_rings": rings}


def composite_character_checks():
    orthogonality = tagged = 0
    moduli = (3, 9, 15, 21, 35, 131073)
    for modulus in moduli:
        for k in range(5):
            size = 1 << k
            normalization = pow(size, -1, modulus)
            for actual, bucket in product(range(size), repeat=2):
                value = sum((-1 if ((character & (actual ^ bucket)).bit_count() & 1) else 1)
                            for character in range(size)) % modulus
                assert value * normalization % modulus == int(actual == bucket)
                orthogonality += 1
        for number in (6, 12, 15, 21, 35):
            bits = number.bit_length()
            witnesses = [(factor, number // factor) for factor in range(1, number + 1)
                         if number % factor == 0]
            size = 4
            basis_hashes = tuple((position * position + position + 1) % size
                                for position in range(bits))
            hashes = []
            for first, _ in witnesses:
                value = 0
                for position, basis in enumerate(basis_hashes):
                    if first >> position & 1:
                        value ^= basis
                hashes.append(value)
            normalization = pow(size, -1, modulus)
            for bucket, family, target in product(range(size), range(2), range(bits)):
                scalar = tag = 0
                for character in range(size):
                    weight = -1 if ((character & bucket).bit_count() & 1) else 1
                    evaluated_scalar = evaluated_tag = 0
                    for pair in witnesses:
                        phase = 1
                        for position, basis in enumerate(basis_hashes):
                            if pair[0] >> position & 1:
                                if (character & basis).bit_count() & 1:
                                    phase = -phase
                        evaluated_scalar += phase
                        evaluated_tag += phase * ((pair[family] >> target) & 1)
                    scalar += weight * evaluated_scalar
                    tag += weight * evaluated_tag
                scalar = scalar * normalization % modulus
                tag = tag * normalization % modulus
                expected_scalar = sum(value == bucket for value in hashes) % modulus
                expected_tag = sum(((pair[family] >> target) & 1)
                                   for pair, value in zip(witnesses, hashes) if value == bucket) % modulus
                assert (scalar, tag) == (expected_scalar, expected_tag)
                tagged += 1
    return {"binary_character_orthogonality": orthogonality,
            "paired_tag_bucket_projections": tagged, "moduli": moduli}


def number_field_checks():
    """Q(sqrt(2)) and R[X]/(X^2-2), with no target-ring inversion."""
    def add(left, right):
        return left[0] + right[0], left[1] + right[1]

    def multiply(left, right):
        return (left[0] * right[0] + 2 * left[1] * right[1],
                left[0] * right[1] + left[1] * right[0])

    def inverse(value):
        norm = value[0] * value[0] - 2 * value[1] * value[1]
        assert norm != 0
        return value[0] / norm, -value[1] / norm

    def reduce(value, modulus):
        return tuple(reduce_fraction(coordinate, modulus) for coordinate in value)

    def ring_multiply(left, right, modulus):
        return tuple(coordinate % modulus for coordinate in multiply(left, right))

    rng = Random(202610300807)
    identities = zero_divisors = character_checks = 0
    for _ in range(400):
        left = tuple(Fraction(rng.randrange(-8, 9), rng.randrange(1, 7)) for _ in range(2))
        right = tuple(Fraction(rng.randrange(-8, 9), rng.randrange(1, 7)) for _ in range(2))
        if left == (0, 0):
            left = (Fraction(1), Fraction(1))
        reciprocal = inverse(left)
        values = left + right + reciprocal
        modulus, _ = safe_modulus(value.denominator for value in values)
        a, b, ai = reduce(left, modulus), reduce(right, modulus), reduce(reciprocal, modulus)
        assert reduce(add(left, right), modulus) == tuple(c % modulus for c in add(a, b))
        assert reduce(multiply(left, right), modulus) == ring_multiply(a, b, modulus)
        assert ring_multiply(a, ai, modulus) == (1, 0)
        identities += 3

    for modulus, root in ((7, 3), (49, 10)):
        # f splits in the target ring, even though X^2-2 is irreducible
        # over Q; the second coefficient modulus is itself composite.
        left, right = ((-root) % modulus, 1), (root % modulus, 1)
        assert left != (0, 0) and right != (0, 0)
        assert ring_multiply(left, right, modulus) == (0, 0)
        zero_divisors += 1
        value = (2, 1)
        normalization = pow(4, -1, modulus)
        for actual, bucket in product(range(4), repeat=2):
            scalar = sum((-1 if ((character & (actual ^ bucket)).bit_count() & 1) else 1)
                         for character in range(4))
            got = tuple(scalar * normalization * coordinate % modulus for coordinate in value)
            assert got == (value if actual == bucket else (0, 0))
            character_checks += 1
    return {"rational_coordinate_identities": identities,
            "split_polynomial_zero_divisor_examples": zero_divisors,
            "quotient_ring_character_identities": character_checks}


def ceil_log2(value):
    return (value - 1).bit_length()


def block_size(label, constant):
    return max(1, min(label, ceil_log2(label + 2) // constant))


def cutoff_checks():
    rng = Random(202610300756)
    cases = block_cases = envelope_cases = 0
    inputs = list(range(2, 3000)) + [1 << power for power in (32, 64, 128, 512, 2048, 8192)]
    for n in inputs:
        volume_bound = (1 << 16) * (n + 2) ** 3
        for height_constant in (1, 3, 7):
            blocking_constant = 96 * height_constant
            logarithm = ceil_log2(volume_bound + 2)
            cutoff = max(1, (logarithm + blocking_constant - 1) // blocking_constant)
            denominator_bits_bound = 1 << (height_constant * cutoff)
            # This is the exact integer form of b(n)=O(n^(1/32)).
            assert denominator_bits_bound ** 32 <= (1 << (32 * height_constant + 6)) * (n + 2)
            for label in (1, 2, volume_bound // 2, volume_bound):
                assert block_size(label, blocking_constant) <= cutoff
                block_cases += 1
            cases += 1
    for ell in range(17):
        n = 1 << (1 << ell)
        coarse_volume = 1296 * n * 4 ** ell
        assert coarse_volume <= (1 << 16) * (n + 2) ** 3
        envelope_cases += 1
    return {"cutoff_and_height_bounds": cases, "block_coverage_cases": block_cases,
            "crude_volume_envelopes": envelope_cases}


def verify():
    return {"status": "pass", "odd_modulus": modulus_checks(),
            "rational_type_reduction": type_checks(), "integer_bounds": cutoff_checks(),
            "composite_character_projection": composite_character_checks(),
            "fixed_number_field": number_field_checks(),
            "scope": "Finite interface checks; the contracting rational family remains a theorem hypothesis."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    evidence = verify()
    rendered = json.dumps(evidence, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
