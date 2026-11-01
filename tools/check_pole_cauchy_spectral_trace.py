#!/usr/bin/env python3
"""Exact finite controls for spectral projection and leading inverse structure."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random

from check_pole_cauchy_conditioning import inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from interval_unit_jets import QuadraticJets
from pole_cauchy_spectral_trace import spectral_trace_contract, count_from_supplied_spectral_traces
from pole_cauchy_spectral_trace import leading_inverse_partner
from pole_cauchy_spectral_trace import padded_trace_contract
from pole_cauchy_residue_blocks import selector_rounding_coefficients

ROOT = Path(__file__).resolve().parents[1]


def identity(size):
    return [[(int(i == j), 0) for j in range(size)] for i in range(size)]


def matrix_add(ring, a, b, scalar=1):
    return [[tuple((x + scalar * y) % ring.modulus for x, y in zip(u, v))
             for u, v in zip(row_a, row_b)] for row_a, row_b in zip(a, b)]


def matrix_multiply(ring, a, b):
    size = len(a)
    result = [[(0, 0)] * size for _ in range(size)]
    for i in range(size):
        for j in range(size):
            real = imaginary = 0
            for ell in range(size):
                x, y = ring.pair_product(a[i][ell], b[ell][j])
                real += x
                imaginary += y
            result[i][j] = real % ring.modulus, imaginary % ring.modulus
    return result


def trace(ring, matrix):
    return tuple(sum(matrix[i][i][coordinate] for i in range(len(matrix))) % ring.modulus
                 for coordinate in (0, 1))


def projector(ring, matrix, precision):
    size = len(matrix)
    value = [[(0, 0)] * size for _ in range(size)]
    power = identity(size)
    coefficients = selector_rounding_coefficients(precision)
    traces = {}
    for degree in range(1, 2 * precision):
        power = matrix_multiply(ring, power, matrix)
        if degree in coefficients:
            traces[degree] = trace(ring, power)
            value = matrix_add(ring, value, power, coefficients[degree])
    return value, traces


def controls():
    rng = random.Random(2026110119)
    arbitrary_projectors = power_agreements = pole_projectors = displacement_entries = 0
    for index in range(48):
        size, rank, precision = 1 + index % 5, index % (2 + index % 5), 1 + index % 6
        rank = min(rank, size)
        ring = QuadraticJets(precision, 0)
        base = [[((int(i == j and i < rank) + 2 * rng.randrange(-3, 4)) % ring.modulus,
                  2 * rng.randrange(-3, 4) % ring.modulus) for j in range(size)] for i in range(size)]
        nilpotent = [[((rng.randrange(4) if i < j else 0), (rng.randrange(4) if i < j else 0))
                      for j in range(size)] for i in range(size)]
        change = matrix_add(ring, identity(size), nilpotent)
        inverse_change, power = identity(size), identity(size)
        for degree in range(1, size):
            power = matrix_multiply(ring, power, nilpotent)
            inverse_change = matrix_add(ring, inverse_change, power, (-1) ** degree)
        matrix = matrix_multiply(ring, matrix_multiply(ring, change, base), inverse_change)
        projected, _ = projector(ring, matrix, precision)
        assert matrix_multiply(ring, projected, projected) == projected
        assert trace(ring, projected) == (rank % ring.modulus, 0)
        powered = matrix
        for _ in range(precision - 1):
            powered = matrix_multiply(ring, powered, powered)
        assert powered == projected
        arbitrary_projectors += 1
        power_agreements += 1
    examples = []
    padded_trace_controls = padded_entry_sum_controls = 0
    for index in range(36):
        size = 1 + index % 5
        aa = tuple(range(1, 2 * size, 2))
        bb = tuple(range(3, 2 * size + 2, 2))
        n = rng.randrange(15, 100, 2)
        k = max(abs(n - a * b) for a in aa for b in bb).bit_length() + 1
        matrix = matrix_at(n, aa, bb, k)
        exact_s = [[sum((matrix[i][ell] * matrix[j][ell] for ell in range(size)), Quad())
                    for j in range(size)] for i in range(size)]
        unit_vector = [sum((matrix[i][j] / b for j, b in enumerate(bb)), Quad())
                       for i in range(size)]
        t = Quad(Fraction(0), Fraction(1 << k))
        for i, a in enumerate(aa):
            for j, other in enumerate(aa):
                assert (a - other) * exact_s[i][j] == t * (unit_vector[i] - unit_vector[j])
                displacement_entries += 1
        precision = size.bit_length()
        ring = QuadraticJets(precision, 0)
        modular_s = [[mod_pair(value, precision) for value in row] for row in exact_s]
        projected, traces = projector(ring, modular_s, precision)
        assert matrix_multiply(ring, projected, projected) == projected
        actual = count_from_supplied_spectral_traces(size, traces)
        count = sum(a * b == n for a in aa for b in bb)
        assert actual["count"] == count
        nonzero_differences = [abs(n - a * b) for a in aa for b in bb if n != a * b]
        bound = max(((value & -value).bit_length() - 1 for value in nonzero_differences), default=0)
        padded = padded_trace_contract(size, bound)
        padded_matrix = matrix_at(n, aa, bb, padded["one_trace_displacement"])
        first_trace = sum((value * value for row in padded_matrix for value in row), Quad())
        assert mod_pair(first_trace, precision) == (count, 0)
        sum_matrix = matrix_at(n, aa, bb, padded["entry_sum_displacement"])
        entry_sum = sum((value for row in sum_matrix for value in row), Quad())
        assert mod_pair(entry_sum, precision) == (count, 0)
        padded_trace_controls += 1
        padded_entry_sum_controls += 1
        pole_projectors += 1
        if len(examples) < 8:
            examples.append({"n": n, "size": size, "k": k, "precision": precision,
                             "factor_pairs": count, "maximum_matrix_power": actual["maximum_matrix_power"],
                             "padded_one_trace": padded})
    leading_matrices = leading_entries = involutions = resolvent_controls = 0
    leading_examples = []
    for size in (1, 2, 4, 8, 16):
        for n in (15, 21, 31, 105):
            aa = tuple(range(1, 2 * size, 2))
            k = max(abs(n - a * b) for a in aa for b in aa).bit_length() + 2
            s = size.bit_length() - 1
            loss = k - s - 1
            dense_inverse = inverse(matrix_at(n, aa, aa, k))
            ring = QuadraticJets(1, 0)
            leading = [[ring.pair_product((0, 1), mod_pair(value * (1 << loss), 1))
                        for value in row] for row in dense_inverse]
            predicted = [[(0, 0)] * size for _ in range(size)]
            for i, a in enumerate(aa):
                edge = leading_inverse_partner(n, size, a)
                if edge["retained"]:
                    predicted[(edge["partner"] - 1) // 2][i] = (1, 0)
            assert leading == predicted
            square = matrix_multiply(ring, leading, leading)
            assert matrix_multiply(ring, square, leading) == leading
            assert all(square[i][j] == ((int(any(value != (0, 0) for value in leading[i])), 0) if i == j else (0, 0))
                       for i in range(size) for j in range(size))
            # At x=w, the scalar inverse of 1-x^2 exists in F4.
            x = (0, 1)
            x2 = ring.pair_product(x, x)
            scalar_inverse = ring.inverse(ring.constant(1 ^ x2[0], x2[1]))[0]
            first = ring.pair_product(x, scalar_inverse)
            second = ring.pair_product(x2, scalar_inverse)
            candidate = identity(size)
            for i in range(size):
                for j in range(size):
                    a, b = ring.pair_product(first, leading[i][j]), ring.pair_product(second, square[i][j])
                    candidate[i][j] = (candidate[i][j][0] ^ a[0] ^ b[0], candidate[i][j][1] ^ a[1] ^ b[1])
            pencil = [[((int(i == j) ^ ring.pair_product(x, leading[i][j])[0]),
                         ring.pair_product(x, leading[i][j])[1]) for j in range(size)] for i in range(size)]
            assert matrix_multiply(ring, pencil, candidate) == identity(size)
            retained = sum(any(value != (0, 0) for value in row) for row in leading)
            leading_examples.append({"n": n, "size": size, "k": k, "inverse_loss": loss,
                                     "leading_inverse_rank": retained,
                                     "true_factor_pairs": sum(a * b == n for a in aa for b in aa)})
            leading_matrices += 1
            leading_entries += size * size
            involutions += 1
            resolvent_controls += 1
    large_contract = spectral_trace_contract(1 << 200)
    assert large_contract["precision"] == 201 and large_contract["maximum_matrix_power"] == 401
    return {"arbitrary_local_ring_projectors": arbitrary_projectors,
            "power_of_two_projection_agreements": power_agreements,
            "pole_matrix_count_projections": pole_projectors,
            "padded_one_trace_controls": padded_trace_controls,
            "padded_entry_sum_controls": padded_entry_sum_controls,
            "exact_displacement_entries": displacement_entries,
            "leading_inverse_matrices": leading_matrices,
            "leading_inverse_entries": leading_entries,
            "partial_involution_identities": involutions, "resolvent_identities": resolvent_controls,
            "examples": examples, "leading_examples": leading_examples,
            "large_supplied_trace_contract": large_contract}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"claim": "A degree-linear spectral trace interface and an explicit leading inverse partial involution",
              "controls": controls(), "factoring_exponent_improved": False,
              "source_sha256": {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "pole_cauchy_spectral_trace.py",
                                             "pole_cauchy_residue_blocks.py")},
              "scope": "Supplied global traces recover the count; evaluating those traces implicitly remains open."}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report["controls"].items() if isinstance(value, int)}, indent=2))


if __name__ == "__main__":
    main()
