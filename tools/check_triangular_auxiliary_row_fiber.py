#!/usr/bin/env python3
"""Exact common-height auxiliary fiber for the four-letter reflection seed.

The linear classification concerns two specified row lattices.  The word
argument applies to arbitrary auxiliary dimension and noncommuting constant
auxiliary coefficients.  It is not a general lower bound on crossing rank.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

import check_polynomial_crossing_interpolation as interpolation
from check_schur_gaussian_reduction import polynomial_multiply as multiply


ROOT = Path(__file__).resolve().parent.parent
HEIGHT = (0, 3, 1, 4)
LABELS = tuple(tuple((mask >> (2 * leg)) & 3 for leg in range(4))
               for mask in range(256))
ZERO_HEIGHTS = (HEIGHT, HEIGHT, (0, 0, 0, 0), (0, 0, 0, 0))
INFINITY_HEIGHTS = ((0, 0, 0, 0), (0, 0, 0, 0),
                    tuple(value + 2 for value in HEIGHT),
                    tuple(value + 2 for value in HEIGHT))
EXPECTED_ZERO = (4, 5, 7, 8, 9, 10, 11, 12,
                 12, 13, 14, 15, 16, 17, 19, 20)
EXPECTED_INFINITY = tuple(value - 4 for value in EXPECTED_ZERO)


def reflection_seed():
    """The M=4 reflection Gaussian, after removing even ramification.

    This is a direct eight-variable construction, so the checker does not
    need a different branch's reflection helper or a saved coefficient file.
    """
    kernel = [[0] * 8 for _ in range(8)]
    mode_heights = (3, 1, 3, 1, 0, 0, 0, 0)

    def put(leg1, i, leg2, j, powers):
        a, b = 2 * leg1 + i, 2 * leg2 + j
        assert a != b and kernel[a][b] == 0
        value = 0
        for power in powers:
            exponent = power + mode_heights[a] + mode_heights[b]
            assert exponent >= 0
            value ^= 1 << exponent
        kernel[a][b] = kernel[b][a] = value

    for leg in range(4):
        put(leg, 0, leg, 1, (0, 4 if leg >= 2 else -4))
    for i, j in product(range(2), repeat=2):
        exponent = 2 * (j & ~i)
        put(3, i, 1, j, (exponent, exponent + 1))
        put(2, i, 0, j, (exponent, exponent + 1))
        exponent = 2 * (1 - (i | j))
        put(3, i, 2, j, (exponent, exponent + 1))
        exponent = -4 + 2 * (i & j)
        put(1, i, 0, j, (exponent, exponent + 1))
        if i == j:
            put(3, i, 0, j, (0,))
            put(1, i, 2, j, (0,))

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return 1
        if mask.bit_count() % 2:
            return 0
        first = mask & -mask
        i = first.bit_length() - 1
        rest = mask ^ first
        value = 0
        for j in range(8):
            if rest >> j & 1:
                value ^= multiply(kernel[i][j], pfaffian(rest ^ (1 << j)))
        return value

    values = []
    for mask, labels in enumerate(LABELS):
        shift = sum(height for i, height in enumerate(mode_heights)
                    if mask >> i & 1)
        value = pfaffian(mask) << 4
        assert not value & ((1 << shift) - 1)
        value >>= shift
        x, y, u, v = labels
        assert value.bit_count() % 2 == int(x == v and y == u)
        values.append(value)
    assert max(value.bit_length() - 1 for value in values) == 12
    return tuple(values)


def encode_column_targets(row, shift, targets):
    result = offset = 0
    for value, width in zip(row, targets):
        result |= ((value << shift) & ((1 << width) - 1)) << offset
        offset += width
    return result


def row_lattice_system(values):
    endpoint_data = []
    for heights, cut, infinity, expected in (
            (ZERO_HEIGHTS, 0, False, EXPECTED_ZERO),
            (INFINITY_HEIGHTS, 1, True, EXPECTED_INFINITY)):
        matrix = interpolation.gauged_matrix(values, heights, cut, infinity)
        orders, left, right = interpolation.smith_jet_changes(matrix, 32)
        assert orders == expected
        endpoint_data.append((orders, right))

    variables, columns = [], []
    first_width = 16 * sum(endpoint_data[0][0])
    for mask, labels in enumerate(LABELS):
        h0 = interpolation.height(labels, ZERO_HEIGHTS)
        hi = interpolation.height(labels, INFINITY_HEIGHTS)
        i0, j0 = interpolation.flatten_position(labels, 0)
        ii, ji = interpolation.flatten_position(labels, 1)
        row0 = tuple(value ^ (value << 1)
                     for value in endpoint_data[0][1][j0])
        rowi = tuple(value ^ (value << 1)
                     for value in endpoint_data[1][1][ji])
        for power in range(-h0, hi):
            first = encode_column_targets(row0, power + h0, endpoint_data[0][0])
            second = encode_column_targets(rowi, hi - power - 1, endpoint_data[1][0])
            first <<= i0 * sum(endpoint_data[0][0])
            second <<= first_width + ii * sum(endpoint_data[1][0])
            variables.append((mask, power))
            columns.append(first | second)
    basis, kernel = interpolation.column_elimination(columns)
    assert (len(variables), len(basis), len(kernel)) == (3072, 3058, 14)
    return variables, columns, kernel


def classify_fiber(values, variables, columns, kernel):
    lookup = {variable: i for i, variable in enumerate(variables)}
    generators, witnesses = [], []
    for target, source in product(range(4), repeat=2):
        for power in range(HEIGHT[source] - HEIGHT[target], 0):
            witness = 0
            for mask, labels in enumerate(LABELS):
                if labels[0] != target:
                    continue
                value = values[(mask & ~3) | source]
                for degree in range(value.bit_length()):
                    if value >> degree & 1:
                        witness ^= 1 << lookup[mask, power + degree]
            image = 0
            for i, column in enumerate(columns):
                if witness >> i & 1:
                    image ^= column
            assert image == 0
            generators.append((target, source, power))
            witnesses.append(witness)
    assert len(generators) == 14
    basis, dependencies = interpolation.column_elimination(witnesses)
    assert len(basis) == 14 and not dependencies
    for witness in kernel:
        assert interpolation.solve_in_basis(basis, witness) is not None
    return generators


def auxiliary_product(left, right):
    return tuple(tuple(multiply(left[i][0], right[0][j]) ^
                       multiply(left[i][1], right[1][j]) for j in range(2))
                 for i in range(2))


IDENTITY_AUXILIARY = ((1, 0), (0, 1))
ZERO_AUXILIARY = ((0, 0), (0, 0))


def noncommuting_leg_map(infinity):
    """A two-direction control with a genuine auxiliary return path."""
    table = {(a, b): IDENTITY_AUXILIARY if a == b else ZERO_AUXILIARY
             for a, b in product(range(4), repeat=2)}
    for target, source, power, i, j in ((2, 0, -1, 0, 1),
                                      (1, 2, -2, 1, 0)):
        entry = [list(row) for row in table[target, source]]
        for exponent in (power, power + 1):
            degree = -exponent if infinity else exponent + HEIGHT[target] - HEIGHT[source]
            assert degree >= 0
            entry[i][j] ^= 1 << degree
        table[target, source] = tuple(tuple(row) for row in entry)
    first, second = table[2, 0], table[1, 2]
    assert auxiliary_product(first, second) != auxiliary_product(second, first)
    return table


def word_operator(table, length):
    words = tuple(product(range(4), repeat=length))
    result = []
    for output in words:
        row = []
        for input_word in words:
            value = IDENTITY_AUXILIARY
            for a, b in zip(output, input_word):
                value = auxiliary_product(value, table[a, b])
            scalar = value[0][0]
            assert scalar == 1 if output == input_word else (
                not scalar or sum(HEIGHT[a] for a in output) >
                sum(HEIGHT[b] for b in input_word))
            row.append(scalar)
        result.append(row)
    return words, result


def verify_word_units(table):
    entries = nilpotence_checks = 0
    for length in range(1, 4):
        words, matrix = word_operator(table, length)
        n = len(words)
        identity = [[int(i == j) for j in range(n)] for i in range(n)]
        nilpotent = [[matrix[i][j] ^ identity[i][j] for j in range(n)]
                     for i in range(n)]
        if length == 2:
            assert any(value for row in nilpotent for value in row)
        power = identity
        inverse = [[0] * n for _ in range(n)]
        precision = 16 * length * length + 1
        for _ in range(4 * length + 1):
            inverse = [[a ^ b for a, b in zip(rowa, rowb)]
                       for rowa, rowb in zip(inverse, power)]
            power = interpolation.polynomial_matrix_product(power, nilpotent, precision)
        assert not any(value for row in power for value in row)
        assert interpolation.polynomial_matrix_product(matrix, inverse, precision) == identity
        entries += n * n
        nilpotence_checks += 1
    return entries, nilpotence_checks


def verify_two_site_identity(values, table, infinity):
    heights = INFINITY_HEIGHTS if infinity else ZERO_HEIGHTS
    base = {}
    for labels, value in zip(LABELS, values):
        h = interpolation.height(labels, heights)
        transformed = 0
        for degree in range(value.bit_length()):
            if value >> degree & 1:
                exponent = h - degree if infinity else h + degree
                assert exponent >= 0
                transformed ^= 1 << exponent
        base[labels] = transformed
    local = {}
    for x, y, u, v in LABELS:
        matrix = [[0, 0], [0, 0]]
        for b in range(4):
            for i, j in product(range(2), repeat=2):
                matrix[i][j] ^= multiply(table[x, b][i][j], base[b, y, u, v])
        local[x, y, u, v] = tuple(tuple(row) for row in matrix)
    words, operator = word_operator(table, 2)
    sparse_rows = {word: [(other, value) for other, value in zip(words, row) if value]
                   for word, row in zip(words, operator)}
    checks = 0
    for first, second in product(LABELS, repeat=2):
        actual = auxiliary_product(local[first], local[second])[0][0]
        expected = 0
        for (a, b), coefficient in sparse_rows[first[0], second[0]]:
            tensor_value = multiply(base[(a,) + first[1:]], base[(b,) + second[1:]])
            expected ^= multiply(coefficient, tensor_value)
        assert actual == expected
        checks += 1
    return checks


def run():
    values = reflection_seed()
    variables, columns, kernel = row_lattice_system(values)
    generators = classify_fiber(values, variables, columns, kernel)
    word_entries = unit_checks = identities = 0
    for infinity in (False, True):
        table = noncommuting_leg_map(infinity)
        entries, checks = verify_word_units(table)
        word_entries += entries
        unit_checks += checks
        identities += verify_two_site_identity(values, table, infinity)
    return {
        "scope": "two specified row lattices; arbitrary auxiliary dimension by the triangular word proof",
        "seed": {"alphabet": 4, "degree_after_removing_even_ramification": 12,
                 "strict_endpoint_entries": 256},
        "zero_orders": EXPECTED_ZERO,
        "infinity_orders_with_common_shift_four": EXPECTED_INFINITY,
        "endpoint_preserving_laurent_variables": len(variables),
        "constraint_rank": len(variables) - len(kernel),
        "fiber_dimension": len(kernel),
        "complete_leg_generators_target_source_power": generators,
        "noncommuting_control": {
            "auxiliary_dimension": 2,
            "word_operator_entries_through_length_three": word_entries,
            "exact_polynomial_word_inverse_checks": unit_checks,
            "literal_two_site_identity_entries": identities,
        },
        "factoring_exponent_improved": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "evidence" / "triangular-auxiliary-row-fiber.json")
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
