#!/usr/bin/env python3
"""Exact polynomial crossing interpolation in prescribed local row lattices.

The finite search is complete for its specified gauges and row lattices.
It is not a general obstruction to improving the factoring exponent.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_block_smith_crossings import binary_rank, series_inverse, smith_profile, toeplitz_rank
from check_schur_gaussian_reduction import polynomial_multiply, seed_data


ROOT = Path(__file__).resolve().parent.parent
ZERO_HEIGHTS = ((0, 0, 0, 0), (0, 0, 0, 0), (0, 1, 1, 2), (0, -1, -1, -2))
INFINITY_HEIGHTS = ((0, 1, 1, 2), (0, 1, 1, 2), (0, 0, 0, 0), (0, 2, 2, 4))
EXPECTED = (0, 1, 1, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 7, 7, 8)
LABELS = tuple(tuple((mask >> (2 * i)) & 3 for i in range(4)) for mask in range(256))


def polynomial_matrix_product(left, right, precision):
    mask = (1 << precision) - 1
    output = [[0] * len(right[0]) for _ in left]
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if value:
                for j, other in enumerate(right[k]):
                    if other:
                        output[i][j] ^= polynomial_multiply(value, other) & mask
    return output


def smith_jet_changes(matrix, precision):
    """Return explicit polynomial unit matrices and visible diagonal orders."""
    n = len(matrix)
    mask = (1 << precision) - 1
    work = [[value & mask for value in row] for row in matrix]
    left = [[int(i == j) for j in range(n)] for i in range(n)]
    right = [[int(i == j) for j in range(n)] for i in range(n)]
    orders = []
    for k in range(n):
        order, pi, pj = min(((value & -value).bit_length() - 1, i, j)
                            for i in range(k, n) for j in range(k, n)
                            if (value := work[i][j]))
        orders.append(order)
        work[k], work[pi] = work[pi], work[k]
        left[k], left[pi] = left[pi], left[k]
        if pj != k:
            for row in work:
                row[k], row[pj] = row[pj], row[k]
            for row in right:
                row[k], row[pj] = row[pj], row[k]
        inverse = series_inverse(work[k][k] >> order, precision)
        work[k] = [polynomial_multiply(value, inverse) & mask for value in work[k]]
        left[k] = [polynomial_multiply(value, inverse) & mask for value in left[k]]
        for i in range(k + 1, n):
            factor = work[i][k] >> order
            if factor:
                work[i] = [a ^ (polynomial_multiply(factor, b) & mask)
                           for a, b in zip(work[i], work[k])]
                left[i] = [a ^ (polynomial_multiply(factor, b) & mask)
                           for a, b in zip(left[i], left[k])]
        for j in range(k + 1, n):
            factor = work[k][j] >> order
            if factor:
                for i in range(n):
                    work[i][j] ^= polynomial_multiply(factor, work[i][k]) & mask
                    right[i][j] ^= polynomial_multiply(factor, right[i][k]) & mask
    expected = [[(1 << orders[i]) if i == j else 0 for j in range(n)] for i in range(n)]
    assert work == expected
    assert polynomial_matrix_product(polynomial_matrix_product(left, matrix, precision), right, precision) == expected
    for change in (left, right):
        assert binary_rank([sum((value & 1) << j for j, value in enumerate(row)) for row in change]) == n
    return tuple(orders), left, right


def flatten_position(labels, cut):
    x, y, u, v = labels
    return (4 * x + y, 4 * u + v) if cut == 0 else (4 * x + u, 4 * y + v)


def height(labels, heights):
    return sum(row[label] for row, label in zip(heights, labels))


def gauged_matrix(values, heights, cut, infinity):
    output = [[0] * 16 for _ in range(16)]
    for labels, polynomial in zip(LABELS, values):
        i, j = flatten_position(labels, cut)
        h = height(labels, heights)
        while polynomial:
            bit = polynomial & -polynomial
            degree = bit.bit_length() - 1
            exponent = h - degree if infinity else h + degree
            assert exponent >= 0
            output[i][j] ^= 1 << exponent
            polynomial ^= bit
    return output


def encode_rows(matrix, shift, targets):
    result = offset = 0
    for row, width in zip(matrix, targets):
        mask = (1 << width) - 1
        for value in row:
            result |= ((value << shift) & mask) << offset
            offset += width
    return result, offset


def column_elimination(columns):
    basis, kernel = {}, []
    for index, column in enumerate(columns):
        combination = 1 << index
        while column:
            lead = column.bit_length() - 1
            if lead in basis:
                old, witness = basis[lead]
                column ^= old
                combination ^= witness
            else:
                basis[lead] = column, combination
                break
        if not column:
            kernel.append(combination)
    return basis, kernel


def solve_in_basis(basis, target):
    witness = 0
    while target:
        lead = target.bit_length() - 1
        if lead not in basis:
            return None
        column, combination = basis[lead]
        target ^= column
        witness ^= combination
    return witness


def interpolation_columns(changes, targets, even=True):
    """All Laurent perturbations vanishing at one, with no degree cutoff."""
    variables, columns = [], []
    first_width = 16 * sum(targets[0])
    for mask, labels in enumerate(LABELS):
        if even and mask.bit_count() % 2:
            continue
        h0, hi = height(labels, ZERO_HEIGHTS), height(labels, INFINITY_HEIGHTS)
        blocks = []
        for cut, (left, right) in enumerate(changes):
            row, column = flatten_position(labels, cut)
            matrix = [[polynomial_multiply(left[i][row], right[column][j])
                       for j in range(16)] for i in range(16)]
            blocks.append([[value ^ (value << 1) for value in row] for row in matrix])
        for k in range(-h0, hi):
            first, _ = encode_rows(blocks[0], k + h0, targets[0])
            second, _ = encode_rows(blocks[1], hi - k - 1, targets[1])
            variables.append((mask, k))
            columns.append(first | (second << first_width))
    return variables, columns


def baseline_target(orders, targets):
    result = offset = 0
    for sigma, tau in zip(orders, targets):
        for i, width in enumerate(tau):
            for j in range(16):
                if i == j and sigma[i] < width:
                    result ^= 1 << (offset + sigma[i])
                offset += width
    return result


def apply_perturbation(base, variables, witness):
    # Return a Laurent table as exponent -> coefficient dictionaries.
    values = [{i: 1 for i in range(value.bit_length()) if value >> i & 1} for value in base]
    for index, (mask, k) in enumerate(variables):
        if witness >> index & 1:
            for exponent in (k, k + 1):
                if exponent in values[mask]:
                    del values[mask][exponent]
                else:
                    values[mask][exponent] = 1
    common_shift = max(0, -min((min(value) for value in values if value), default=0))
    return [sum(1 << (exponent + common_shift) for exponent in value) for value in values], common_shift


def wick_failures(values, common_shift=0):
    cache = {0: 1}

    def pfaffian(mask):
        if mask not in cache:
            first = mask & -mask
            remaining = mask ^ first
            value = 0
            for i in range(8):
                second = 1 << i
                if remaining & second:
                    value ^= polynomial_multiply(values[first | second], pfaffian(remaining ^ second))
            cache[mask] = value
        return cache[mask]

    assert values[0] == 1 << common_shift
    return [mask for mask in range(1, 256) if mask.bit_count() % 2 == 0
            and (values[mask] << (common_shift * (mask.bit_count() // 2 - 1))) != pfaffian(mask)]


def gauged_leg_operator(matrix, heights, infinity):
    """Gauge a physical leg map, with its row the target and column the source."""
    result = [[0] * 4 for _ in range(4)]
    for target, source in product(range(4), repeat=2):
        polynomial = matrix[target][source]
        while polynomial:
            bit = polynomial & -polynomial
            degree = bit.bit_length() - 1
            exponent = heights[target] - heights[source] + (-degree if infinity else degree)
            assert exponent >= 0
            result[target][source] ^= 1 << exponent
            polynomial ^= bit
    return result


def check_unrestricted_fiber(base, changes, orders):
    """Classify every Laurent perturbation, including odd physical parity."""
    variables, columns = interpolation_columns(changes, orders, even=False)
    basis, kernel = column_elimination(columns)
    assert len(variables) == 1024 and len(basis) == 1018 and len(kernel) == 6
    lookup = {variable: i for i, variable in enumerate(variables)}
    generators = ((1, 0, 1), (2, 0, 1), (3, 1, 1), (3, 2, 1),
                  (3, 0, 2), (3, 0, 3))
    known = []
    for target, source, shift in generators:
        witness = 0
        for mask in range(256):
            if mask >> 6 == target:
                value = base[(mask & 63) | (source << 6)]
                for exponent in range(value.bit_length()):
                    if value >> exponent & 1:
                        witness ^= 1 << lookup[mask, shift + exponent]
        image = 0
        for i, column in enumerate(columns):
            if witness >> i & 1:
                image ^= column
        assert image == 0
        known.append(witness)
    assert binary_rank(known) == 6
    identity = [[int(i == j) for j in range(4)] for i in range(4)]
    records, toeplitz_checks = [], 0
    for selector in range(64):
        witness = 0
        nilpotent = [[0] * 4 for _ in range(4)]
        for i, (target, source, shift) in enumerate(generators):
            if selector >> i & 1:
                witness ^= known[i]
                nilpotent[target][source] ^= 3 << shift
        square = polynomial_matrix_product(nilpotent, nilpotent, 20)
        assert polynomial_matrix_product(square, nilpotent, 20) == [[0] * 4 for _ in range(4)]
        leg_map = [[identity[i][j] ^ nilpotent[i][j] for j in range(4)] for i in range(4)]
        inverse = [[leg_map[i][j] ^ square[i][j] for j in range(4)] for i in range(4)]
        assert polynomial_matrix_product(leg_map, inverse, 20) == identity
        for heights, infinity in ((ZERO_HEIGHTS[3], False), (INFINITY_HEIGHTS[3], True)):
            gauged = gauged_leg_operator(leg_map, heights, infinity)
            gauged_inverse = gauged_leg_operator(inverse, heights, infinity)
            assert polynomial_matrix_product(gauged, gauged_inverse, 20) == identity
            assert binary_rank([sum((value & 1) << j for j, value in enumerate(row))
                                for row in gauged]) == 4
        values, common_shift = apply_perturbation(base, variables, witness)
        assert common_shift == 0
        for mask, (labels, value) in enumerate(zip(LABELS, values)):
            reconstructed = 0
            for source in range(4):
                reconstructed ^= polynomial_multiply(
                    leg_map[labels[3]][source], base[(mask & 63) | (source << 6)])
            assert value == reconstructed
            assert value.bit_count() % 2 == int(labels[0] == labels[3] and labels[1] == labels[2])
        profiles = []
        for heights, cut, infinity in ((ZERO_HEIGHTS, 0, False), (INFINITY_HEIGHTS, 1, True)):
            matrix = gauged_matrix(values, heights, cut, infinity)
            degree = max(value.bit_length() - 1 for row in matrix for value in row)
            actual, _ = smith_profile(matrix, degree)
            assert actual == EXPECTED
            if selector == 1:  # independent Toeplitz audit of an odd-parity member
                for precision in range(1, max(actual) + 2):
                    assert toeplitz_rank(matrix, precision) == sum(
                        max(0, precision - order) for order in actual)
                    toeplitz_checks += 1
            profiles.append(actual)
        records.append({"selector": selector, "zero_orders": profiles[0],
                        "infinity_orders": profiles[1], "endpoint_entries": 256,
                        "regular_leg_map_and_inverse_at_both_ends": True})
    return {
        "all_laurent_endpoint_variables": len(variables),
        "unrestricted_constraint_rank": len(basis),
        "unrestricted_fiber_dimension": len(kernel),
        "unrestricted_fiber_size": 64,
        "local_v_generators": [{"target": target, "source": source,
                                "coefficient": f"(1+z)*z^{shift}"}
                               for target, source, shift in generators],
        "all_parity_fiber_checks": records,
        "additional_independent_toeplitz_ranks": toeplitz_checks,
    }


def check():
    base = [0] * 256
    for entry in seed_data()["local"]:
        mask = sum(label << (2 * i) for i, label in enumerate(entry["indices"]))
        base[mask] = entry["matrix"][0][0]
    changes, orders, base_matrices = [], [], []
    for heights, cut, infinity in ((ZERO_HEIGHTS, 0, False), (INFINITY_HEIGHTS, 1, True)):
        matrix = gauged_matrix(base, heights, cut, infinity)
        sigma, left, right = smith_jet_changes(matrix, 20)
        assert sigma == EXPECTED
        changes.append((left, right))
        orders.append(sigma)
        base_matrices.append(matrix)
    variables, columns = interpolation_columns(changes, orders)
    basis, kernel = column_elimination(columns)
    assert len(variables) == 512 and len(basis) == 510 and len(kernel) == 2
    lookup = {variable: i for i, variable in enumerate(variables)}
    known = []
    for k in (2, 3):
        witness = 0
        for mask in range(256):
            if mask & 192 == 192:
                value = base[mask ^ 192]
                for exponent in range(value.bit_length()):
                    if value >> exponent & 1:
                        witness ^= 1 << lookup[(mask, k + exponent)]
        assert witness
        image = 0
        for i, column in enumerate(columns):
            if witness >> i & 1:
                image ^= column
        assert image == 0
        known.append(witness)
    assert binary_rank(known) == 2
    records, toeplitz_checks = [], 0
    for a, b in product(range(2), repeat=2):
        values, common_shift = apply_perturbation(base, variables, (known[0] if a else 0) ^ (known[1] if b else 0))
        assert common_shift == 0
        assert not wick_failures(values)
        for labels, value in zip(LABELS, values):
            assert value.bit_count() % 2 == int(labels[0] == labels[3] and labels[1] == labels[2])
        profiles = []
        for heights, cut, infinity in ((ZERO_HEIGHTS, 0, False), (INFINITY_HEIGHTS, 1, True)):
            matrix = gauged_matrix(values, heights, cut, infinity)
            degree = max(value.bit_length() - 1 for row in matrix for value in row)
            actual, _ = smith_profile(matrix, degree)
            assert actual == EXPECTED
            for precision in range(1, max(actual) + 2):
                assert toeplitz_rank(matrix, precision) == sum(max(0, precision - order) for order in actual)
                toeplitz_checks += 1
            profiles.append(actual)
        records.append({"a": a, "b": b, "zero_orders": profiles[0], "infinity_orders": profiles[1],
                        "endpoint_entries": 256, "wick_identities": 128})

    # A non-Gaussian endpoint-preserving perturbation is excluded by the flags.
    quartic = base[:]
    for mask in range(256):
        if mask & 15 == 15:
            value = base[mask ^ 15]
            quartic[mask] ^= value ^ (value << 1)
    assert 15 in wick_failures(quartic)
    assert all(value.bit_count() % 2 == int(labels[0] == labels[3] and labels[1] == labels[2])
               for value, labels in zip(quartic, LABELS))
    violated_flags = 0
    for side, (heights, cut, infinity) in enumerate(((ZERO_HEIGHTS, 0, False), (INFINITY_HEIGHTS, 1, True))):
        matrix = gauged_matrix(quartic, heights, cut, infinity)
        changed = polynomial_matrix_product(polynomial_matrix_product(changes[side][0], matrix, 20), changes[side][1], 20)
        violated_flags += bool(encode_rows(changed, 0, orders[side])[0])
    assert violated_flags

    # A relaxed row lattice produces a genuine non-Gaussian crossing.
    relaxed = [[0] * 16, [0] * 16]
    relaxed[0][3] = 3
    rvariables, rcolumns = interpolation_columns(changes, relaxed)
    rbasis, _ = column_elimination(rcolumns)
    witness = solve_in_basis(rbasis, baseline_target(orders, relaxed))
    assert witness is not None
    values, common_shift = apply_perturbation(base, rvariables, witness)
    failures = wick_failures(values, common_shift)
    assert failures
    assert all(value.bit_count() % 2 == int(labels[0] == labels[3] and labels[1] == labels[2])
               for value, labels in zip(values, LABELS))
    assert common_shift == 0
    for side, (heights, cut, infinity) in enumerate(((ZERO_HEIGHTS, 0, False), (INFINITY_HEIGHTS, 1, True))):
        matrix = gauged_matrix(values, heights, cut, infinity)
        changed = polynomial_matrix_product(polynomial_matrix_product(changes[side][0], matrix, 20), changes[side][1], 20)
        assert encode_rows(changed, 0, relaxed[side])[0] == 0

    unrestricted = check_unrestricted_fiber(base, changes, orders)
    return {
        "status": "exact complete interpolation fiber for specified local row lattices; no improved exponent",
        "zero_heights": ZERO_HEIGHTS, "infinity_heights": INFINITY_HEIGHTS,
        "zero_cut": 0, "infinity_cut": 1, "row_targets": EXPECTED,
        "even_laurent_endpoint_variables": len(variables), "constraint_rank": len(basis),
        "fiber_dimension": len(kernel), "fiber_formula": "G * (1 + (a*z^2+b*z^3)*(1+z)*theta_v_plus*theta_v_minus)",
        "fiber_checks": records, "independent_toeplitz_ranks": toeplitz_checks,
        "unrestricted_parity_classification": unrestricted,
        "quartic_wick_negative_control": 15, "quartic_violated_row_lattices": violated_flags,
        "relaxed_non_gaussian_control": {
            "row_targets": relaxed, "common_monomial_shift": common_shift,
            "perturbation_terms": [variable for i, variable in enumerate(rvariables) if witness >> i & 1],
            "violated_wick_masks": failures,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = check()
    report["source_sha256"] = {name: sha256((ROOT / name).read_bytes()).hexdigest()
                               for name in ("tools/check_polynomial_crossing_interpolation.py",
                                            "tools/check_block_smith_crossings.py",
                                            "tools/check_schur_gaussian_reduction.py")}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
