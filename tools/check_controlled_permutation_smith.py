#!/usr/bin/env python3
"""Exact controls for the controlled-permutation Smith obstruction."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path


def binary_rank(rows, low=False):
    pivots = {}
    for row in rows:
        while row:
            column = ((row & -row).bit_length() - 1
                      if low else row.bit_length() - 1)
            if column in pivots:
                row ^= pivots[column]
            else:
                pivots[column] = row
                break
    return len(pivots)


def permutation_mask(p):
    q = len(p)
    return sum(1 << (q * p[v] + v) for v in range(q))


def permutation_power(p, exponent):
    answer = list(range(len(p)))
    for _ in range(exponent):
        answer = [p[a] for a in answer]
    return answer


def phi_rows(p, q):
    rows = [0] * (q * q)
    for source, target in enumerate(p):
        v, y = divmod(source, q)
        x, u = divmod(target, q)
        rows[q * x + y] ^= 1 << (q * u + v)
    return rows


def toeplitz_rank(coefficients, precision):
    n = len(coefficients[0])
    rows = []
    for i in range(precision):
        for row in range(n):
            value = 0
            for j in range(i + 1):
                if i - j < len(coefficients):
                    value ^= coefficients[i - j][row] << (n * j)
            rows.append(value)
    high = binary_rank(rows)
    assert high == binary_rank(rows, low=True)
    return high


def binomial_polynomial(exponent):
    """The polynomial (1+z)^exponent over F2."""
    answer, submask = 0, exponent
    while True:
        answer ^= 1 << submask
        if not submask:
            return answer
        submask = (submask - 1) & exponent


def check_table(table, full_polynomial=False):
    q, n = len(table), len(table) ** 2
    assert q >= 2 and not q & (q - 1)
    assert all(sorted(p) == list(range(q)) for p in table)
    inverse = []
    for p in table:
        inv = [0] * q
        for i, image in enumerate(p):
            inv[image] = i
        inverse.append(inv)

    U = [q * j + table[j][i] for i in range(q) for j in range(q)]
    assert sorted(U) == list(range(n))
    orbit, current = [], 0
    while current not in orbit:
        orbit.append(current)
        current = U[current]
    assert current == 0 and len(orbit) == n

    r = binary_rank([permutation_mask(p) for p in table])
    r1 = phi_rows(U, q)
    r2 = phi_rows(permutation_power(U, 2), q)
    r3 = phi_rows(permutation_power(U, 3), q)
    assert binary_rank(r1) == r == binary_rank(r1, low=True)
    expected_r1 = [0] * n
    expected_r2 = []
    fibers = Counter()
    for x in range(q):
        expected_r1[q * x + x] = permutation_mask(table[x])
        for y in range(q):
            pair = (table[x][y], inverse[y][x])
            expected_r2.append(1 << (q * pair[0] + pair[1]))
            fibers[pair] += 1
    assert r1 == expected_r1 and r2 == expected_r2
    image_size = len(fibers)
    assert binary_rank(r2) == image_size
    assert max(fibers.values()) <= r
    assert image_size * r >= n

    for a, b in fibers:
        possible_y = {inverse[x][a] for x in range(q)}
        evaluation_rows = [1 << y for y in possible_y]
        assert binary_rank(evaluation_rows) == len(possible_y) <= r
        recovered = {(table[y][b], y) for y in possible_y
                     if table[table[y][b]][y] == a}
        actual = {(x, y) for x in range(q) for y in range(q)
                  if (table[x][y], inverse[y][x]) == (a, b)}
        assert recovered == actual

    coefficients = [r1, [a ^ b for a, b in zip(r1, r2)],
                    [a ^ b for a, b in zip(r1, r3)]]
    t1 = toeplitz_rank(coefficients, 1)
    t2 = toeplitz_rank(coefficients, 2)
    t3 = toeplitz_rank(coefficients, 3)
    multiplicities = (t1, t2 - 2 * t1, t3 - 2 * t2 + t1)
    assert multiplicities[0] == r
    assert multiplicities[1] >= max(0, image_size - 2 * r)
    assert multiplicities[1] >= max(0, (n + r - 1) // r - 2 * r)

    position = {value: i for i, value in enumerate(orbit)}
    powers = [binomial_polynomial(i) for i in range(n)]
    row_polynomials = [[0] * n for _ in range(n)]
    column_polynomials = [[0] * n for _ in range(n)]
    endpoint_entries = 0
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    distance = (position[q * x + u] - position[q * v + y]) % n
                    exponent = 0 if not distance else n - distance
                    value = powers[exponent]
                    assert value.bit_count() % 2 == int(x == v and y == u)
                    row_polynomials[q * x + y][q * u + v] = value
                    column_polynomials[q * x + u][q * y + v] = value
                    for j in range(3):
                        assert (value >> (n - 1 - j) & 1) == (
                            coefficients[j][q * x + y] >> (q * u + v) & 1)
                    endpoint_entries += 1
    column_checks = []
    if full_polynomial:
        column_coefficients = [
            [sum(((column_polynomials[i][j] >> k) & 1) << j
                 for j in range(n)) for i in range(n)]
            for k in range(n)]
        for precision in sorted({1, 2, n, n + 1}):
            actual = toeplitz_rank(column_coefficients, precision)
            expected = precision + (n - 1) * max(0, precision - n)
            assert actual == expected
            column_checks.append({"precision": precision, "rank": actual})
    report = {
        "q": q, "distinct_controls": len(set(map(tuple, table))),
        "control_span_rank": r, "cycle_length": n,
        "second_reshuffling_rank": image_size,
        "maximum_fiber": max(fibers.values()),
        "row_infinity_order_multiplicities_0_1_2": list(multiplicities),
        "complete_row_profile_through_order_2": sum(multiplicities) == n,
        "endpoint_entries": endpoint_entries,
        "column_toeplitz_controls": column_checks,
    }
    return report


def two_control_table(q, control, T):
    shift = tuple((i + 1) % q for i in range(q))
    changed = tuple(T[shift[i]] for i in range(q))
    return [changed if control >> i & 1 else shift for i in range(q)]


def f4_multiply(a, b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return value


F4 = tuple(tuple(f4_multiply(a, b) for b in range(4)) for a in range(4))


def f4_binary_rows(matrix):
    rows = []
    for values in matrix:
        pair = [0, 0]
        for j, value in enumerate(values):
            for k, image in enumerate((value, F4[value][2])):
                for output in range(2):
                    pair[output] ^= ((image >> output) & 1) << (2 * j + k)
        rows.extend(pair)
    return rows


def f4_rank(matrix):
    binary = f4_binary_rows(matrix)
    result = binary_rank(binary)
    assert result % 2 == 0 and result == binary_rank(binary, low=True)
    return result // 2


def check_weighted_table(table):
    q, n = len(table), len(table) ** 2
    weights = []
    for j in range(q):
        exponents = [(1 + (j & 1)) * i + (j & 1) for i in range(q)]
        exponents[-1] -= sum(exponents)
        values = [(1, 2, 3)[e % 3] for e in exponents]
        product = 1
        for value in values:
            product = F4[product][value]
        assert product == 1
        weights.append(values)
    U = [q * j + table[j][i] for i in range(q) for j in range(q)]
    amplitudes = [weights[j][i] for i in range(q) for j in range(q)]
    orbit, current, product = [], 0, 1
    while current not in orbit:
        orbit.append(current)
        product = F4[product][amplitudes[current]]
        current = U[current]
    assert current == 0 and len(orbit) == n and product == 1
    matrices = []
    positions, values = list(range(n)), [1] * n
    for _ in range(3):
        values = [F4[value][amplitudes[position]]
                  for position, value in zip(positions, values)]
        positions = [U[position] for position in positions]
        matrix = [[0] * n for _ in range(n)]
        for source, (target, value) in enumerate(zip(positions, values)):
            v, y = divmod(source, q)
            x, u = divmod(target, q)
            matrix[q * x + y][q * u + v] = value
        matrices.append(matrix)
    leading, square, cube = matrices
    # Use the q^2-entry field vectors, then restrict scalars to F2.
    control_vectors = []
    for j in range(q):
        vector = [0] * n
        for v in range(q):
            vector[q * table[j][v] + v] = weights[j][v]
        control_vectors.append(vector)
    r = f4_rank(control_vectors)
    assert f4_rank(leading) == r
    image_size = f4_rank(square)
    assert image_size * r >= n
    coefficients = [leading,
                    [[a ^ b for a, b in zip(ar, br)]
                     for ar, br in zip(leading, square)],
                    [[a ^ b for a, b in zip(ar, br)]
                     for ar, br in zip(leading, cube)]]
    binary_coefficients = [f4_binary_rows(a) for a in coefficients]
    binary_ranks = [toeplitz_rank(binary_coefficients, p) for p in (1, 2, 3)]
    assert all(value % 2 == 0 for value in binary_ranks)
    ranks = [value // 2 for value in binary_ranks]
    multiplicities = [ranks[0], ranks[1] - 2 * ranks[0],
                      ranks[2] - 2 * ranks[1] + ranks[0]]
    assert multiplicities[1] >= max(0, (n + r - 1) // r - 2 * r)
    powers = [binomial_polynomial(i) for i in range(n)]
    checked = 0
    for source in range(n):
        target, amplitude = source, 1
        for distance in range(1, n + 1):
            amplitude = F4[amplitude][amplitudes[target]]
            target = U[target]
            polynomial = powers[n - distance]
            v, y = divmod(source, q)
            x, u = divmod(target, q)
            endpoint = amplitude if polynomial.bit_count() % 2 else 0
            assert endpoint == int(x == v and y == u)
            for j in range(3):
                coefficient = amplitude if polynomial >> (n - 1 - j) & 1 else 0
                assert coefficient == coefficients[j][q * x + y][q * u + v]
            checked += 1
    return {"q": q, "field": "F4=F2[w]/(w^2+w+1)",
            "cycle_length": n, "cycle_weight": product,
            "control_span_rank": r, "second_reshuffling_rank": image_size,
            "row_infinity_order_multiplicities_0_1_2": multiplicities,
            "literal_endpoint_and_three_jet_entries": checked}


def check():
    small_count, small_histogram = 0, Counter()
    q = 4
    for T in permutations(range(q)):
        for control in range(1, (1 << q) - 1):
            table = two_control_table(q, control, T)
            U = [q * j + table[j][i] for i in range(q) for j in range(q)]
            seen, current = set(), 0
            while current not in seen:
                seen.add(current)
                current = U[current]
            if len(seen) != q * q:
                continue
            record = check_table(table)
            small_count += 1
            small_histogram[record["row_infinity_order_multiplicities_0_1_2"][1]] += 1

    controls = []
    for q in (2, 4, 8, 16):
        T = list(range(q))
        T[0], T[1] = T[1], T[0]
        control = sum(1 << i for i in range(1, q // 2))
        record = check_table(two_control_table(q, control, T), q <= 4)
        expected = [1, 2, 1] if q == 2 else [2, q * q - 3, 1]
        assert record["row_infinity_order_multiplicities_0_1_2"] == expected
        controls.append(record)

    T = (1, 5, 4, 7, 0, 3, 2, 6)
    table = two_control_table(8, 173, T)
    noncommuting = check_table(table)
    assert noncommuting["row_infinity_order_multiplicities_0_1_2"][1] == 44
    first, second = table[0], table[1]
    assert any(first[second[i]] != second[first[i]] for i in range(8))

    shift = tuple((i + 1) % 8 for i in range(8))
    permutations_four = []
    for bits in range(4):
        T = list(range(8))
        if bits & 1:
            T[0], T[1] = T[1], T[0]
        if bits & 2:
            T[2], T[3] = T[3], T[2]
        permutations_four.append(tuple(T[shift[i]] for i in range(8)))
    word = (0, 0, 0, 0, 1, 2, 1, 3)
    dependent = check_table([permutations_four[i] for i in word])
    assert dependent["distinct_controls"] == 4
    assert dependent["control_span_rank"] == 3
    weighted = [
        check_weighted_table(two_control_table(4, 2, (1, 0, 2, 3))),
        check_weighted_table(two_control_table(8, 173, (1, 5, 4, 7, 0, 3, 2, 6))),
    ]

    potential_checks = 0
    for q in (2, 4, 8, 16):
        n = q * q
        for s in map(Fraction, ("1", "5/4", "3/2", "2", "3")):
            a0 = s ** (n - 1) + (n - 1) / s
            a1 = s ** (n - 2) + (n - 1) / (s * s)
            assert a0 >= n and a1 >= n - 1
            potential_checks += 2
        for r in range(1, n + 1):
            lower = n * r + (n - 1) * max(0, Fraction(n, r) - 2 * r)
            assert 2 * lower * lower >= n ** 3
            potential_checks += 1
    assert 2 * 17 ** 2 > 24 ** 2
    return {
        "status": "pass", "field": "F2 and F4; proof over every characteristic-two field",
        "small_exhaustive_scope": "q=4, cyclic L0, every permutation T and nontrivial binary control subset",
        "small_regular_cases": small_count,
        "small_order_one_histogram": dict(sorted(small_histogram.items())),
        "controlled_transposition_records": controls,
        "noncommuting_two_control_record": noncommuting,
        "dependent_four_control_record": dependent,
        "weighted_monomial_F4_records": weighted,
        "rational_potential_checks": potential_checks,
        "proved_basic_Smith_floor": "B(s) >= q^3/sqrt(2); rho >= 17/24 for q>=8",
        "scope": "Single-cycle controlled monomial nilpotent paths with total cycle weight one, and the bare two-spectrum Smith certificate. No factoring lower bound or improved factoring exponent.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = check()
    report["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
