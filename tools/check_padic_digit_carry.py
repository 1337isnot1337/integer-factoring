#!/usr/bin/env python3
"""Exact finite checks of the base-prime/carry-lift research note.

The reference implementations use ordinary Python lists. They check
identities and algorithms, not the packed-workspace machine itself.
"""

from __future__ import annotations

import argparse
from itertools import product
import json
from math import factorial
from pathlib import Path
import random

from check_hankel_digit_box import det_gauss, hankel


def digits(value, base, length):
    answer = []
    for _ in range(length):
        answer.append(value % base)
        value //= base
    assert value == 0
    return answer


def convolution(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return result


def rank_pivots(matrix, p):
    if not matrix:
        return [], []
    a = [list(row) for row in matrix]
    original_rows = list(range(len(a)))
    chosen_rows, chosen_columns = [], []
    row = 0
    for col in range(len(a[0])):
        selected = next((j for j in range(row, len(a)) if a[j][col] % p), None)
        if selected is None:
            continue
        a[row], a[selected] = a[selected], a[row]
        original_rows[row], original_rows[selected] = original_rows[selected], original_rows[row]
        chosen_rows.append(original_rows[row])
        chosen_columns.append(col)
        inv = pow(a[row][col] % p, -1, p)
        for j in range(row + 1, len(a)):
            multiplier = a[j][col] * inv % p
            for k in range(col, len(a[0])):
                a[j][k] = (a[j][k] - multiplier * a[row][k]) % p
        row += 1
        if row == len(a):
            break
    return chosen_rows, chosen_columns


def digit_lift_solutions(matrix, rhs, p, precision):
    size = len(matrix[0])
    ii, jj = rank_pivots(matrix, p)
    free = [j for j in range(size) if j not in jj]
    minor = [[matrix[i][j] % p for j in jj] for i in ii]
    denominator = det_gauss(minor, p)
    modulus = p ** precision
    result = []
    candidates = 0
    for values in product(range(p), repeat=len(free)):
        candidates += 1
        x = [0] * size
        for j, value in zip(free, values):
            x[j] = value
        target = [(rhs[i] - sum(matrix[i][j] * x[j] for j in free)) % p for i in ii]
        for k, col in enumerate(jj):
            replaced = [[target[a] if b == k else minor[a][b] for b in range(len(jj))] for a in range(len(ii))]
            x[col] = det_gauss(replaced, p) * pow(denominator, -1, p) % p
        if all((sum(a * b for a, b in zip(row, x)) - y) % modulus == 0 for row, y in zip(matrix, rhs)):
            result.append(tuple(x))
    return sorted(result), len(jj), candidates


class PrimePowerRoot:
    def __init__(self, base, precision, coefficient_prime):
        self.base = base
        self.precision = precision
        self.modulus = base ** precision
        self.step = self.modulus // base
        self.degree = self.modulus - self.step
        self.field = coefficient_prime
        self.zero = (0,) * self.degree

    def phase(self, exponent):
        r = exponent % self.modulus
        if r < self.degree:
            return tuple(int(i == r) for i in range(self.degree))
        offset = r - self.degree
        return tuple(self.field - 1 if i % self.step == offset else 0 for i in range(self.degree))

    def add(self, a, b):
        return tuple((x + y) % self.field for x, y in zip(a, b))

    def phase_sum(self, exponents):
        answer = self.zero
        for exponent in exponents:
            answer = self.add(answer, self.phase(exponent))
        return answer

    def defining_polynomial(self):
        return [int(i % self.step == 0) for i in range(self.degree + 1)]


def trim(poly):
    poly = list(poly)
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def poly_mod(a, b, p):
    a = trim(a)
    b = trim(b)
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b):
        factor = a[-1] * inv % p
        shift = len(a) - len(b)
        for i, value in enumerate(b):
            a[shift + i] = (a[shift + i] - factor * value) % p
        a = trim(a)
    return a


def is_unit(poly, defining, p):
    a, b = trim(poly), trim(defining)
    while b:
        a, b = b, poly_mod(a, b, p)
    return len(a) == 1


def partitions(n, minimum=1):
    if n == 0:
        yield ()
    else:
        for first in range(minimum, n + 1):
            for rest in partitions(n - first, first):
                yield (first,) + rest


def pfree(value, p):
    valuation = 0
    while value % p == 0:
        valuation += 1
        value //= p
    return valuation, value


def factorial_unit(n, p, modulus):
    valuation, unit = 0, 1
    for value in range(1, n + 1):
        v, u = pfree(value, p)
        valuation += v
        unit = unit * u % modulus
    return valuation, unit


def matrix_product(a, b, modulus=None):
    n = len(a)
    c = [[sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    if modulus:
        c = [[x % modulus for x in row] for row in c]
    return c


def trace_power(matrix, k, modulus):
    n = len(matrix)
    result = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(k):
        result = matrix_product(result, matrix, modulus)
    return sum(result[i][i] for i in range(n)) % modulus


def determinant_partition_lift(matrix, p):
    n = len(matrix)
    if n == 0:
        return 1
    valuation = sum(pfree(i, p)[0] for i in range(1, n + 1))
    modulus = p ** (valuation + 1)
    v, unit = factorial_unit(n, p, modulus)
    assert v == valuation
    total = 0
    for parts in partitions(n):
        counts = {k: parts.count(k) for k in set(parts)}
        zvaluation, zunit = 0, 1
        z_exact = 1
        for k, count in counts.items():
            kv, ku = pfree(k, p)
            fv, fu = factorial_unit(count, p, modulus)
            zvaluation += count * kv + fv
            zunit = zunit * pow(ku, count, modulus) * fu % modulus
            z_exact *= k ** count * factorial(count)
        assert factorial(n) % z_exact == 0
        assert zvaluation <= valuation
        coefficient = pow(p, valuation - zvaluation, modulus) * unit * pow(zunit, -1, modulus) % modulus
        assert coefficient == factorial(n) // z_exact % modulus
        term = coefficient
        for k, count in counts.items():
            term = term * pow(trace_power(matrix, k, modulus), count, modulus) % modulus
        total = (total + (-1 if (n - len(parts)) % 2 else 1) * term) % modulus
    assert total % (p ** valuation) == 0
    return (total // (p ** valuation)) * pow(unit % p, -1, p) % p


def sum_bit(term_bit, count, bit):
    carry = 0
    answer = 0
    for position in range(bit + 1):
        column = carry + sum(term_bit(i, position) for i in range(count))
        answer = column & 1
        carry = column >> 1
    return answer


def multiply_bit(left_bit, right_bit, bit):
    carry = 0
    answer = 0
    for position in range(bit + 1):
        column = carry + sum(left_bit(i) * right_bit(position - i) for i in range(position + 1))
        answer = column & 1
        carry = column >> 1
    return answer


def integer_power_bit(matrix, k, i, j, bit):
    if k == 0:
        return int(i == j and bit == 0)
    if k == 1:
        return matrix[i][j] >> bit & 1
    left = k // 2
    return sum_bit(
        lambda pivot, position: multiply_bit(
            lambda b: integer_power_bit(matrix, left, i, pivot, b),
            lambda b: integer_power_bit(matrix, k - left, pivot, j, b),
            position,
        ),
        len(matrix),
        bit,
    )


def run_checks():
    rng = random.Random(20261031)
    reduced_collisions = []
    for p, m, value, good, bad in (
        (2, 3, 5, (1, 5), (3, 3)),
        (3, 3, 133, (7, 19), (5, 11)),
        (5, 2, 143, (11, 13), (19, 17)),
    ):
        first = convolution(digits(good[0], p, m), digits(good[1], p, m))
        second = convolution(digits(bad[0], p, m), digits(bad[1], p, m))
        assert [c % p for c in first] == [c % p for c in second]
        assert good[0] * good[1] == value and bad[0] * bad[1] != value
        assert all(0 < x <= value for x in good + bad)
        reduced_collisions.append({"p": p, "N": value, "true_pair": good, "false_pair": bad, "common_reduced_convolution": [c % p for c in first]})
    relaxed_cases, relaxed_examples = 0, []
    for p, m in ((2, 2), (2, 3), (3, 2), (5, 2)):
        for value in range(1, p ** m):
            if value % p == 0:
                continue
            target = digits(value, p, 2 * m - 1)
            count = 0
            exact = 0
            first_false = None
            for a, b in product(range(p ** m), repeat=2):
                c = convolution(digits(a, p, m), digits(b, p, m))
                carry = [0] + [(target[j] - c[j]) % p for j in range(1, len(c))] + [0]
                valid = all((c[j] + carry[j] - target[j] - p * carry[j + 1]) % p == 0 for j in range(len(c)))
                if valid:
                    count += 1
                    exact += int(a * b == value)
                    if a * b != value and first_false is None:
                        first_false = {"a": a, "b": b, "convolution": c, "relaxed_carries": carry}
            assert count == (p - 1) * p ** (2 * m - 2)
            relaxed_examples.append({"p": p, "m": m, "N": value, "relaxed": count, "actual_factor_pairs": exact, "first_false": first_false})
            relaxed_cases += 1

    nowrap_checks = 0
    for p, m in ((2, 1), (2, 2), (3, 3), (5, 2)):
        cmax, tmax = m * (p - 1) ** 2, m * (p - 1)
        modulus = p
        while modulus <= (m * p + 1) * (p - 1):
            modulus *= p
        for c, t, u, digit in product(range(cmax + 1), range(tmax + 1), range(tmax + 1), range(p)):
            discrepancy = c + t - digit - p * u
            assert abs(discrepancy) < modulus
            assert (discrepancy % modulus == 0) == (discrepancy == 0)
            nowrap_checks += 1

    lift_cases = 0
    for p, precision in ((2, 2), (2, 3), (3, 2), (5, 2)):
        modulus = p ** precision
        for rows, columns in ((1, 2), (2, 2), (2, 3), (3, 2)):
            for case in range(20):
                matrix = [[0 if case == 0 else rng.randrange(modulus) for _ in range(columns)] for _ in range(rows)]
                rhs = [rng.randrange(modulus) for _ in range(rows)]
                if case % 3 == 0:
                    witness = [rng.randrange(p) for _ in range(columns)]
                    rhs = [sum(a * b for a, b in zip(row, witness)) % modulus for row in matrix]
                expected = [x for x in product(range(p), repeat=columns) if all((sum(a * b for a, b in zip(row, x)) - y) % modulus == 0 for row, y in zip(matrix, rhs))]
                result, rank, candidates = digit_lift_solutions(matrix, rhs, p, precision)
                assert result == expected
                assert candidates == p ** (columns - rank)
                lift_cases += 1

    support_cases = []
    for p, precision, coefficient_prime in ((2, 1, 3), (2, 2, 3), (2, 3, 3), (3, 2, 2), (5, 2, 2)):
        ring = PrimePowerRoot(p, precision, coefficient_prime)
        support = 0
        units = 0
        for frequency in range(ring.modulus):
            value = ring.phase_sum(frequency * digit for digit in range(p))
            expected_zero = frequency != 0 and frequency % ring.step == 0
            assert (value == ring.zero) == expected_zero
            if not expected_zero:
                support += 1
                assert is_unit(value, ring.defining_polynomial(), coefficient_prime)
                units += 1
        assert support == ring.modulus - p + 1
        support_cases.append({"p": p, "precision": precision, "coefficient_prime": coefficient_prime, "nonzero_frequencies": support, "unit_frequencies": units})

    carry_cases = []
    for p, precision, m, coefficient_prime in ((2, 2, 2, 3), (2, 3, 3, 3), (3, 2, 2, 2)):
        ring = PrimePowerRoot(p, precision, coefficient_prime)
        for value in range(1, p ** m):
            legal = [c for c in product(range(ring.modulus), repeat=m) if sum(c[i] * p ** i for i in range(m)) == value]
            selected = None
            for initial in product(range(ring.modulus), repeat=m):
                if initial[-1] % p == 0:
                    continue
                transformed = ring.phase_sum(sum(x * y for x, y in zip(initial, c)) for c in legal)
                if transformed != ring.zero:
                    selected = initial
                    break
            assert selected is not None
            extensions = 0
            for tail in product(range(0, ring.modulus, p), repeat=m - 1):
                assert det_gauss(hankel(selected + tail, m), p) != 0
                extensions += 1
            assert extensions == (ring.modulus // p) ** (m - 1)
            carry_cases.append({"p": p, "precision": precision, "m": m, "N": value, "selected_nonzero_lower_frequency": selected, "invertible_nonzero_extensions": extensions})

    determinant_cases = 0
    for p in (2, 3, 5, 7):
        for n in range(7):
            for _ in range(20):
                matrix = [[rng.randrange(p) for _ in range(n)] for _ in range(n)]
                assert determinant_partition_lift(matrix, p) == det_gauss(matrix, p)
                determinant_cases += 1
    bit_query_cases = 0
    for matrix in ([[1, 2], [3, 1]], [[0, 1], [1, 1]]):
        power = [[int(i == j) for j in range(2)] for i in range(2)]
        for k in range(1, 4):
            power = matrix_product(power, matrix)
            for i, j in product(range(2), repeat=2):
                for bit in range(power[i][j].bit_length() + 1):
                    assert integer_power_bit(matrix, k, i, j, bit) == (power[i][j] >> bit & 1)
                    bit_query_cases += 1

    return {"status": "passed", "scope": "Finite exact identities; no factoring exponent improvement and no measured bit-space claim.", "reduced_convolution_collisions": reduced_collisions, "relaxed_carry_cases": relaxed_cases, "relaxed_carry_examples": relaxed_examples, "no_wrap_checks": nowrap_checks, "digit_lift_cases": lift_cases, "digit_fourier_support": support_cases, "nonzero_carry_hankel_cases": carry_cases, "small_characteristic_determinant_cases": determinant_cases, "virtual_integer_power_bit_cases": bit_query_cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_checks()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
