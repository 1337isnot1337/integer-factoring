#!/usr/bin/env python3
"""Exact constant-field, constructor-budget, and rank-two weight controls."""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import random

from check_carry_fourier_elimination import BinaryField, construct_field_and_root, field_rank, geometric
from check_frobenius_sparse_recovery import (
    TensorAlgebra, binary_coefficient_evaluate, construct_field, degree_family,
)
from check_weighted_dyadic_carry import digits, genuine_carries


def adjusted_width(preliminary):
    return preliminary + (2 - preliminary) % 6


def cube_root(field):
    for candidate in range(2, field.size):
        beta = field.power(candidate, (field.size - 1) // 3)
        if beta != 1:
            assert field.power(beta, 3) == 1
            assert field.mul(beta, beta) ^ beta ^ 1 == 0
            return beta, candidate
    raise AssertionError("missing nontrivial cube root")


def check_field_interfaces(counts):
    rng = random.Random(26022026)
    records = []
    for d, width in ((2, 2), (4, 2), (8, 2), (2, 8)):
        t = d.bit_length() - 1
        m = d * ((2 * width + t + d - 1) // d)
        modulus, base = (1 << m) + 1, 1 << width
        field, root, construction = construct_field_and_root(m)
        beta, candidate = cube_root(field)
        gamma = field.power(beta, 1 << (t % 2))
        assert m % 2 == 0 and modulus % 3 == 2
        assert base % 9 == 4 and field.power(beta, base) == beta
        assert field.power(gamma, d) == beta
        assert 1 ^ field.power(beta, modulus) == beta
        frequencies = sorted({0, 1, 2, modulus - 1,
                              *(rng.randrange(modulus) for _ in range(6))})
        for scalar, length in ((beta, base), (gamma, d * base)):
            for frequency in frequencies:
                y = field.mul(scalar, field.power(root, frequency))
                denominator_inverse = field.mul(geometric(field, y, modulus),
                                                 field.inverse(1 ^ field.power(scalar, modulus)))
                assert field.mul(1 ^ y, denominator_inverse) == 1
                value = geometric(field, y, length)
                inverse = field.mul(1 ^ y, field.power(denominator_inverse, length))
                assert field.mul(value, inverse) == 1
                counts["fixed_field_box_unit_inverses"] += 1
        for ell in (3, 5, 7):
            recovery = construct_field(ell)
            algebra = TensorAlgebra(field, recovery)
            alpha = 2
            for position in range(4):
                positional = recovery.power(alpha, base ** position)
                assert positional not in (0, 1)
                for frequency in (0, 1, modulus - 1):
                    y = algebra.mul(algebra.weight_scalar(positional),
                                    algebra.root_scalar(field.power(root, frequency)))
                    inverse = algebra.mul(
                        algebra.geometric(y, modulus),
                        algebra.weight_scalar(recovery.inverse(1 ^ recovery.power(positional, modulus))))
                    assert algebra.mul(algebra.add(algebra.one, y), inverse) == algebra.one
                    value = algebra.geometric(y, base)
                    value_inverse = algebra.mul(algebra.add(algebra.one, y), algebra.power(inverse, base))
                    assert algebra.mul(value, value_inverse) == algebra.one
                    counts["positional_recovery_box_units"] += 1
            operand_digits = [1, 1, min(3, base - 1)]
            operand = sum(digit * base ** i for i, digit in enumerate(operand_digits))
            product_weight = 1
            for position, digit in enumerate(operand_digits):
                product_weight = recovery.mul(product_weight,
                                               recovery.power(recovery.power(alpha, base ** position), digit))
            assert product_weight == recovery.power(alpha, operand)
            counts["positional_operand_weight_identities"] += 1
            if ell == 5:
                assert recovery.power(alpha, base) != alpha
                counts["recovery_degree_not_dividing_width_controls"] += 1
        records.append({**construction, "D": d, "digit_width": width,
                        "beta": beta, "gamma": gamma, "cube_projection_candidate": candidate,
                        "tested_recovery_degrees": [3, 5, 7]})
    return records


def network_parameters(p, n, operand_length):
    h = 1 << ((p ** 9 - 1).bit_length())
    a, d = 2 * h, 2 * h * h
    total, t, g = 2 * d, d.bit_length() - 1, 3 * p ** 8
    preliminary = ((n + 1 + total - a - 2) // (total - a - 1)
                   if operand_length else (n + d - 1) // d)
    width = adjusted_width(preliminary)
    m = d * ((2 * width + t + d - 1) // d)
    f, c = m + 1, width + t
    u = total * (8 * f + 2 * c)
    s = total * (width + 4 * f + c)
    ua, ut = a // g + 2, total // g + 2
    q = 7 * ua * a * f + ut * (f + 2 * c)
    w = ua * a * (width + 3 * f) + ut * (f + c)
    lam = max(1, (p * u + g * q) // (p ** 3))
    table = w * (2 * g // p) + 1
    return {"D": d, "A": a, "t": t, "preliminary_width": preliminary,
            "width": width, "m": m, "U": u, "S": s,
            "top_volume": p * p * lam, "scaled_table_bound": p * table}


def check_constructor_budgets(counts):
    records = []
    for p in list(range(2, 65)) + [127, 128, 255, 256]:
        n0 = p ** 40
        reference = network_parameters(p, n0, False)
        d, a = reference["D"], reference["A"]
        denominator = 2 * d - a - 1
        lengths = {n0, 2 * n0, 8 * n0, n0 + 1, n0 + d - 1}
        lengths.update(n0 + k * d for k in range(7))
        lengths.update(n0 + k * denominator for k in range(7))
        for n in sorted(lengths):
            for operand_length in (False, True):
                v = network_parameters(p, n, operand_length)
                assert p ** 40 <= n and 2 * d * d <= n and d >= 1 << 19
                assert v["t"] % 2 == 1
                assert 0 <= v["width"] - v["preliminary_width"] <= 5
                assert v["width"] % 6 == 2 and v["m"] % d == 0
                assert 2 * v["width"] + v["t"] <= v["m"]
                assert 216 + 20 * v["t"] < d
                assert 120 + 10 * v["t"] < d
                if operand_length:
                    assert 100000 * d * v["width"] < 50101 * n
                    assert (1 << 15) * d * v["t"] < n
                    if p == 2:
                        assert 10 * v["U"] < 221 * n and 10 * v["S"] < 121 * n
                    else:
                        assert 50 * d * d < 19 * n
                        assert v["U"] < 25 * n and v["S"] < 14 * n
                    assert v["top_volume"] <= 49 * n
                    assert v["scaled_table_bound"] <= 57 * n
                else:
                    assert v["U"] < 45 * n and v["S"] < 25 * n
                    assert v["top_volume"] <= 99 * n
                    assert v["scaled_table_bound"] <= 117 * n
                counts["uniform_constructor_budget_controls"] += 1
                if p in (2, 3, 4, 16, 128) and n == n0:
                    records.append({"P": p, "n_bit_length": n.bit_length(),
                                    "operand_length_constructor": operand_length,
                                    "D_bit_length": d.bit_length(),
                                    "width_increase": v["width"] - v["preliminary_width"],
                                    "all_uniform_bounds_passed": True})
    return records


def coefficient_kernel(field, target_sum, first_sum, second_sum):
    c = target_sum % 3
    assert c != 0
    if first_sum * second_sum % 3 != c:
        return 0
    residue = (first_sum * second_sum - target_sum) % 9
    assert residue % 3 == 0
    return field.power(2, (second_sum + 2 * (residue // 3)) % 3)


def separated_kernel(field, target_sum, first_sum, second_sum):
    beta, c = 2, target_sum % 3
    k = (target_sum - c) // 3
    answer = 0
    for u in (1, 2):
        v = c * u % 3
        assert v in (1, 2) and u * v % 3 == c
        if first_sum % 3 == u and second_sum % 3 == v:
            i, j = (first_sum - u) // 3, (second_sum - v) // 3
            scalar = v + 2 * ((u * v - c) // 3) - 2 * k
            answer ^= field.power(beta, (scalar + 2 * v * i + 2 * u * j) % 3)
    return answer


def check_rank_and_phases(counts):
    field4, field64 = BinaryField(0b111), construct_field(6)
    for candidate in range(2, field64.size):
        eta = field64.power(candidate, (field64.size - 1) // 9)
        if field64.power(eta, 3) != 1:
            break
    else:
        raise AssertionError("missing primitive ninth root")
    assert field64.power(eta, 9) == 1
    beta64 = field64.power(eta, 3)
    assert field64.mul(beta64, beta64) ^ beta64 ^ 1 == 0

    def embed4(value):
        return (value & 1) ^ (beta64 if value & 2 else 0)

    records = []
    for target_sum in (1, 2, 4, 5, 7, 8):
        matrix = [[coefficient_kernel(field4, target_sum, first, second)
                   for second in range(9)] for first in range(9)]
        assert field_rank(field4, matrix) == 2
        counts["exact_rank_two_digit_sum_kernels"] += 1
        c, k = target_sum % 3, (target_sum - target_sum % 3) // 3
        for first, second in product(range(9), repeat=2):
            expected = matrix[first][second]
            assert expected == separated_kernel(field4, target_sum, first, second)
            phase_sum = 0
            for u in (1, 2):
                v = c * u % 3
                exponent = v + 2 * ((u * v - c) // 3) - 2 * k
                scalar = field64.power(beta64, exponent % 3)
                for a, b in product(range(3), repeat=2):
                    first_phase = field64.power(eta, 2 * v + 3 * a)
                    second_phase = field64.power(eta, 2 * u + 3 * b)
                    assert field64.power(first_phase, 3) != 1
                    assert field64.power(second_phase, 3) != 1
                    phase_sum ^= field64.product(scalar,
                                                  field64.power(first_phase, (first - u) % 9),
                                                  field64.power(second_phase, (second - v) % 9))
            assert phase_sum == embed4(expected)
            counts["eighteen_phase_kernel_entries"] += 1
        records.append({"target_digit_sum_mod_9": target_sum, "exact_rank": 2})
    return {"kernels": records, "field64_polynomial_binary": bin(field64.polynomial),
            "primitive_ninth_root": eta, "embedded_cube_root": beta64}


def variable_digits(value, base):
    result = []
    while value:
        result.append(value % base)
        value //= base
    return result or [0]


def check_literal_carry_weights(counts):
    rng, field = random.Random(4102026), BinaryField(0b111)
    samples = []
    for width in (2, 8, 14):
        base = 1 << width
        assert base % 9 == 4
        values = {1, 2, 3, 4, base - 1, base, base + 1, 2 * base + 7, base * base - 1}
        values.update(rng.randrange(1, base * base) for _ in range(12))
        for p, q in product(sorted(values), repeat=2):
            target = p * q
            if target % 3 == 0:
                continue
            left, right = variable_digits(p, base), variable_digits(q, base)
            target_digits = digits(target, base, len(left) + len(right))
            rho = sum(genuine_carries(left, right, base, target_digits))
            sp, sq, sn = sum(left), sum(right), sum(target_digits)
            assert (base - 1) * rho == sp * sq - sn
            expected = field.power(2, (q + 2 * rho) % 3)
            assert expected == coefficient_kernel(field, sn % 9, sp % 9, sq % 9)
            assert expected == separated_kernel(field, sn, sp, sq)
            counts["literal_product_carry_weight_identities"] += 1
            counts["literal_products_with_nonzero_carries"] += rho != 0
            if rho and len(samples) < 6:
                samples.append({"digit_width": width, "p": p, "q": q,
                                "carry_sum": rho, "digit_sums": [sp, sq, sn],
                                "coefficient_in_F4": expected})
    return samples


def check_constant_field_recovery(counts):
    field4 = BinaryField(0b111)
    fields = {ell: construct_field(ell) for ell in (3, 5, 7, 11, 13)}
    records = []
    for initial in (7, 13, 25, 35, 49, 77, 121):
        target, factors = initial, []
        while target > 1:
            n, d = target.bit_length(), 8
            width = adjusted_width((n + d - 1) // d)
            base = 1 << width
            divisors = [value for value in range(1, target + 1) if target % value == 0]
            family = degree_family(n, len(divisors))
            assert set(family) <= set(fields)
            terms = []
            for divisor in divisors:
                rho = sum(genuine_carries(digits(divisor, base, d),
                                           digits(target // divisor, base, d), base,
                                           digits(target, base, 2 * d)))
                coefficient = field4.power(2, (target // divisor + 2 * rho) % 3)
                assert coefficient != 0
                if divisor == 1:
                    assert coefficient == field4.power(2, target % 3)
                else:
                    terms.append((coefficient, divisor, 0))

            def nonzero_prefix(k, prefix):
                counts["constant_field_prefix_queries"] += 1
                selected = [term for term in terms if term[1] >> k == prefix]
                for ell in family:
                    field = fields[ell]
                    for alpha in range(2, field.size):
                        counts["constant_field_recovery_evaluations"] += 1
                        if binary_coefficient_evaluate(field, selected, alpha, 2) != (0, 0):
                            return True
                return False

            prefix = 0
            for k in range(n - 1, -1, -1):
                prefix *= 2
                if not nonzero_prefix(k, prefix):
                    prefix += 1
            assert prefix == min(value for value in divisors if value > 1)
            factors.append(prefix)
            target //= prefix
        assert all(all(prime % divisor for divisor in range(2, prime)) for prime in factors)
        reconstructed = 1
        for prime in factors:
            reconstructed *= prime
        assert reconstructed == initial
        counts["constant_field_complete_factorizations"] += 1
        records.append({"N": initial, "recovered_primes": factors})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts = Counter()
    fields = check_field_interfaces(counts)
    budgets = check_constructor_budgets(counts)
    kernels = check_rank_and_phases(counts)
    carries = check_literal_carry_weights(counts)
    recovery = check_constant_field_recovery(counts)
    result = {"status": "all exact controls passed", "score_improvement_claimed": False,
              "summary": dict(sorted(counts.items())), "root_field_models": fields,
              "constructor_budget_samples": budgets, "kernel_and_phase_models": kernels,
              "literal_carry_samples": carries, "complete_recovery_samples": recovery}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "summary": result["summary"]}, indent=2))


if __name__ == "__main__":
    main()
