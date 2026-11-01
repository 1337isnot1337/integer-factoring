#!/usr/bin/env python3
"""Exact first-jet carry identities and the eighteen F4 input phases."""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import random

from check_carry_fourier_elimination import BinaryField, field_rank, geometric
from check_constant_field_carry_weights import network_parameters, variable_digits
from check_frobenius_sparse_recovery import binary_coefficient_evaluate, construct_field
from check_weighted_dyadic_carry import convolution, digits, genuine_carries


def tags(digit_vector):
    return sum(digit_vector) % 3, sum(i * digit for i, digit in enumerate(digit_vector)) % 3


def kernel(field, c, h, left, right):
    u, r = left
    v, s = right
    if u * v % 3 != c:
        return 0
    return field.power(2, (v + 2 * h + v * r + u * s) % 3)


def phase_expansion(field, c, h, left_digits, right_digits):
    answer = 0
    for u in (1, 2):
        v = c * u % 3
        for a, d in product(range(3), repeat=2):
            exponent = v + 2 * h - a * u - d * v
            exponent += sum((a + v * i) * digit for i, digit in enumerate(left_digits))
            exponent += sum((d + u * j) * digit for j, digit in enumerate(right_digits))
            answer ^= field.power(2, exponent % 3)
    return answer


def check_kernels(counts):
    field = BinaryField(0b111)
    tag_values = list(product(range(3), repeat=2))
    records = []
    for c, h in product((1, 2), range(3)):
        matrix = [[kernel(field, c, h, left, right) for right in tag_values] for left in tag_values]
        assert field_rank(field, matrix) == 2
        counts["exact_first_jet_rank_two_kernels"] += 1
        for left, right in product(tag_values, repeat=2):
            # Two positions realize each tag: z1=derivative, z0=value-z1.
            left_digits = [(left[0] - left[1]) % 3, left[1]]
            right_digits = [(right[0] - right[1]) % 3, right[1]]
            assert tags(left_digits) == left and tags(right_digits) == right
            assert phase_expansion(field, c, h, left_digits, right_digits) == kernel(field, c, h, left, right)
            counts["all_tag_eighteen_phase_identities"] += 1
        records.append({"target_residue": c, "target_derivative_tag": h, "rank": 2})
    return records


def check_literal_products(counts):
    rng, field = random.Random(7312026), BinaryField(0b111)
    records = []
    for width in (2, 4, 6, 8, 10, 12):
        base = 1 << width
        assert base % 3 == 1
        values = {1, 2, 3, 4, base - 1, base, base + 1, 2 * base + 5, base * base - 1}
        values.update(rng.randrange(1, base * base) for _ in range(10))
        for p, q in product(sorted(values), repeat=2):
            target = p * q
            left, right = variable_digits(p, base), variable_digits(q, base)
            length = len(left) + len(right)
            target_digits = digits(target, base, length)
            carries = genuine_carries(left, right, base, target_digits)
            rho = sum(carries)
            coefficients = convolution(left, right, length)
            padded_carries = [0, *carries, 0]
            assert all(coefficients[j] == target_digits[j] + base * padded_carries[j + 1] - padded_carries[j]
                       for j in range(length))
            derivative_product = sum(i * digit for i, digit in enumerate(coefficients))
            derivative_target = sum(i * digit for i, digit in enumerate(target_digits))
            derivative_carry = sum(i * carry for i, carry in enumerate(carries))
            assert derivative_product == derivative_target + (base - 1) * derivative_carry - rho
            p0, p1 = tags(left)
            q0, q1 = tags(right)
            n0, n1 = tags(target_digits)
            assert rho % 3 == (n1 - p1 * q0 - p0 * q1) % 3
            counts["literal_integer_derivative_identities"] += 1
            if target % 3 == 0:
                continue
            expected = field.power(2, (q + 2 * rho) % 3)
            assert expected == kernel(field, n0, n1, (p0, p1), (q0, q1))
            assert expected == phase_expansion(field, n0, n1, left, right)
            counts["literal_carry_phase_equalities"] += 1
            counts["nonzero_carry_phase_equalities"] += rho != 0
            if rho and len([r for r in records if r["digit_width"] == width]) < 2:
                records.append({"digit_width": width, "width_mod_6": width % 6,
                                "p": p, "q": q, "carry_sum": rho,
                                "first_jet_tags": [[p0, p1], [q0, q1], [n0, n1]]})
    return records


def check_prefix_contractions(counts):
    field4 = BinaryField(0b111)
    recovery_fields = {ell: construct_field(ell) for ell in (3, 5, 7)}
    records = []
    for target in (5, 7, 25, 35, 49, 77):
        for width in (2, 4, 6):
            base, n = 1 << width, target.bit_length()
            divisors = [p for p in range(1, target + 1) if target % p == 0]
            length = 2 * max(1, (n + width - 1) // width)
            target_digits = digits(target, base, length)
            c, h = tags(target_digits)
            for k in range(n + 1):
                for prefix in sorted({0, 1, *(p >> k for p in divisors)}):
                    selected = [p for p in divisors if p >> k == prefix]
                    direct_terms = []
                    phase_terms = [[] for _ in range(18)]
                    for p in selected:
                        q = target // p
                        left, right = digits(p, base, length // 2), digits(q, base, length // 2)
                        rho = sum(genuine_carries(left, right, base, target_digits))
                        direct_terms.append((field4.power(2, (q + 2 * rho) % 3), p, 0))
                        index = 0
                        for u in (1, 2):
                            v = c * u % 3
                            for a, d in product(range(3), repeat=2):
                                exponent = v + 2 * h - a * u - d * v
                                exponent += sum((a + v * i) * digit for i, digit in enumerate(left))
                                exponent += sum((d + u * j) * digit for j, digit in enumerate(right))
                                phase_terms[index].append((field4.power(2, exponent % 3), p, 0))
                                index += 1
                    for ell, field in recovery_fields.items():
                        for alpha in (2, field.size - 1):
                            direct = binary_coefficient_evaluate(field, direct_terms, alpha, 2)
                            combined = [0, 0]
                            for terms in phase_terms:
                                value = binary_coefficient_evaluate(field, terms, alpha, 2)
                                combined[0] ^= value[0]
                                combined[1] ^= value[1]
                            assert direct == tuple(combined)
                            counts["complete_prefix_eighteen_contractions"] += 1
            records.append({"N": target, "digit_width": width,
                            "target_first_jet": [c, h], "tested_recovery_degrees": [3, 5, 7]})
    return records


def check_even_width_budgets(counts):
    records = []
    for p in list(range(2, 65)) + [127, 128, 255, 256]:
        n0 = p ** 40
        h = 1 << ((p ** 9 - 1).bit_length())
        a, d = 2 * h, 2 * h * h
        total, t, g = 2 * d, d.bit_length() - 1, 3 * p ** 8
        for n in (n0, n0 + 1, n0 + d - 1, n0 + d, 2 * n0, 8 * n0):
            for operand_length in (False, True):
                earlier = network_parameters(p, n, operand_length)
                preliminary = earlier["preliminary_width"]
                width = preliminary + preliminary % 2
                m = d * ((2 * width + t + d - 1) // d)
                f, c = m + 1, width + t
                u, s = total * (8 * f + 2 * c), total * (width + 4 * f + c)
                ua, ut = a // g + 2, total // g + 2
                q = 7 * ua * a * f + ut * (f + 2 * c)
                w = ua * a * (width + 3 * f) + ut * (f + c)
                lam = max(1, (p * u + g * q) // p ** 3)
                table = w * (2 * g // p) + 1
                assert 0 <= width - preliminary <= 1
                assert width % 2 == 0 and width <= earlier["width"]
                assert m <= earlier["m"] and u <= earlier["U"] and s <= earlier["S"]
                assert p * p * lam <= (49 if operand_length else 99) * n
                assert p * table <= (57 if operand_length else 117) * n
                counts["smaller_even_width_constructor_controls"] += 1
                if p in (2, 3, 16, 128) and n == n0:
                    records.append({"P": p, "n_bit_length": n.bit_length(),
                                    "operand_length_constructor": operand_length,
                                    "width_increase": width - preliminary,
                                    "width_mod_6": width % 6})
    field = BinaryField(0b111)
    for width in (2, 4, 6, 8):
        base = 1 << width
        assert geometric(field, 1, base) == 0
        for u, d in product((1, 2), range(3)):
            phases = [field.power(2, (d + u * j) % 3) for j in range(3)]
            assert phases.count(1) == 1
            assert [geometric(field, phase, base) for phase in phases].count(0) == 1
            counts["alternate_zero_factor_controls"] += 1
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts = Counter()
    kernels = check_kernels(counts)
    products = check_literal_products(counts)
    prefixes = check_prefix_contractions(counts)
    budgets = check_even_width_budgets(counts)
    result = {"status": "all exact controls passed", "score_improvement_claimed": False,
              "summary": dict(sorted(counts.items())), "first_jet_kernels": kernels,
              "literal_product_samples": products, "prefix_contraction_models": prefixes,
              "even_width_budget_samples": budgets}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "summary": result["summary"]}, indent=2))


if __name__ == "__main__":
    main()
