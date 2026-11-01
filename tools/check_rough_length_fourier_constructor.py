#!/usr/bin/env python3
"""Exact controls for the rough-cofactor, shorter-transform constructor.

The complete time bounds use the proofs in the two companion notes. These
finite controls do not benchmark the asymptotic prime search or factoring bound.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from math import isqrt, prod
from pathlib import Path
import json

from check_small_word_fourier_constructor import (
    Parameters, carry_digits, dft_product, full_dft, integer_convolution,
    interpolate, is_prime, make_crt_parameters,
)


def ceildiv(a, b):
    return (a + b - 1) // b


def branches(m, s):
    return [(k, l) for k in range(s + 1, m - s + 1)
            for l in (m - k, m + 1 - k) if s + 1 <= l <= m - s]


def length_controls():
    parameter_cases = branch_cases = beyond_half_rows = 0
    for n in range(2, 129):
        for a in (2, 4, 8):
            t = a * a
            if t > n + 1:
                continue
            d = ceildiv(n + 1, t)
            vertices = 2
            while vertices <= a:
                block = a // vertices
                s = t * d // vertices
                assert s * vertices == t * d
                for m in range(1, n + 1):
                    choices = branches(m, s)
                    assert len(choices) == len(set(choices)) <= 2 * m
                    for k, l in choices:
                        hx, hy = ceildiv(k, d), ceildiv(l, d)
                        assert k + l in (m, m + 1)
                        assert hx + hy - 1 <= ceildiv(k + l, d) <= t
                        assert min(hx, hy) <= t // 2
                        assert max(hx, hy) <= t - t // vertices
                        assert max(ceildiv(hx, a), ceildiv(hy, a)) <= a - block
                        assert ((1 << k) - 1) * ((1 << l) - 1) < 1 << (t * d)
                        beyond_half_rows += max(hx, hy) > t // 2
                        branch_cases += 1
                    if m <= 2 * s:
                        assert not choices
                    parameter_cases += 1
                vertices *= 2
    assert beyond_half_rows > 0
    return {"padded_length_parameter_cases": parameter_cases,
            "exact_length_branches": branch_cases,
            "branches_with_input_longer_than_old_half_capacity": beyond_half_rows}


def input_word(value, bits, p):
    assert value.bit_length() == bits
    result = []
    for _ in range(ceildiv(bits, p.d)):
        value, digit = divmod(value, p.base)
        result.append(digit)
    assert not value
    return result


def variable_input_messages(p, word):
    return [[pow(p.omega, r * k, p.q) * sum(
        word[j] * pow(p.xi, ((j - r) // p.a) * k, p.q)
        for j in range(r, len(word), p.a)
    ) % p.q for k in range(p.a)] for r in range(p.a)]


def variable_grid_multiply(p, x, y):
    x_messages, y_messages = variable_input_messages(p, x), variable_input_messages(p, y)
    x_hat, y_hat = [0] * p.t, [0] * p.t
    inverse_a = pow(p.a, -1, p.q)
    z = [[0] * p.a for _ in range(p.a)]
    for k in range(p.a):
        products = []
        for v in range(p.a):
            x_hat[k + p.a * v] = sum(x_messages[r][k] * pow(p.xi, r * v, p.q)
                                        for r in range(p.a)) % p.q
            y_hat[k + p.a * v] = sum(y_messages[r][k] * pow(p.xi, r * v, p.q)
                                        for r in range(p.a)) % p.q
            products.append(x_hat[k + p.a * v] * y_hat[k + p.a * v] % p.q)
        for r in range(p.a):
            z[k][r] = inverse_a * sum(products[v] * pow(p.xi, (-r * v) % p.a, p.q)
                                      for v in range(p.a)) % p.q
    coefficients = [0] * p.t
    for r in range(p.a):
        for i in range(p.a):
            coefficients[r + p.a * i] = inverse_a * sum(
                z[k][r] * pow(p.omega, (-r * k) % p.t, p.q)
                * pow(p.xi, (-i * k) % p.a, p.q) for k in range(p.a)
            ) % p.q
    return coefficients, x_hat, y_hat


def recover_variable_row(p, r, length, boundary):
    if len(boundary) < length or any(not 0 <= value < p.q for value in boundary.values()):
        return None
    chosen = sorted(boundary)[:length]
    points = [pow(p.xi, k, p.q) for k in chosen]
    values = [boundary[k] * pow(pow(p.omega, r * k, p.q), -1, p.q) % p.q for k in chosen]
    coefficients = interpolate(points, values, p.q)
    if any(coefficient >= p.base for coefficient in coefficients):
        return None
    for k, value in boundary.items():
        expected = pow(p.omega, r * k, p.q) * sum(
            coefficient * pow(p.xi, i * k, p.q) for i, coefficient in enumerate(coefficients)
        ) % p.q
        if value != expected:
            return None
    return coefficients


def transform_interpolation_controls():
    pairs = rows = longer_rows = 0
    records = []
    for a, n, vertices in ((4, 27, 4), (8, 120, 4), (8, 120, 8), (16, 500, 4)):
        d = ceildiv(n + 1, a * a)
        p, primes = make_crt_parameters(a, d)
        s, block = p.t * d // vertices, a // vertices
        chosen_lengths = sorted({(s + 1, n - s), (n - s, s + 1), (n // 2, n - n // 2)})
        local_pairs = 0
        for k, l in chosen_lengths:
            assert (k, l) in branches(n, s)
            candidates = [((1 << (k - 1)), (1 << (l - 1))),
                          ((1 << k) - 1, (1 << l) - 1),
                          ((1 << (k - 1)) + (1 << (k - 2)) + 1,
                           (1 << (l - 1)) + 1)]
            for xv, yv in candidates:
                x, y = input_word(xv, k, p), input_word(yv, l, p)
                coefficients, x_hat, y_hat = variable_grid_multiply(p, x, y)
                assert x_hat == full_dft(p, x) and y_hat == full_dft(p, y)
                assert coefficients == dft_product(p, x, y)
                assert coefficients == integer_convolution(x, y, p.t)
                assert max(coefficients) <= min(len(x), len(y)) * (p.base - 1)**2
                assert max(coefficients) <= p.half * (p.base - 1)**2 < p.q
                output = carry_digits(p, coefficients)
                assert sum(digit * p.base**j for j, digit in enumerate(output)) == xv * yv
                for word in (x, y):
                    messages = variable_input_messages(p, word)
                    for r in range(a):
                        wanted = word[r::a]
                        external = {k: messages[r][k] for k in range(a) if k // block != r // block}
                        assert len(wanted) <= len(external) == a - block
                        assert recover_variable_row(p, r, len(wanted), external) == wanted
                        longer_rows += len(wanted) > p.h
                        rows += 1
                pairs += 1
                local_pairs += 1
        records.append({"A": a, "padded_input_bits": n, "vertices": vertices,
                        "digit_bits": d, "threshold_exponent": s,
                        "arithmetic_modulus": str(p.q), "crt_primes": primes,
                        "length_pairs": chosen_lengths, "transform_pairs": local_pairs})
    assert longer_rows > 0
    p, _ = make_crt_parameters(4, 1)
    assert recover_variable_row(p, 0, 0, {1: 0, 2: 0, 3: 0}) == []
    assert recover_variable_row(p, 0, 0, {1: 0, 2: 1, 3: 0}) is None
    return {"transform_pairs": pairs, "recovered_input_rows": rows,
            "recovered_rows_longer_than_old_half_row": longer_rows,
            "empty_row_controls": 2, "cases": records}


def divisor_pairs(value):
    divisors = set()
    for trial in range(1, isqrt(value) + 1):
        if value % trial == 0:
            divisors.update((trial, value // trial))
    return [(divisor, value // divisor) for divisor in sorted(divisors)]


def rough_coverage_controls():
    values = list(range(2, 4097)) + [257**2, 257 * 263, 257 * 1000003, 3 * 257 * 1000003]
    cases = pair_cases = composite_cofactors = 0
    selected = []
    for value in values:
        n = value.bit_length()
        a, vertices = 4, 4
        d = ceildiv(n + 1, a * a)
        s = a * a * d // vertices
        primes = [p for p in range(2, (1 << s) + 1) if is_prime(p)]
        remaining = value
        removed = []
        for prime in primes:
            while remaining % prime == 0:
                remaining //= prime
                removed.append(prime)
        assert prod(removed) * remaining == value
        assert all(remaining % prime != 0 for prime in primes)
        expected = divisor_pairs(remaining)
        if remaining == 1:
            found = [(1, 1)]
            proper = []
        else:
            m = remaining.bit_length()
            choices = branches(m, s)
            counts = Counter(choices)
            proper = [(x, y) for x, y in expected if x > 1 and y > 1]
            for x, y in proper:
                assert x > 1 << s and y > 1 << s
                assert counts[x.bit_length(), y.bit_length()] == 1
                assert max(ceildiv(x.bit_length(), d), ceildiv(y.bit_length(), d)) <= a * a - a * a // vertices
            found = [(1, remaining), (remaining, 1)] + proper
            if m <= 2 * s:
                assert not proper and is_prime(remaining)
        assert sorted(found) == expected
        assert sum((x + 3 * y) % 101 for x, y in found) % 101 == sum(
            (x + 3 * y) % 101 for x, y in expected
        ) % 101
        composite_cofactors += bool(proper)
        cases += 1
        pair_cases += len(proper)
        if value in values[-4:]:
            selected.append({"input": value, "threshold_exponent": s, "rough_cofactor": remaining,
                             "removed_prime_factors": removed, "proper_pairs": proper})
    assert composite_cofactors > 10
    return {"inputs": cases, "proper_factor_pairs": pair_cases,
            "composite_rough_cofactors": composite_cofactors, "selected_cases": selected}


def carry_pair_controls():
    p, _ = make_crt_parameters(4, 1)
    vertices, s = 4, 4
    cases = accepted = 0
    records = []
    for target in (323, 361, 437, 629, 667, 899):
        n = target.bit_length()
        assert p.d == ceildiv(n + 1, p.t) and s == p.t * p.d // vertices
        assert all(target % prime for prime in range(2, 17) if is_prime(prime))
        found = []
        for k, l in branches(n, s):
            for x in range(1 << (k - 1), 1 << k):
                for y in range(1 << (l - 1), 1 << l):
                    coefficients = integer_convolution(input_word(x, k, p), input_word(y, l, p), p.t)
                    states = [0] * (p.carry_cap + 1)
                    states[0] = 1
                    rest = target
                    for coefficient in coefficients:
                        rest, target_digit = divmod(rest, p.base)
                        following = [0] * len(states)
                        for incoming, multiplicity in enumerate(states):
                            outgoing, digit = divmod(coefficient + incoming, p.base)
                            if outgoing <= p.carry_cap and digit == target_digit:
                                following[outgoing] += multiplicity
                        states = following
                    assert states[0] == int(x * y == target)
                    if states[0]:
                        found.append((x, y))
                        accepted += 1
                    cases += 1
        assert sorted(found) == [(x, y) for x, y in divisor_pairs(target) if x > 1 and y > 1]
        records.append({"target": target, "proper_pairs": found})
    return {"candidate_pair_target_cases": cases, "accepted_proper_pairs": accepted, "cases": records}


def actual_policy_controls():
    lengths = groupings = 0
    for n in range(128, 65537):
        mu = (n + 1).bit_length()
        g = (mu + 1).bit_length()
        a = 1
        while (2 * a)**2 * mu * g <= n:
            a *= 2
        assert a >= 2
        t, d = a * a, ceildiv(n + 1, a * a)
        assert t * mu * g <= n < 4 * t * mu * g
        assert n + 1 <= t * d <= n + t
        assert mu * g <= d <= 4 * mu * g + 1
        vertices = 2
        while vertices <= a:
            s = t * d // vertices
            assert n <= s * vertices and s <= n
            assert t * d == s * vertices
            assert (t - t // vertices) * d == t * d - s > n - s
            groupings += 1
            vertices *= 2
        lengths += 1
    metadata = []
    for n in (128, 256, 1024, 4096, 65536, 1048576):
        mu = (n + 1).bit_length()
        g = (mu + 1).bit_length()
        a = 1
        while (2 * a)**2 * mu * g <= n:
            a *= 2
        d = ceildiv(n + 1, a * a)
        p, primes = make_crt_parameters(a, d)
        vertices = 2
        while (2 * vertices)**4 <= a:
            vertices *= 2
        block = a // vertices
        extra = 2 * a - 1 if vertices == 2 else a
        width = (6 * block * block + extra) * p.f
        assert vertices**2 * width == 6 * p.t * p.f + vertices**2 * extra * p.f
        assert p.f - p.carry_bits >= d + 2
        metadata.append({"input_bits": n, "A": a, "digit_bits": d, "residue_bits": p.f,
                         "transform_word_volume": p.t * p.f,
                         "vertices": vertices, "threshold_exponent": p.t * d // vertices,
                         "column_clique_volume": vertices**2 * width,
                         "number_of_crt_primes": len(primes), "largest_crt_prime": primes[-1],
                         "prime_candidates_tested": p.candidates,
                         "space_efficient_grouping": vertices**2 <= a})
    return {"input_lengths": lengths, "threshold_groupings": groupings,
            "actual_metadata_cases": len(metadata), "metadata": metadata,
            "scope": "Metadata samples only; the constant twelve is an asymptotic proof, not a finite measured bound."}


def insufficient_interpolation_control():
    p, _ = make_crt_parameters(2, 2)
    first, second = [0, 0, 2], [1, 0, 3]
    assert sum(digit * p.base**j for j, digit in enumerate(first)).bit_length() == 6
    assert sum(digit * p.base**j for j, digit in enumerate(second)).bit_length() == 6
    first_messages = variable_input_messages(p, first)
    second_messages = variable_input_messages(p, second)
    assert first != second and first_messages[0][1] == second_messages[0][1]
    assert first_messages[0][0] != second_messages[0][0]
    assert recover_variable_row(p, 0, 2, {1: first_messages[0][1]}) is None
    return {"A": 2, "vertices": 2, "digit_bits": 2,
            "same_exact_factor_bit_length": 6,
            "first_digit_word": first, "second_digit_word": second,
            "same_external_message": first_messages[0][1],
            "reason": "Two input-row coefficients but only one prescribed external output."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    dependency = Path(__file__).with_name("check_small_word_fourier_constructor.py")
    result = {"status": "passed", "scope": "Exact constructor controls; factoring score remains one half.",
              "length_capacities": length_controls(),
              "transforms_and_interpolation": transform_interpolation_controls(),
              "rough_cofactor_coverage": rough_coverage_controls(),
              "carry_pair_bijection": carry_pair_controls(),
              "actual_parameter_policy": actual_policy_controls(),
              "missing_capacity_negative_control": insufficient_interpolation_control(),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "arithmetic_helpers_sha256": hashlib.sha256(dependency.read_bytes()).hexdigest()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"],
                      "length_branches": result["length_capacities"]["exact_length_branches"],
                      "transform_pairs": result["transforms_and_interpolation"]["transform_pairs"],
                      "interpolated_rows": result["transforms_and_interpolation"]["recovered_input_rows"],
                      "rough_inputs": result["rough_cofactor_coverage"]["inputs"],
                      "carry_pair_targets": result["carry_pair_bijection"]["candidate_pair_target_cases"],
                      "policy_lengths": result["actual_parameter_policy"]["input_lengths"],
                      "metadata_cases": result["actual_parameter_policy"]["actual_metadata_cases"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
