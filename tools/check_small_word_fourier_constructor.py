#!/usr/bin/env python3
"""Exact small controls for the logarithmic-word Fourier constructor.

The transform tests use small independently chosen A,d values. They check
the identities and message layout, not the asymptotic running time of a
Linnik prime search. The latter is a cited theorem, not a finite benchmark.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from math import gcd, isqrt, prod
from pathlib import Path


@lru_cache(maxsize=None)
def is_prime(value: int) -> bool:
    return value >= 2 and all(value % trial for trial in range(2, isqrt(value) + 1))


@dataclass(frozen=True)
class Parameters:
    a: int
    d: int
    q: int
    omega: int
    xi: int
    candidates: int
    nonresidue: int

    @property
    def t(self) -> int:
        return self.a * self.a

    @property
    def h(self) -> int:
        return self.a // 2

    @property
    def half(self) -> int:
        return self.t // 2

    @property
    def base(self) -> int:
        return 1 << self.d

    @property
    def f(self) -> int:
        return (self.q - 1).bit_length()

    @property
    def carry_cap(self) -> int:
        return self.half * (self.base - 1)

    @property
    def carry_bits(self) -> int:
        return (self.half - 1).bit_length() + self.d


def make_parameters(a: int, d: int) -> Parameters:
    assert a >= 2 and a & (a - 1) == 0 and d >= 1
    t = a * a
    modulus = t << (2 * d)
    candidate, candidates = modulus + 1, 1
    while not is_prime(candidate):
        candidate += modulus
        candidates += 1
    nonresidue = next(g for g in range(2, candidate)
                      if pow(g, (candidate - 1) // 2, candidate) == candidate - 1)
    omega = pow(nonresidue, (candidate - 1) // t, candidate)
    assert pow(omega, t, candidate) == 1
    assert pow(omega, t // 2, candidate) == candidate - 1
    result = Parameters(a, d, candidate, omega, pow(omega, a, candidate), candidates, nonresidue)
    assert candidate > modulus
    assert result.carry_cap < 1 << result.carry_bits
    assert result.f - result.carry_bits >= d + 2
    return result


def make_crt_parameters(a: int, d: int) -> tuple[Parameters, list[int]]:
    assert a >= 2 and a & (a - 1) == 0 and d >= 1
    t = a * a
    threshold = t << (2 * d)
    modulus, omega, candidate, trials = 1, 0, 1 + t, 0
    primes = []
    while modulus <= threshold:
        trials += 1
        if is_prime(candidate):
            nonresidue = next(g for g in range(2, candidate)
                              if pow(g, (candidate - 1) // 2, candidate) == candidate - 1)
            local_root = pow(nonresidue, (candidate - 1) // t, candidate)
            assert pow(local_root, t // 2, candidate) == candidate - 1
            adjustment = (local_root - omega) * pow(modulus, -1, candidate) % candidate
            omega += modulus * adjustment
            modulus *= candidate
            assert 0 <= omega < modulus
            primes.append(candidate)
        candidate += t
    result = Parameters(a, d, modulus, omega, pow(omega, a, modulus), trials, 0)
    assert modulus == prod(primes)
    assert modulus // primes[-1] <= threshold < modulus
    assert pow(omega, t // 2, modulus) == modulus - 1
    assert result.f - result.carry_bits >= d + 2
    return result, primes


def digits(value: int, p: Parameters) -> list[int]:
    assert 0 <= value < p.base**p.half
    result = []
    for _ in range(p.half):
        value, digit = divmod(value, p.base)
        result.append(digit)
    assert value == 0
    return result


def input_messages(p: Parameters, word: list[int]) -> list[list[int]]:
    return [
        [pow(p.omega, r * k, p.q) * sum(
            word[r + p.a * i] * pow(p.xi, i * k, p.q)
            for i in range(p.h)
        ) % p.q for k in range(p.a)]
        for r in range(p.a)
    ]


def grid_multiply(p: Parameters, x: list[int], y: list[int]) -> tuple[list[int], list[int], list[int]]:
    x_messages, y_messages = input_messages(p, x), input_messages(p, y)
    x_hat, y_hat = [0] * p.t, [0] * p.t
    z = [[0] * p.a for _ in range(p.a)]
    inverse_a = pow(p.a, -1, p.q)
    for k in range(p.a):
        products = []
        for v in range(p.a):
            x_hat[k + p.a * v] = sum(
                x_messages[r][k] * pow(p.xi, r * v, p.q) for r in range(p.a)
            ) % p.q
            y_hat[k + p.a * v] = sum(
                y_messages[r][k] * pow(p.xi, r * v, p.q) for r in range(p.a)
            ) % p.q
            products.append(x_hat[k + p.a * v] * y_hat[k + p.a * v] % p.q)
        for r in range(p.a):
            z[k][r] = inverse_a * sum(
                products[v] * pow(p.xi, (-r * v) % p.a, p.q) for v in range(p.a)
            ) % p.q
    coefficients = [0] * p.t
    for r in range(p.a):
        for i in range(p.a):
            coefficients[r + p.a * i] = inverse_a * sum(
                z[k][r] * pow(p.omega, (-r * k) % p.t, p.q)
                * pow(p.xi, (-i * k) % p.a, p.q) for k in range(p.a)
            ) % p.q
    return coefficients, x_hat, y_hat


def full_dft(p: Parameters, word: list[int]) -> list[int]:
    return [sum(digit * pow(p.omega, j * frequency, p.q)
                for j, digit in enumerate(word)) % p.q for frequency in range(p.t)]


def dft_product(p: Parameters, x: list[int], y: list[int]) -> list[int]:
    x_hat, y_hat = full_dft(p, x), full_dft(p, y)
    products = [a * b % p.q for a, b in zip(x_hat, y_hat)]
    inverse_t = pow(p.t, -1, p.q)
    return [inverse_t * sum(products[k] * pow(p.omega, (-j * k) % p.t, p.q)
                            for k in range(p.t)) % p.q for j in range(p.t)]


def integer_convolution(x: list[int], y: list[int], length: int) -> list[int]:
    result = [0] * length
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            result[i + j] += xi * yj
    return result


def carry_digits(p: Parameters, coefficients: list[int]) -> list[int]:
    carry = 0
    output = []
    for coefficient in coefficients:
        carry, digit = divmod(coefficient + carry, p.base)
        assert 0 <= carry <= p.carry_cap
        output.append(digit)
    assert carry == 0
    return output


def word_patterns(p: Parameters) -> list[list[int]]:
    zero = [0] * p.half
    one = [1] + [0] * (p.half - 1)
    top = [0] * (p.half - 1) + [p.base - 1]
    maximum = [p.base - 1] * p.half
    alternating = [p.base - 1 if j % 2 == 0 else 0 for j in range(p.half)]
    mixed = [(3 * j * j + j + 1) % p.base for j in range(p.half)]
    return list({tuple(word): word for word in (zero, one, top, maximum, alternating, mixed)}.values())


def transform_controls() -> dict:
    cases = 0
    records = []
    for a, d in ((2, 1), (2, 2), (2, 3), (4, 1), (4, 2), (8, 1), (8, 2), (16, 1)):
        p = make_parameters(a, d)
        patterns = word_patterns(p)
        pairs = [(word, patterns[(i + 1) % len(patterns)]) for i, word in enumerate(patterns)]
        pairs += [(patterns[-2], patterns[-2]), (patterns[-1], patterns[-1])]
        for x, y in pairs:
            coefficients, x_hat, y_hat = grid_multiply(p, x, y)
            assert x_hat == full_dft(p, x)
            assert y_hat == full_dft(p, y)
            assert coefficients == dft_product(p, x, y)
            integer = integer_convolution(x, y, p.t)
            assert max(integer) <= p.half * (p.base - 1)**2 < p.q
            assert coefficients == integer
            output = carry_digits(p, coefficients)
            xv = sum(digit * p.base**j for j, digit in enumerate(x))
            yv = sum(digit * p.base**j for j, digit in enumerate(y))
            assert sum(digit * p.base**j for j, digit in enumerate(output)) == xv * yv
            cases += 1
        records.append({
            "A": a, "digit_bits": d, "q": p.q, "omega": p.omega, "xi": p.xi,
            "residue_bits": p.f, "carry_cap": p.carry_cap,
            "spare_control_bits": p.f - p.carry_bits,
            "prime_candidates_tested": p.candidates, "nonresidue": p.nonresidue,
            "word_pairs": len(pairs),
        })
    return {"word_pairs": cases, "parameters": records}


def interpolate(points: list[int], values: list[int], q: int) -> list[int]:
    result = [0] * len(points)
    for i, point in enumerate(points):
        polynomial, denominator = [1], 1
        for j, other in enumerate(points):
            if j == i:
                continue
            product_coefficients = [0] * (len(polynomial) + 1)
            for power, coefficient in enumerate(polynomial):
                product_coefficients[power] -= other * coefficient
                product_coefficients[power + 1] += coefficient
            polynomial = [coefficient % q for coefficient in product_coefficients]
            denominator = denominator * (point - other) % q
        scale = values[i] * pow(denominator, -1, q) % q
        result = [(old + scale * coefficient) % q for old, coefficient in zip(result, polynomial)]
    return result


def recover_input(p: Parameters, r: int, boundary: dict[int, int]) -> list[int] | None:
    if len(boundary) < p.h or any(not 0 <= value < p.q for value in boundary.values()):
        return None
    chosen = sorted(boundary)[:p.h]
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


def interpolation_controls() -> dict:
    cases = 0
    for a, d in ((2, 1), (2, 2), (4, 1), (4, 2), (8, 1), (8, 2), (16, 1)):
        p = make_parameters(a, d)
        subsets = list(combinations(range(a), p.h)) if a <= 8 else [
            tuple(range(p.h)), tuple(range(p.h, a)), tuple(range(0, a, 2)),
        ]
        for word in word_patterns(p):
            messages = input_messages(p, word)
            for r in sorted({0, 1, a - 1}):
                expected = [word[r + a * i] for i in range(p.h)]
                for selected in subsets:
                    assert recover_input(p, r, {k: messages[r][k] for k in selected}) == expected
                    cases += 1

    p = make_parameters(4, 2)
    invalid_digit = {k: p.base for k in range(p.a)}
    assert recover_input(p, 0, invalid_digit) is None
    invalid_residue = {k: 0 for k in range(p.a)}
    invalid_residue[0] = p.q
    assert recover_input(p, 0, invalid_residue) is None
    inconsistent_extra = {k: 0 for k in range(p.a)}
    inconsistent_extra[p.a - 1] = 1
    assert recover_input(p, 0, inconsistent_extra) is None
    return {
        "recovered_input_gates": cases,
        "invalid_digit_rejected": True,
        "invalid_residue_rejected": True,
        "inconsistent_extra_boundary_rejected": True,
    }


def boundary_controls() -> dict:
    cases = codec_cases = 0
    for a, d in ((2, 2), (4, 2), (8, 1), (16, 1), (32, 1)):
        p = make_parameters(a, d)
        vertices = 2
        while vertices <= a:
            block = a // vertices
            crossings: Counter[tuple[int, int]] = Counter()

            def edge(left: int, right: int) -> None:
                if left != right:
                    crossings[tuple(sorted((left, right)))] += 1

            for r in range(a):
                for k in range(a):
                    edge(r // block, k // block)
                    edge(r // block, k // block)
                    edge(k // block, r // block)
                for i in range(a):
                    edge(r // block, i // block)
            for j in range(p.t - 1):
                edge((j // a) // block, ((j + 1) // a) // block)
            for left in range(vertices):
                for right in range(left + 1, vertices):
                    assert crossings[left, right] == 8 * block * block + int(right == left + 1)
            width = (8 * block * block + 1) * p.f
            assert vertices * vertices * width <= 9 * p.t * p.f
            assert (vertices - 1) * width <= 9 * p.t * p.f // vertices
            cases += 1
            vertices *= 2

        for carry in range(p.carry_cap + 1):
            for state in (0, 1, 2, (1 << (p.f - p.carry_bits)) - 1):
                encoded = carry + (state << p.carry_bits)
                assert encoded < 1 << p.f
                assert encoded & ((1 << p.carry_bits) - 1) == carry
                assert encoded >> p.carry_bits == state
                codec_cases += 1
    return {"groupings": cases, "carry_control_codecs": codec_cases}


def column_boundary_controls() -> dict:
    cases = segment_cases = efficient_cases = 0
    records = []
    for a, d in ((2, 2), (4, 2), (8, 1), (16, 1), (32, 1), (64, 1)):
        p, _ = make_crt_parameters(a, d)
        vertices = 2
        while vertices <= a:
            block = a // vertices
            crossings: Counter[tuple[int, int]] = Counter()
            carry_crossings: Counter[tuple[int, int]] = Counter()

            def edge(left: int, right: int, carry: bool = False) -> None:
                if left != right:
                    pair = tuple(sorted((left, right)))
                    crossings[pair] += 1
                    if carry:
                        carry_crossings[pair] += 1

            for r in range(a):
                for k in range(a):
                    edge(r // block, k // block)
                    edge(r // block, k // block)
                    edge(k // block, r // block)
                for _ in range(a):
                    edge(r // block, r // block)
            for j in range(p.t - 1):
                edge((j % a) // block, ((j + 1) % a) // block, carry=True)
            assert sum(carry_crossings.values()) == vertices * a - 1
            for left in range(vertices):
                for right in range(left + 1, vertices):
                    expected_carries = a * int(right == left + 1) + (a - 1) * int(
                        left == 0 and right == vertices - 1
                    )
                    assert carry_crossings[left, right] == expected_carries
                    assert crossings[left, right] == 6 * block * block + expected_carries
            extra = 2 * a - 1 if vertices == 2 else a
            width = (6 * block * block + extra) * p.f
            mass = vertices * vertices * width
            assert mass == 6 * p.t * p.f + vertices * vertices * extra * p.f
            if vertices * vertices <= a:
                assert (vertices - 1) * width <= 8 * p.t * p.f // vertices
                assert mass <= (8 if vertices == 2 else 7) * p.t * p.f
                efficient_cases += 1
            for group in range(vertices):
                segments = [list(range(row * a + group * block, row * a + (group + 1) * block))
                            for row in range(a)]
                members = [j for j in range(p.t) if (j % a) // block == group]
                assert [j for segment in segments for j in segment] == members
                for segment in segments:
                    assert segment[0] == 0 or ((segment[0] - 1) % a) // block != group
                    assert segment[-1] == p.t - 1 or ((segment[-1] + 1) % a) // block != group
                    segment_cases += 1
            for group in range(vertices - 1):
                source = (group + 1) * block - 1
                assert source // a == 0
                assert (source % a) // block == group
                assert ((source + 1) % a) // block == group + 1
            records.append({"A": a, "vertices": vertices, "width": width, "volume": mass,
                            "carry_crossings": sum(carry_crossings.values()),
                            "P2_at_most_A": vertices * vertices <= a})
            cases += 1
            vertices *= 2
    return {"groupings": cases, "carry_segments": segment_cases,
            "space_efficient_groupings": efficient_cases, "cases": records}


def factor_pair_controls() -> dict:
    p = make_parameters(2, 2)
    word_limit = p.base**p.half
    counts = Counter()
    cases = 0
    for x in range(word_limit):
        for y in range(word_limit):
            coefficients, _, _ = grid_multiply(p, digits(x, p), digits(y, p))
            for value in range(2, word_limit):
                states = [0] * (p.carry_cap + 1)
                states[0] = 1
                remaining = value
                for coefficient in coefficients:
                    remaining, target_digit = divmod(remaining, p.base)
                    next_states = [0] * (p.carry_cap + 1)
                    for incoming, multiplicity in enumerate(states):
                        outgoing, digit = divmod(coefficient + incoming, p.base)
                        if digit == target_digit and outgoing <= p.carry_cap:
                            next_states[outgoing] += multiplicity
                    states = next_states
                assert states[0] == int(x * y == value)
                counts[value] += states[0]
                cases += 1
    for value in range(2, word_limit):
        expected = sum(value % divisor == 0 for divisor in range(1, value + 1))
        assert counts[value] == expected
    return {"candidate_pair_target_cases": cases, "complete_divisor_counts": dict(sorted(counts.items()))}


def parameter_controls() -> dict:
    cases = 0
    examples = []
    for n in range(64, 65537):
        m = (n + 1).bit_length()
        a = 1
        while (2 * a)**2 * m <= n:
            a *= 2
        t = a * a
        half = t // 2
        d = (n + half - 1) // half
        assert a >= 2
        assert t * m <= n < 4 * t * m
        assert 2 * m <= d < 8 * m + 1
        assert n <= half * d < n + half
        assert (t.bit_length() - 1) + 2 * d < 17 * m + 2
        cases += 1
        if n in (64, 128, 256, 1024, 4096, 16384, 65536):
            examples.append({"input_bits": n, "m": m, "A": a, "digit_bits": d,
                             "input_capacity_bits": half * d,
                             "log2_progression_modulus": (t.bit_length() - 1) + 2 * d})
    return {"input_lengths": cases, "examples": examples,
            "prime_search_scope": "No asymptotic-policy prime searches are benchmarked here."}


def crt_controls() -> dict:
    transform_cases = interpolation_cases = unit_cases = 0
    records = []
    for a, d in ((2, 2), (2, 32), (4, 2), (4, 32), (8, 12), (16, 4), (32, 1)):
        p, primes = make_crt_parameters(a, d)
        for prime in primes:
            assert prime % p.t == 1
            assert pow(p.omega % prime, p.t // 2, prime) == prime - 1
        for k in range(a):
            for ell in range(k + 1, a):
                assert gcd(pow(p.xi, k, p.q) - pow(p.xi, ell, p.q), p.q) == 1
                unit_cases += 1
        patterns = word_patterns(p)
        for x, y in ((patterns[-1], patterns[-1]), (patterns[-2], patterns[-1])):
            coefficients, x_hat, y_hat = grid_multiply(p, x, y)
            assert x_hat == full_dft(p, x) and y_hat == full_dft(p, y)
            assert coefficients == dft_product(p, x, y)
            assert coefficients == integer_convolution(x, y, p.t)
            carry_digits(p, coefficients)
            transform_cases += 1
        for word in patterns:
            messages = input_messages(p, word)
            for r in sorted({0, 1, a - 1}):
                expected = [word[r + a * i] for i in range(p.h)]
                boundary = {k: messages[r][k] for k in range(p.h)}
                assert recover_input(p, r, boundary) == expected
                interpolation_cases += 1
        threshold = p.t << (2 * d)
        power, count = p.t + 1, 1
        while power <= threshold:
            power *= p.t + 1
            count += 1
        auxiliary = 1
        while auxiliary < 2 * count + 3:
            auxiliary *= 3
        classes = [1 + j * p.t for j in range(auxiliary)
                   if gcd(1 + j * p.t, p.t * auxiliary) == 1]
        assert len(classes) == 2 * auxiliary // 3 >= count >= len(primes)
        records.append({
            "A": a, "digit_bits": d, "modulus": str(p.q), "residue_bits": p.f,
            "primes": primes, "omega": str(p.omega), "prime_candidates_tested": p.candidates,
            "sufficient_prime_count": count, "auxiliary_power_of_three": auxiliary,
            "available_unit_classes": len(classes),
        })

    aspect_cases = slower_policy_cases = 0
    for n in range(64, 4097):
        m = (n + 1).bit_length()
        a = 2
        while a * a * m <= n:
            t = a * a
            d = (2 * n + t - 1) // t
            assert t * d < 2 * n + t
            assert t * m <= n
            assert a * d <= (2 * n + t - 1) // a
            aspect_cases += 1
            a *= 2
    for n in range(256, 65537):
        m = (n + 1).bit_length()
        g = (m + 1).bit_length()
        a = 1
        while (2 * a)**2 * m * g <= n:
            a *= 2
        assert a >= 2
        t = a * a
        d = (2 * n + t - 1) // t
        assert t * m * g <= n < 4 * t * m * g
        assert 2 * m * g <= d < 8 * m * g + 1
        assert 2 * t * d < 4 * n + 2 * t
        slower_policy_cases += 1
    return {
        "transform_word_pairs": transform_cases,
        "recovered_input_gates": interpolation_cases,
        "unit_vandermonde_differences": unit_cases,
        "aspect_ratio_parameters": aspect_cases,
        "smaller_volume_policy_lengths": slower_policy_cases,
        "parameters": records,
    }


def actual_policy_metadata_controls() -> dict:
    records = []
    for larger_digits, lengths in (
        (False, (64, 127, 128, 256, 512, 1024, 4096, 16384, 65536, 262144, 1048576)),
        (True, (256, 1024, 4096, 65536, 1048576)),
    ):
        for n in lengths:
            m = (n + 1).bit_length()
            g = (m + 1).bit_length() if larger_digits else 1
            a = 1
            while (2 * a)**2 * m * g <= n:
                a *= 2
            assert a >= 2
            d = (2 * n + a * a - 1) // (a * a)
            p, primes = make_crt_parameters(a, d)
            assert p.q > p.t << (2 * d)
            assert p.q // primes[-1] <= p.t << (2 * d)
            assert pow(p.omega, p.t // 2, p.q) == p.q - 1
            for prime in primes:
                assert pow(p.omega % prime, p.t // 2, prime) == prime - 1
            vertices = 2
            block = a // vertices
            column_width = (6 * block * block + 2 * a - 1) * p.f
            records.append({
                "input_bits": n, "larger_digit_policy": larger_digits,
                "A": a, "digit_bits": d, "arithmetic_word_bits": p.f,
                "number_of_crt_primes": len(primes), "largest_crt_prime": primes[-1],
                "prime_candidates_tested": p.candidates,
                "all_primitive_root_components_checked": True,
                "row_carry_maximum_P_volume": 9 * p.t * p.f,
                "column_carry_P2_volume": vertices * vertices * column_width,
            })
    return {
        "cases": len(records), "samples": records,
        "scope": "Actual parameter-policy modulus/root construction only; no large full transform or factoring benchmark.",
    }


def negative_controls() -> dict:
    a, base, q = 4, 4, 17
    half = a * a // 2
    truth = integer_convolution([base - 1] * half, [base - 1] * half, 2 * half)
    assert max(truth) == 72 and max(truth) >= q
    assert [value % q for value in truth] != truth
    omega = q - 1
    xi = pow(omega, a, q)
    assert xi == 1 and len({pow(xi, k, q) for k in range(a)}) == 1
    composite_modulus, composite_root, composite_t = 15, 2, 4
    assert pow(composite_root, composite_t, composite_modulus) == 1
    assert pow(composite_root, composite_t // 2, composite_modulus) != composite_modulus - 1
    nonunit_difference = pow(composite_root, 2, composite_modulus) - 1
    assert gcd(nonunit_difference, composite_modulus) == 3
    return {
        "modulus_too_small": {"q": q, "max_integer_coefficient": max(truth),
                              "canonical_residue_at_maximum": max(truth) % q},
        "root_order_too_small": {"q": q, "omega": omega, "required_order": a * a,
                                 "actual_order": 2, "distinct_interpolation_points": 1},
        "global_composite_order_is_insufficient": {
            "modulus": composite_modulus, "root": composite_root, "global_order": composite_t,
            "nonunit_vandermonde_difference": nonunit_difference, "gcd": 3,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = {
        "status": "passed",
        "scope": "Exact finite transform, carry, interpolation, layout and codec identities; no factoring exponent claim.",
        "transforms": transform_controls(),
        "interpolation": interpolation_controls(),
        "boundaries": boundary_controls(),
        "column_carry_grouping": column_boundary_controls(),
        "factor_pairs": factor_pair_controls(),
        "parameter_policy": parameter_controls(),
        "crt_refinement": crt_controls(),
        "actual_policy_metadata": actual_policy_metadata_controls(),
        "negative_controls": negative_controls(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "transform_word_pairs": result["transforms"]["word_pairs"],
        "interpolated_gates": result["interpolation"]["recovered_input_gates"],
        "groupings": result["boundaries"]["groupings"],
        "column_groupings": result["column_carry_grouping"]["groupings"],
        "carry_control_codecs": result["boundaries"]["carry_control_codecs"],
        "factor_pair_targets": result["factor_pairs"]["candidate_pair_target_cases"],
        "parameter_lengths": result["parameter_policy"]["input_lengths"],
        "crt_transform_pairs": result["crt_refinement"]["transform_word_pairs"],
        "crt_interpolated_gates": result["crt_refinement"]["recovered_input_gates"],
        "crt_unit_differences": result["crt_refinement"]["unit_vandermonde_differences"],
        "actual_metadata_cases": result["actual_policy_metadata"]["cases"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
