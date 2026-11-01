#!/usr/bin/env python3
"""Exact finite controls for scalar recovery with a virtual rough cofactor.

This is not an implementation of the asymptotic clique evaluator or of the
Hesse--Allender--Barrington logspace arithmetic algorithms. Integer reference
calculations check the identities which those algorithms are used to execute.
The phase controls compare independently computed bucket counts, deterministic
state paths, and small transfer-matrix contractions.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import gcd, isqrt, prod
from pathlib import Path


@lru_cache(maxsize=None)
def is_prime(value: int) -> bool:
    return value >= 2 and all(value % d for d in range(2, isqrt(value) + 1))


def factor_reference(value: int) -> list[tuple[int, int]]:
    result = []
    d = 2
    while d * d <= value:
        exponent = 0
        while value % d == 0:
            value //= d
            exponent += 1
        if exponent:
            result.append((d, exponent))
        d += 1
    if value > 1:
        result.append((value, 1))
    return result


def divisor_reference(value: int) -> list[int]:
    small = [d for d in range(1, isqrt(value) + 1) if value % d == 0]
    return sorted(set(small + [value // d for d in small]))


def prime_power_slot(value: int, candidate: int) -> int:
    if not is_prime(candidate):
        return 1
    power = candidate
    slot = 1
    while power <= value and value % power == 0:
        slot = power
        power *= candidate
    return slot


def factor_via_divisor_products(value: int) -> list[tuple[int, int]]:
    """Finite reference for the descriptor-side prime and power scans."""
    divisors = divisor_reference(value)
    answer = []
    for d in divisors:
        if d == 1:
            continue
        if any(a > 1 and c > 1 and a * c == d for a in divisors for c in divisors):
            continue
        exponent, power = 1, d
        while any(candidate == power * d for candidate in divisors):
            power *= d
            exponent += 1
        answer.append((d, exponent))
    return answer


def ceil_log2(value: int) -> int:
    assert value >= 1
    return (value - 1).bit_length()


def rough_controls() -> dict:
    cases = 0
    examples = []
    for value in range(2, 512):
        n = value.bit_length()
        for s in range(1, min(n, 6) + 1):
            threshold = 1 << s
            slots = [prime_power_slot(value, candidate) for candidate in range(2, threshold + 1)]
            smooth = prod(slots)
            assert value % smooth == 0
            rough = value // smooth
            assert gcd(smooth, rough) == 1
            assert all(p > threshold for p, _ in factor_reference(rough))
            assert all(p <= threshold for p, _ in factor_reference(smooth))
            q = (n + s - 1) // s
            divisors = divisor_reference(rough)
            assert len(divisors) <= 1 << q
            omega = sum(e for _, e in factor_reference(rough))
            assert rough == 1 or s * omega < n
            small_factors = [
                (p, e) for p, e in factor_reference(value) if p <= threshold
            ]
            recovered = small_factors + factor_via_divisor_products(rough)
            assert recovered == factor_reference(value)
            assert prod(p**e for p, e in recovered) == value
            cases += 1
            if (value, s) in {(12, 2), (45, 1), (210, 2), (385, 2), (511, 2)}:
                examples.append({
                    "N": value, "s": s, "threshold": threshold,
                    "smooth": smooth, "rough": rough, "rough_divisors": len(divisors),
                    "q": q, "factorization": recovered,
                })

    value, s = 1 << 20, 10
    false_q = (value.bit_length() + s - 1) // s
    assert len(divisor_reference(value)) > 1 << false_q
    wrong_slots = prod(
        max([1] + [d**e for e in range(1, (12).bit_length() + 1) if 12 % (d**e) == 0])
        for d in range(2, 5)
    )
    assert wrong_slots == 48 and 12 % wrong_slots != 0
    return {
        "cases": cases,
        "examples": examples,
        "nonrough_support_negative_control": {
            "N": value, "s": s, "tau": len(divisor_reference(value)),
            "incorrect_unconditional_bound": 1 << false_q,
        },
        "composite_slot_negative_control": {"N": 12, "B": 4, "incorrect_product": wrong_slots},
    }


def adaptive_threshold_controls() -> dict:
    cases = 0
    examples = []
    inputs = list(range(2, 1025)) + [
        2**128, 3**80, 2**128 * 1009,
        2**80 * 3**40 * 5**20 * 1009 * 1013,
    ]
    for value in inputs:
        n = value.bit_length()
        factors = factor_reference(value)
        rough_parts = [
            prod(prime**exponent for prime, exponent in factors if prime > 1 << s)
            for s in range(1, n + 1)
        ]
        lengths = [0 if rough == 1 else rough.bit_length() for rough in rough_parts]
        assert all(a >= b for a, b in zip(lengths, lengths[1:]))
        stop = next(s for s, length in enumerate(lengths, 1) if length <= s * s)
        assert stop <= isqrt(n - 1) + 1
        for t, length in enumerate(lengths, 1):
            assert stop * t <= t * t + length
        assert lengths[stop - 1] <= stop * stop
        emitted = []
        for s in range(1, stop + 1):
            emitted.extend((prime, exponent) for prime, exponent in factors
                           if 1 << (s - 1) < prime <= 1 << s)
        emitted.extend(factor_reference(rough_parts[stop - 1]))
        assert emitted == factors
        if value in (12, 231, 511, 1024) or n > 20:
            best_t = min(range(1, n + 1), key=lambda t: Fraction(t * t + lengths[t - 1], t))
            examples.append({
                "N": str(value), "input_bits": n, "stopping_threshold": stop,
                "rough_bits": lengths[stop - 1], "rough_part": str(rough_parts[stop - 1]),
                "optimal_visible_cost": {
                    "numerator": best_t * best_t + lengths[best_t - 1], "denominator": best_t,
                },
            })
        cases += 1
    return {"inputs": cases, "examples": examples}


def digits(value: int, base: int, length: int) -> tuple[int, ...]:
    output = []
    for _ in range(length):
        value, digit = divmod(value, base)
        output.append(digit)
    assert value == 0
    return tuple(output)


def trim(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def polynomial_remainder(poly: tuple[int, ...], divisor: tuple[int, ...], p: int) -> tuple[int, ...]:
    remainder = trim(list(poly))
    assert divisor[-1] == 1
    while len(remainder) >= len(divisor) and any(remainder):
        shift = len(remainder) - len(divisor)
        leading = remainder[-1]
        for j, coefficient in enumerate(divisor):
            remainder[shift + j] = (remainder[shift + j] - leading * coefficient) % p
        trim(remainder)
    return tuple(remainder)


@lru_cache(maxsize=None)
def irreducible_polynomial(p: int, degree: int) -> tuple[int, ...]:
    if degree == 1:
        return (0, 1)
    for low_code in range(p**degree):
        low = digits(low_code, p, degree)
        if low[0] == 0:
            continue
        candidate = low + (1,)
        reducible = False
        for factor_degree in range(1, degree // 2 + 1):
            for factor_code in range(p**factor_degree):
                factor = digits(factor_code, p, factor_degree) + (1,)
                if not any(polynomial_remainder(candidate, factor, p)):
                    reducible = True
                    break
            if reducible:
                break
        if not reducible:
            return candidate
    raise AssertionError("finite field polynomial not found")


class ControlField:
    def __init__(self, p: int, degree: int):
        self.p = p
        self.degree = degree
        self.order = p**degree
        self.modulus = irreducible_polynomial(p, degree)
        self.coordinates = [digits(value, p, degree) for value in range(self.order)]

    def encode(self, coefficients: tuple[int, ...] | list[int]) -> int:
        return sum((coefficient % self.p) * self.p**j for j, coefficient in enumerate(coefficients))

    def add(self, a: int, b: int) -> int:
        return self.encode([(x + y) % self.p for x, y in zip(self.coordinates[a], self.coordinates[b])])

    def multiply(self, a: int, b: int) -> int:
        coefficients = [0] * (2 * self.degree - 1)
        for i, ai in enumerate(self.coordinates[a]):
            for j, bj in enumerate(self.coordinates[b]):
                coefficients[i + j] = (coefficients[i + j] + ai * bj) % self.p
        remainder = polynomial_remainder(tuple(coefficients), self.modulus, self.p)
        return self.encode(remainder)

    def powers(self, value: int, length: int) -> list[int]:
        result, current = [], 1
        for _ in range(length):
            result.append(current)
            current = self.multiply(current, value)
        return result

    def dot(self, a: int, b: int) -> int:
        return sum(x * y for x, y in zip(self.coordinates[a], self.coordinates[b])) % self.p

    def hash_word(self, value: int, powers: list[int]) -> int:
        answer = 0
        for j, power in enumerate(powers):
            if (value >> j) & 1:
                answer = self.add(answer, power)
        return answer


def scalar_weight(coefficient_prime: int, phase: int) -> int:
    if coefficient_prime == 2:
        return (1, 0, 1)[phase % 3]
    return 1 if phase % 2 == 0 else coefficient_prime - 1


def phase_matrix_chain(a: int, phases: list[int], offset: int, coefficient_prime: int) -> int:
    """Independent finite-state contraction, including an invalid fourth state."""
    ell = 3 if coefficient_prime == 2 else 2
    width = 4 if ell == 3 else 2
    vector = [1] + [0] * (width - 1)
    for j, coefficient in enumerate(phases):
        bit = (a >> j) & 1
        matrix = [
            [int(q < ell and r < ell and r == (q + coefficient * bit) % ell)
             for r in range(width)]
            for q in range(width)
        ]
        vector = [sum(vector[q] * matrix[q][r] for q in range(width)) % coefficient_prime
                  for r in range(width)]
    if ell == 3:
        assert vector[3] == 0
    assert sum(vector) == 1
    return sum(vector[q] * scalar_weight(coefficient_prime, q - offset)
               for q in range(ell)) % coefficient_prime


def scalar_phase_vector(
    n: int, coefficient_prime: int, pairs: list[tuple[int, int]],
    phases: list[int], offset: int, check_chain: bool,
) -> tuple[list[int], int]:
    ell = 3 if coefficient_prime == 2 else 2
    values = [0] * (2 * n + 1)
    chain_controls = 0
    for a, b in pairs:
        state = 0
        for j, coefficient in enumerate(phases):
            state = (state + coefficient * ((a >> j) & 1)) % ell
        weight = scalar_weight(coefficient_prime, state - offset)
        if check_chain:
            assert phase_matrix_chain(a, phases, offset, coefficient_prime) == weight
            chain_controls += 1
        values[0] += weight
        for j in range(n):
            values[1 + j] += weight * ((a >> j) & 1)
            values[1 + n + j] += weight * ((b >> j) & 1)
    return [value % coefficient_prime for value in values], chain_controls


def true_bucket_vector(
    n: int, coefficient_prime: int, pairs: list[tuple[int, int]],
    field: ControlField, powers: list[int], bucket: int,
) -> list[int]:
    values = [0] * (2 * n + 1)
    for a, b in pairs:
        if field.hash_word(a, powers) != bucket:
            continue
        values[0] += 1
        for j in range(n):
            values[1 + j] += (a >> j) & 1
            values[1 + n + j] += (b >> j) & 1
    return [value % coefficient_prime for value in values]


def phase_controls() -> dict:
    results = []
    projection_cases = chain_cases = isolation_cases = 0
    for coefficient_prime in (2, 3, 5, 7):
        ell = 3 if coefficient_prime == 2 else 2
        for value in (6, 12, 15, 35, 45, 49, 77, 105, 225, 231):
            n = value.bit_length()
            divisors = divisor_reference(value)
            pairs = [(d, value // d) for d in divisors]
            degree = 1
            while ell**degree <= n * (len(divisors) - 1):
                degree += 1
            field = ControlField(ell, degree)
            powers_by_t = [field.powers(t, n) for t in range(field.order)]
            first_isolator = {}
            for t, powers in enumerate(powers_by_t):
                buckets: dict[int, list[int]] = defaultdict(list)
                for d in divisors:
                    buckets[field.hash_word(d, powers)].append(d)
                for bucket, members in buckets.items():
                    if len(members) == 1:
                        first_isolator.setdefault(members[0], (t, bucket))
            assert set(first_isolator) == set(divisors)
            isolation_cases += len(divisors)
            selected_t = {0, 1, first_isolator[divisors[1]][0], first_isolator[value][0]}
            local_projection_cases = 0
            for t in sorted(selected_t):
                powers = powers_by_t[t]
                selected_buckets = {0, field.order - 1}
                selected_buckets.update(field.hash_word(d, powers) for d in (1, divisors[1], value))
                for bucket in sorted(selected_buckets):
                    transformed = [0] * (2 * n + 1)
                    for v in range(field.order):
                        phases = [field.dot(v, power) for power in powers]
                        offset = field.dot(v, bucket)
                        scalar, checked = scalar_phase_vector(
                            n, coefficient_prime, pairs, phases, offset,
                            v in {0, 1, field.order - 1},
                        )
                        chain_cases += checked
                        transformed = [(x + y) % coefficient_prime for x, y in zip(transformed, scalar)]
                    inverse = pow(field.order, -1, coefficient_prime)
                    transformed = [inverse * x % coefficient_prime for x in transformed]
                    truth = true_bucket_vector(n, coefficient_prime, pairs, field, powers, bucket)
                    assert transformed == truth
                    projection_cases += 1
                    local_projection_cases += 1
            results.append({
                "coefficient_prime": coefficient_prime,
                "hash_characteristic": ell,
                "N": value, "divisors": len(divisors), "control_order": field.order,
                "control_polynomial": list(field.modulus),
                "projection_cases": local_projection_cases,
                "all_divisors_have_singleton_handles": True,
            })

    value, coefficient_prime = 12, 2
    field = ControlField(3, 3)
    n = value.bit_length()
    pairs = [(d, value // d) for d in divisor_reference(value)]
    bucket = true_bucket_vector(n, coefficient_prime, pairs, field, field.powers(1, n), 1)
    assert bucket[0] == 1 and all(bit in (0, 1) for bit in bucket[1:])
    a = sum(bit << j for j, bit in enumerate(bucket[1:1 + n]))
    b = sum(bit << j for j, bit in enumerate(bucket[1 + n:]))
    assert (a, b) == (7, 9) and a * b != value

    return {
        "projection_cases": projection_cases,
        "independent_transfer_matrix_cases": chain_cases,
        "isolated_divisors": isolation_cases,
        "cases": results,
        "colliding_scalar_one_negative_control": {
            "N": value, "coefficient_prime": coefficient_prime, "t": 1, "bucket": 1,
            "actual_bucket": [1, 2, 4], "candidate_a": a, "candidate_b": b,
            "product": a * b, "multiplication_rejects": True,
        },
    }


def composite_ring_controls() -> dict:
    records = []
    cases = 0
    for modulus in (4, 6, 8, 9, 10, 12, 15):
        ell = next(prime for prime in (2, 3, 5, 7) if gcd(prime, modulus) == 1)
        field = ControlField(ell, 2)
        inverse = pow(field.order, -1, modulus)
        for x in range(field.order):
            total = 0
            for v in range(field.order):
                phase = field.dot(v, x)
                total += int(phase == 0) - int(phase == ell - 1)
            assert (total * inverse) % modulus == int(x == 0)
            cases += 1
        records.append({
            "coefficient_modulus": modulus,
            "hash_characteristic": ell,
            "control_order": field.order,
            "all_linear_forms_checked": True,
        })
    return {"linear_form_cases": cases, "rings": records}


def coefficient_extension_controls() -> dict:
    records = []
    forms = 0
    for p, degree in ((2, 2), (2, 3), (3, 2), (5, 2), (7, 2)):
        coefficients = ControlField(p, degree)
        ell = 3 if p == 2 else 2
        control = ControlField(ell, 2)
        inverse = pow(control.order, -1, p)
        for a in range(p):
            for b in range(p):
                assert coefficients.add(a, b) == (a + b) % p
                assert coefficients.multiply(a, b) == (a * b) % p
        for x in range(control.order):
            accumulator = 0
            for v in range(control.order):
                weight = scalar_weight(p, control.dot(v, x))
                accumulator = coefficients.add(accumulator, weight)
            assert coefficients.multiply(inverse, accumulator) == int(x == 0)
            forms += 1
        cubic_character_roots = [
            x for x in range(coefficients.order)
            if coefficients.add(coefficients.add(coefficients.multiply(x, x), x), 1) == 0
        ]
        if (p, degree) == (2, 3):
            assert not cubic_character_roots
        records.append({
            "characteristic": p, "degree": degree, "order": coefficients.order,
            "control_order": control.order,
            "roots_of_X2_plus_X_plus_1": cubic_character_roots,
            "scalar_filter_works_without_adjoining_roots": True,
        })
    return {"linear_form_cases": forms, "fields": records}


def grouped_phase_controls() -> dict:
    """Check the added group path, without simulating multiplier gates."""
    volume_cases = path_cases = balanced_width_cases = 0
    largest_ratio = (0, 1)
    largest_example = None
    for n in range(2, 16385):
        source_scale = 1
        while n > 4 * source_scale**4:
            source_scale *= 2
        gates = 2 * source_scale
        field_bits = 8 * source_scale**2 + 1
        if n >= 1024:
            threshold = isqrt(4 * n - 1) + 1
            control_bits_upper = (n + threshold - 1) // threshold + ceil_log2(n) + 1
            assert control_bits_upper <= 2 * source_scale**2
            balanced_width_cases += 1
        vertices = 2
        while vertices <= gates:
            block = gates // vertices
            width = (8 * block**2 + 1) * field_bits
            mass = vertices**2 * width
            assert mass <= 1296 * n
            assert gates**2 * field_bits <= 144 * n
            volume_cases += 1
            if mass * largest_ratio[1] > largest_ratio[0] * n:
                largest_ratio = mass, n
                largest_example = {"n": n, "vertices": vertices, "width": width, "mass": mass}
            vertices *= 2

    # This identity holds for any ownership partition. The cyclic assignment
    # here intentionally need not match the multiplier's particular gate codec.
    for coefficient_prime in (2, 3, 5, 7):
        ell = 3 if coefficient_prime == 2 else 2
        for vertices in (2, 4, 8):
            n = 2 * vertices + 3
            for word in (0, 1, (1 << n) - 1, sum(1 << j for j in range(0, n, 2))):
                phases = [(j * j + 2 * j + 1) % ell for j in range(n)]
                local = [0] * vertices
                for j, coefficient in enumerate(phases):
                    owner = (3 * j + 1) % vertices
                    local[owner] = (local[owner] + coefficient * ((word >> j) & 1)) % ell
                for offset in range(ell):
                    paths = value = 0
                    for labels in product(range(4), repeat=vertices - 1):
                        if any(label >= ell for label in labels):
                            continue
                        states = (0,) + labels
                        if any(states[a + 1] != (states[a] + local[a]) % ell
                               for a in range(vertices - 1)):
                            continue
                        paths += 1
                        value += scalar_weight(coefficient_prime, states[-1] + local[-1] - offset)
                    assert paths == 1
                    assert value % coefficient_prime == scalar_weight(
                        coefficient_prime, sum(local) - offset
                    )
                    path_cases += 1

    # If the spare two-bit code is treated as a second encoding of zero,
    # one valid original setting acquires multiplicity two in characteristic two.
    unchecked_encodings = [label for label in range(4) if label % 3 == 0]
    assert unchecked_encodings == [0, 3] and len(unchecked_encodings) % 2 == 0
    return {
        "volume_cases": volume_cases,
        "balanced_hash_width_cases": balanced_width_cases,
        "group_path_contractions": path_cases,
        "largest_checked_mass_ratio": {
            "numerator": largest_ratio[0], "denominator": largest_ratio[1],
            "example": largest_example,
        },
        "invalid_encoding_negative_control": {
            "coefficient_prime": 2, "valid_encoding": 0, "forbidden_encoding": 3,
            "incorrect_modular_path_count": 0, "correct_path_count": 1,
        },
        "scope": "Additional phase channels and exact parameter bounds only; no multiplier table reimplementation.",
    }


def carry_codec_controls() -> dict:
    transition_cases = codec_cases = 0
    for d in (2, 8, 32, 128, 512):
        base = 1 << d
        modulus = base**4 + 1
        cap = base**3 + base**2 + base + 1
        carry_bits = 3 * d + 1
        total_bits = 4 * d + 1
        assert cap < 1 << carry_bits
        assert (modulus - 1 + cap) // base <= cap
        if d == 2:
            carries = list(range(cap + 1))
            states = list(range(1 << d))
            for carry in carries:
                for coefficient in range(modulus):
                    assert (carry + coefficient) // base <= cap
                    transition_cases += 1
        else:
            carries = [0, 1, cap // 2, cap]
            states = [0, 1, (1 << d) - 1]
        for carry in carries:
            for state in states:
                packed = carry + (state << carry_bits)
                assert packed < 1 << total_bits
                assert (packed & ((1 << carry_bits) - 1)) == carry
                assert packed >> carry_bits == state
                codec_cases += 1
    return {
        "carry_transitions": transition_cases,
        "packed_carry_control_values": codec_cases,
        "control_bits_fit_without_increasing_original_width": True,
    }


def direct_hash_path_controls() -> dict:
    bucket_cases = path_extensions = 0
    for degree in (1, 2, 3):
        field = ControlField(2, degree)
        for value in (6, 12, 35, 77):
            n = value.bit_length()
            pairs = [(divisor, value // divisor) for divisor in divisor_reference(value)]
            for vertices in (2, 4):
                for t in sorted({0, 1, field.order - 1}):
                    powers = field.powers(t, n)
                    vectors = [[0] * (2 * n + 1) for _ in range(field.order)]
                    for a, b in pairs:
                        local = [0] * vertices
                        for j, power in enumerate(powers):
                            if (a >> j) & 1:
                                owner = (3 * j + 1) % vertices
                                local[owner] = field.add(local[owner], power)
                        valid_extensions = 0
                        for labels in product(range(field.order), repeat=vertices - 1):
                            states = (0,) + labels
                            if any(states[i + 1] != field.add(states[i], local[i])
                                   for i in range(vertices - 1)):
                                continue
                            valid_extensions += 1
                            bucket = field.add(states[-1], local[-1])
                            vectors[bucket][0] += 1
                            for j in range(n):
                                vectors[bucket][1 + j] += (a >> j) & 1
                                vectors[bucket][1 + n + j] += (b >> j) & 1
                        assert valid_extensions == 1
                        path_extensions += 1
                    for coefficient_prime in (2, 3, 5):
                        for bucket in range(field.order):
                            actual = [entry % coefficient_prime for entry in vectors[bucket]]
                            expected = true_bucket_vector(n, coefficient_prime, pairs, field, powers, bucket)
                            assert actual == expected
                            bucket_cases += 1
    return {
        "scalar_bucket_vectors": bucket_cases,
        "unique_control_extensions": path_extensions,
        "scope": "Direct group hash paths and bit pins; no original multiplier gates simulated.",
    }


def output_controls() -> dict:
    primes = []
    candidate = 2
    while len(primes) < 128:
        if is_prime(candidate):
            primes.append(candidate)
        candidate += 1
    primorial = 1
    records = []
    for k, prime in enumerate(primes, 1):
        primorial *= prime
        n = primorial.bit_length()
        prime_supply_bound = 8 * k * ceil_log2(k + 2)
        assert prime <= prime_supply_bound
        assert k <= n <= 1 + k * ceil_log2(prime_supply_bound)
        assert len(set(primes[:k])) == k
        if k in (1, 2, 4, 8, 16, 32, 64, 128):
            records.append({
                "number_of_prime_factors": k,
                "input_bits": n,
                "largest_prime": prime,
                "prime_supply_bound": prime_supply_bound,
                "divisor_output_count": str(1 << k),
                "log2_divisor_output_count": k,
            })
    return {"primorial_cases": len(primes), "examples": records}


def enumeration_controls() -> dict:
    cases = comparisons = interval_cases = 0
    nonmonotone_example = None
    for value in range(2, 512):
        factors = factor_reference(value)
        count = prod(exponent + 1 for _, exponent in factors)

        def divisor_at(index: int) -> int:
            assert 0 <= index < count
            answer = 1
            remaining = index
            for prime, exponent in factors:
                remaining, digit = divmod(remaining, exponent + 1)
                answer *= prime**digit
            assert remaining == 0
            return answer

        mixed_values = [divisor_at(index) for index in range(count)]
        truth = divisor_reference(value)
        assert sorted(mixed_values) == truth
        if nonmonotone_example is None and mixed_values != truth:
            nonmonotone_example = {"N": value, "mixed_order": mixed_values, "sorted_order": truth}

        previous = None
        emitted = []
        for _ in range(count):
            best = None
            for index in range(count):
                candidate = divisor_at(index)
                if previous is not None and candidate <= divisor_at(previous):
                    comparisons += 1
                    continue
                if best is None or candidate < divisor_at(best):
                    best = index
                comparisons += 1
            assert best is not None
            emitted.append(divisor_at(best))
            previous = best
        assert emitted == truth
        for lower, upper in ((0, value), (2, value - 1), (value // 3, 2 * value // 3)):
            assert sum(lower <= d <= upper for d in emitted) == sum(
                value % d == 0 for d in range(max(1, lower), upper + 1)
            )
            interval_cases += 1
        cases += 1
    assert nonmonotone_example is not None
    return {
        "inputs": cases, "successor_scan_comparisons": comparisons,
        "interval_cases": interval_cases,
        "mixed_radix_is_not_sorted_control": nonmonotone_example,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = {
        "status": "passed",
        "scope": "Exact finite identities; not a space benchmark or a faster factoring evaluator.",
        "rough_preprocessing": rough_controls(),
        "adaptive_thresholds": adaptive_threshold_controls(),
        "scalar_phase_recovery": phase_controls(),
        "composite_coefficient_rings": composite_ring_controls(),
        "coefficient_extensions": coefficient_extension_controls(),
        "explicit_group_phase_path": grouped_phase_controls(),
        "spare_carry_codec": carry_codec_controls(),
        "direct_scalar_hash_paths": direct_hash_path_controls(),
        "output_size": output_controls(),
        "output_sensitive_enumeration": enumeration_controls(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print(json.dumps({
        "status": result["status"],
        "rough_cases": result["rough_preprocessing"]["cases"],
        "adaptive_threshold_inputs": result["adaptive_thresholds"]["inputs"],
        "phase_projection_cases": result["scalar_phase_recovery"]["projection_cases"],
        "matrix_chain_cases": result["scalar_phase_recovery"]["independent_transfer_matrix_cases"],
        "isolated_divisors": result["scalar_phase_recovery"]["isolated_divisors"],
        "composite_ring_forms": result["composite_coefficient_rings"]["linear_form_cases"],
        "extension_field_forms": result["coefficient_extensions"]["linear_form_cases"],
        "group_phase_cases": result["explicit_group_phase_path"]["group_path_contractions"],
        "explicit_volume_cases": result["explicit_group_phase_path"]["volume_cases"],
        "balanced_hash_width_cases": result["explicit_group_phase_path"]["balanced_hash_width_cases"],
        "carry_transition_cases": result["spare_carry_codec"]["carry_transitions"],
        "direct_bucket_vectors": result["direct_scalar_hash_paths"]["scalar_bucket_vectors"],
        "primorial_cases": result["output_size"]["primorial_cases"],
        "sorted_enumeration_inputs": result["output_sensitive_enumeration"]["inputs"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
