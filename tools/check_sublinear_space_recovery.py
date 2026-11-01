#!/usr/bin/env python3
"""Exact controls for paired tags and descriptor-only divisor recovery.

The finite fixture supplies the scalar and bit outputs of an isolated
bucket. Its stored witness/bucket tables are verification data, not part
of the descriptor algorithms and not an implementation of the large
clique evaluator. The control algorithms below retain handles and small
counters; they never decode a candidate into an n-bit integer.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import random

from check_linear_space_recovery import TernaryField, PHASES, f4_mul, ring_mul


@dataclass
class Counts:
    scalar_queries: int = 0
    bit_queries: int = 0
    product_checks: int = 0
    comparisons: int = 0
    max_carry: int = 0
    max_column_sum: int = 0


def product_matches(width, left_bit, right_bit, target_bit, counts=None):
    """Check an exact product using O(log(width+2)) control bits.

    Factors have width bits. The target callback is read through 2*width
    positions; a width-bit target returns zero in its upper positions.
    Callback workspace is separately charged. No product word is built.
    """
    if counts is not None:
        counts.product_checks += 1
    carry = 0
    for column in range(2 * width):
        total = carry
        for i in range(max(0, column - width + 1), min(width, column + 1)):
            a, b = left_bit(i), right_bit(column - i)
            assert a in (0, 1) and b in (0, 1)
            total += a * b
        assert carry <= width and total <= 2 * width
        if counts is not None:
            counts.max_carry = max(counts.max_carry, carry)
            counts.max_column_sum = max(counts.max_column_sum, total)
        target = target_bit(column)
        assert target in (0, 1)
        if total & 1 != target:
            return False
        carry = total >> 1
    assert carry <= width
    return carry == 0


def paired_query(pairs, width, variable):
    """Independent literal weighted-pair polynomial over the 16-element ring."""
    scalar = tag = 0
    for first, second in pairs:
        term = (1, 0)
        for role, witness in enumerate((first, second)):
            for bit in range(width):
                if (witness >> bit) & 1:
                    term = ring_mul(term, variable(role, bit))
        scalar ^= term[0]
        tag ^= term[1]
    return scalar, tag


def projected_pair_bit(pairs, width, point, bucket, role, bit, field):
    scalar = tag = 0
    for character in range(field.size):
        def variable(which, position):
            phase = (PHASES[field.dot(character, field.power(point, position))]
                     if which == 0 else 1)
            return phase, phase if (which, position) == (role, bit) else 0

        a, b = paired_query(pairs, width, variable)
        factor = PHASES[-field.dot(character, bucket) % 3]
        scalar ^= f4_mul(factor, a)
        tag ^= f4_mul(factor, b)
    return scalar, tag


class FiniteBucketFixture:
    """A finite truth-table oracle used only for verification.

    The production evaluator replaces scalar()/bit(). The recovery code
    is not allowed to access pairs or values. In particular, its handles
    represent oracle outputs rather than storing those outputs.
    """

    def __init__(self, number, field_degree=None):
        self.number = number  # Read-only input fixture.
        self.width = number.bit_length()
        self.pairs = tuple((d, number // d) for d in range(1, number + 1)
                           if number % d == 0)
        # Exact support size keeps the fixtures affordable. The theorem
        # replaces it by the uniform bound checked in divisor_bound_checks().
        degree = 1
        while 3 ** degree <= self.width * (len(self.pairs) - 1):
            degree += 1
        self.field = TernaryField(degree if field_degree is None else field_degree)
        self.total_handles = self.field.size ** 2
        # Full buckets and membership counts are independent test fixtures.
        self.values = [[0, 0, 0] for _ in range(self.total_handles)]
        self.members = [[] for _ in range(self.total_handles)]
        for point in range(self.field.size):
            for first, second in self.pairs:
                bucket = self.field.hash(first, point, self.width)
                handle = point * self.field.size + bucket
                value = self.values[handle]
                value[0] ^= 1
                value[1] ^= first
                value[2] ^= second
                self.members[handle].append(first)
        self.counts = Counts()

    def scalar(self, handle):
        self.counts.scalar_queries += 1
        return self.values[handle][0]

    def bit(self, handle, role, position):
        self.counts.bit_queries += 1
        if position < 0 or position >= self.width:
            return 0
        return (self.values[handle][role + 1] >> position) & 1

    def input_bit(self, position):
        return (self.number >> position) & 1

    def independently_isolated(self):
        return {members[0] for members in self.members if len(members) == 1}


class DescriptorRecovery:
    """Recovery control that sees only one-bit oracle answers and handles."""

    def __init__(self, oracle):
        self.oracle = oracle
        self.width = oracle.width

    def handles(self):
        yield from range(self.oracle.total_handles)

    def valid(self, handle):
        if self.oracle.scalar(handle) != 1:
            return False
        return product_matches(
            self.width,
            lambda i: self.oracle.bit(handle, 0, i),
            lambda i: self.oracle.bit(handle, 1, i),
            self.oracle.input_bit,
            self.oracle.counts,
        )

    def nonunit(self, handle):
        # A valid factor is positive. It exceeds one exactly when a higher bit is set.
        for bit in range(1, self.width):
            if self.oracle.bit(handle, 0, bit):
                return True
        return False

    def compare(self, first, second):
        self.oracle.counts.comparisons += 1
        for bit in range(self.width - 1, -1, -1):
            a = self.oracle.bit(first, 0, bit)
            b = self.oracle.bit(second, 0, bit)
            if a != b:
                return 1 if a else -1
        return 0

    def next_after(self, previous=None, exclude_one=False):
        best = None
        for handle in self.handles():
            if not self.valid(handle):
                continue
            if exclude_one and not self.nonunit(handle):
                continue
            if previous is not None and self.compare(handle, previous) <= 0:
                continue
            if best is None or self.compare(handle, best) < 0:
                best = handle
        return best

    def divisors(self):
        previous = None
        while True:
            following = self.next_after(previous)
            if following is None:
                return
            yield following
            previous = following

    def least_factor(self):
        return self.next_after(exclude_one=True)

    def prefix_count(self, free_bits, prefix_bit):
        count = 0
        for handle in self.divisors():
            matches = True
            for position in range(free_bits, self.width):
                if self.oracle.bit(handle, 0, position) != prefix_bit(position - free_bits):
                    matches = False
                    break
            count += matches
        return count

    def multiply_equals(self, first, second, target):
        return product_matches(
            self.width,
            lambda i: self.oracle.bit(first, 0, i),
            lambda i: self.oracle.bit(second, 0, i),
            lambda i: self.oracle.bit(target, 0, i),
            self.oracle.counts,
        )

    def prime(self, target):
        """Target is already a validated divisor; no primality oracle is used."""
        if not self.nonunit(target):
            return False
        for first in self.handles():
            if not self.valid(first) or not self.nonunit(first):
                continue
            if self.compare(first, target) >= 0:
                continue
            for second in self.handles():
                if not self.valid(second) or not self.nonunit(second):
                    continue
                if self.compare(second, target) >= 0:
                    continue
                if self.multiply_equals(first, second, target):
                    return False
        return True

    def prime_power_exponent(self, prime):
        power, exponent = prime, 1
        while True:
            following = None
            for handle in self.handles():
                if self.valid(handle) and self.multiply_equals(power, prime, handle):
                    following = handle
                    break
            if following is None:
                return exponent
            power = following
            exponent += 1
            assert exponent <= self.width

    def factorization(self):
        for handle in self.divisors():
            if self.prime(handle):
                yield handle, self.prime_power_exponent(handle)


def decode_for_test(oracle, handle, role=0):
    """Only the verifier decodes an entire output; the control code does not."""
    return sum(oracle.bit(handle, role, i) << i for i in range(oracle.width))


def elementary_factors(number):
    result, divisor = [], 2
    while divisor * divisor <= number:
        exponent = 0
        while number % divisor == 0:
            exponent += 1
            number //= divisor
        if exponent:
            result.append((divisor, exponent))
        divisor += 1
    if number > 1:
        result.append((number, 1))
    return result


def divisor_bound_checks():
    """Exact integer controls for the uniform bucket-description length."""
    inequality_checks = 0
    for length in range(2, 1025):
        # Square both positive sides of L(L+1) <= 6*2^(L/2).
        assert (length * (length + 1)) ** 2 <= 36 * (1 << length)
        if length >= 5:
            assert (length + 2) * (length + 3) <= 2 * length * (length + 1)
        inequality_checks += 1

    def check(number, divisor_count):
        width = number.bit_length()
        if width < 4:
            return
        length = width.bit_length() - 1
        exponent = min(width, (9 * width + length - 1) // length)
        assert divisor_count <= 1 << exponent

    limit = 8192
    divisor_counts = [0] * (limit + 1)  # Independent verification fixture.
    for divisor in range(1, limit + 1):
        for multiple in range(divisor, limit + 1, divisor):
            divisor_counts[multiple] += 1
    cases = 0
    for number in range(2, limit + 1):
        check(number, divisor_counts[number])
        cases += 1

    rng = random.Random(20261031)
    primes = (2, 3, 5, 7, 11, 13, 17, 19)
    maximum_width = 0
    for _ in range(96):
        number, divisor_count = 1, 1
        for prime in primes:
            exponent = rng.randrange(2048)
            number *= prime ** exponent
            divisor_count *= exponent + 1
        check(number, divisor_count)
        maximum_width = max(maximum_width, number.bit_length())
        cases += 1
    return {"inequality_checks": inequality_checks,
            "divisor_bound_checks": cases,
            "largest_constructed_input_bits": maximum_width}


def verify():
    bound_checks = divisor_bound_checks()
    product_counts = Counts()
    product_cases = 0
    for width in range(1, 5):
        limit = 1 << width
        for first in range(limit):
            for second in range(limit):
                for target in range(limit):
                    answer = product_matches(
                        width, lambda i: (first >> i) & 1,
                        lambda i: (second >> i) & 1,
                        lambda i: (target >> i) & 1, product_counts,
                    )
                    assert answer == (first * second == target)
                    product_cases += 1
                # Also check full 2*width-bit outputs, including final carries.
                assert product_matches(
                    width, lambda i: (first >> i) & 1,
                    lambda i: (second >> i) & 1,
                    lambda i: ((first * second) >> i) & 1, product_counts,
                )
                product_cases += 1
    rng = random.Random(20261030)
    for width in (8, 16, 32, 64):
        for _ in range(24):
            first, second = rng.randrange(1 << width), rng.randrange(1 << width)
            target = first * second
            assert product_matches(width, lambda i: (first >> i) & 1,
                                   lambda i: (second >> i) & 1,
                                   lambda i: (target >> i) & 1, product_counts)
            wrong = target ^ (1 << rng.randrange(2 * width))
            assert not product_matches(width, lambda i: (first >> i) & 1,
                                       lambda i: (second >> i) & 1,
                                       lambda i: (wrong >> i) & 1, product_counts)
            product_cases += 2
    for width in (1, 2, 3, 8, 16, 32, 64, 128):
        first = second = (1 << width) - 1
        target = first * second
        assert product_matches(width, lambda i: (first >> i) & 1,
                               lambda i: (second >> i) & 1,
                               lambda i: (target >> i) & 1, product_counts)
        product_cases += 1

    projection_cases = 0
    for number, degree in ((6, 1), (6, 2), (12, 2), (15, 2)):
        fixture = FiniteBucketFixture(number, degree)
        for handle, expected in enumerate(fixture.values):
            point, bucket = divmod(handle, fixture.field.size)
            for role in (0, 1):
                for bit in range(fixture.width):
                    actual = projected_pair_bit(fixture.pairs, fixture.width, point,
                                                bucket, role, bit, fixture.field)
                    assert actual == (expected[0], (expected[role + 1] >> bit) & 1)
                    projection_cases += 1

    reports = []
    prefix_cases = duplicate_handles = rejected_odd_buckets = 0
    for number in (2, 3, 4, 6, 9, 12, 15, 16, 18, 25, 49, 77):
        fixture = FiniteBucketFixture(number)
        recovery = DescriptorRecovery(fixture)
        expected_divisors = [first for first, _ in fixture.pairs]
        assert fixture.independently_isolated() == set(expected_divisors)
        valid = []  # Independent instrumentation, not retained by recovery.
        for handle in recovery.handles():
            accepted = recovery.valid(handle)
            if accepted:
                value = decode_for_test(fixture, handle)
                assert number % value == 0 and value > 0
                valid.append(value)
            elif fixture.values[handle][0] == 1:
                rejected_odd_buckets += 1
        duplicate_handles += len(valid) - len(set(valid))
        generated = [decode_for_test(fixture, h) for h in recovery.divisors()]
        assert generated == expected_divisors
        least = recovery.least_factor()
        assert least is not None
        assert decode_for_test(fixture, least) == expected_divisors[1]
        factorization = [(decode_for_test(fixture, h), exponent)
                         for h, exponent in recovery.factorization()]
        assert factorization == elementary_factors(number)
        # Nontrivial prefix conditions, plus the all-free endpoint.
        requests = {(0, number), (fixture.width, 0),
                    (fixture.width // 2, number >> (fixture.width // 2)),
                    (1, 1)}
        for free_bits, prefix in sorted(requests):
            if prefix >= (1 << (fixture.width - free_bits)):
                continue
            actual = recovery.prefix_count(free_bits, lambda i: (prefix >> i) & 1)
            assert actual == sum((d >> free_bits) == prefix for d in expected_divisors)
            prefix_cases += 1
        reports.append({"input": number, "control_degree": fixture.field.degree,
                        "raw_descriptors": fixture.total_handles,
                        "divisors": len(expected_divisors),
                        "prime_factorization": factorization,
                        "counts": vars(fixture.counts)})

    assert duplicate_handles > 0 and rejected_odd_buckets > 0
    return {"status": "PASS", "paired_tag_projection_checks": projection_cases,
            "uniform_divisor_bound": bound_checks,
            "streamed_product_checks": product_cases,
            "maximum_product_carry": product_counts.max_carry,
            "maximum_product_column_sum": product_counts.max_column_sum,
            "recovery_instances": reports, "exact_prefix_count_checks": prefix_cases,
            "duplicate_valid_descriptors_removed": duplicate_handles,
            "odd_collision_buckets_rejected": rejected_odd_buckets,
            "coefficient_ring_cardinality": 16,
            "scope": "Finite bucket oracles verify descriptor controls; no huge clique evaluator is implemented.",
            "space_model": "Read-only input, write-only output; handles and counters only in recovery control."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"],
                      "paired_tag_projection_checks": result["paired_tag_projection_checks"],
                      "streamed_product_checks": result["streamed_product_checks"],
                      "recovery_instances": len(result["recovery_instances"]),
                      "exact_prefix_count_checks": result["exact_prefix_count_checks"],
                      "duplicate_valid_descriptors_removed": result["duplicate_valid_descriptors_removed"],
                      "odd_collision_buckets_rejected": result["odd_collision_buckets_rejected"]}))


if __name__ == "__main__":
    main()
