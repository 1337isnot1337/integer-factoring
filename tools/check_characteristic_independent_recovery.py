"""Exact finite checks for characteristic-independent divisor recovery.

The weighted evaluator is supplied by literal divisor sums in this checker.
This validates recovery identities and controls, not an asymptotic factorer.
Only the Python standard library is used.
"""

from __future__ import annotations

import argparse
import json
from functools import lru_cache
from pathlib import Path
from random import Random


def digits(value: int, base: int, length: int) -> tuple[int, ...]:
    answer = []
    for _ in range(length):
        value, digit = divmod(value, base)
        answer.append(digit)
    assert value == 0
    return tuple(answer)


def encode(values, base: int) -> int:
    answer = 0
    for value in reversed(tuple(values)):
        answer = answer * base + value
    return answer


def polynomial_remainder(dividend, divisor, prime: int):
    answer = list(dividend)
    while answer and answer[-1] == 0:
        answer.pop()
    inverse = pow(divisor[-1], -1, prime)
    while len(answer) >= len(divisor):
        shift = len(answer) - len(divisor)
        coefficient = answer[-1] * inverse % prime
        for j, value in enumerate(divisor):
            answer[shift + j] = (answer[shift + j] - coefficient * value) % prime
        while answer and answer[-1] == 0:
            answer.pop()
    return tuple(answer)


@lru_cache(maxsize=None)
def first_irreducible(prime: int, degree: int) -> tuple[int, ...]:
    """The exponential-time, linear-storage construction used in the proof."""
    assert prime in (2, 3, 5, 7) and degree >= 1
    for candidate in range(prime**degree):
        polynomial = digits(candidate, prime, degree) + (1,)
        good = True
        for divisor_degree in range(1, degree // 2 + 1):
            for lower in range(prime**divisor_degree):
                divisor = digits(lower, prime, divisor_degree) + (1,)
                if not polynomial_remainder(polynomial, divisor, prime):
                    good = False
                    break
            if not good:
                break
        if good:
            return polynomial
    raise AssertionError("No irreducible polynomial found")


class SmallField:
    def __init__(self, prime: int, degree: int):
        self.prime = prime
        self.degree = degree
        self.order = prime**degree
        self.modulus = first_irreducible(prime, degree)

    def coordinates(self, value: int):
        return digits(value, self.prime, self.degree)

    def scalar(self, value: int):
        return value % self.prime

    def add(self, left: int, right: int):
        return encode(((a + b) % self.prime for a, b in zip(
            self.coordinates(left), self.coordinates(right))), self.prime)

    def neg(self, value: int):
        return encode(((-a) % self.prime for a in self.coordinates(value)), self.prime)

    def mul(self, left: int, right: int):
        a, b = self.coordinates(left), self.coordinates(right)
        product = [0] * (2 * self.degree - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                product[i + j] = (product[i + j] + x * y) % self.prime
        residue = polynomial_remainder(product, self.modulus, self.prime)
        return encode(residue, self.prime)

    def power(self, value: int, exponent: int):
        answer = 1
        while exponent:
            if exponent & 1:
                answer = self.mul(answer, value)
            value = self.mul(value, value)
            exponent >>= 1
        return answer

    def hash_bits(self, integer: int, width: int, point: int):
        answer = 0
        for bit in range(width - 1, -1, -1):
            answer = self.add(self.mul(answer, point), (integer >> bit) & 1)
        return answer


class CharacterAlgebra:
    """Odd fields, or K[omega]/(omega^2+omega+1) when char K=2."""

    def __init__(self, prime: int, degree: int = 1):
        self.base = SmallField(prime, degree)
        self.prime = prime
        self.radix = 3 if prime == 2 else 2
        self.extended = prime == 2
        self.order = self.base.order ** (2 if self.extended else 1)
        self.omega = self.base.order if self.extended else self.base.neg(1)

    def scalar(self, value: int):
        return self.base.scalar(value)

    def add(self, left: int, right: int):
        if not self.extended:
            return self.base.add(left, right)
        q = self.base.order
        return self.base.add(left % q, right % q) + q * self.base.add(left // q, right // q)

    def mul(self, left: int, right: int):
        if not self.extended:
            return self.base.mul(left, right)
        q = self.base.order
        a, b, c, d = left % q, left // q, right % q, right // q
        bd = self.base.mul(b, d)
        constant = self.base.add(self.base.mul(a, c), bd)
        linear = self.base.add(self.base.add(self.base.mul(a, d), self.base.mul(b, c)), bd)
        return constant + q * linear

    def power(self, value: int, exponent: int):
        answer = 1
        while exponent:
            if exponent & 1:
                answer = self.mul(answer, value)
            value = self.mul(value, value)
            exponent >>= 1
        return answer

    def inverse_integer(self, integer: int):
        value = self.base.scalar(integer)
        assert value != 0
        return self.base.power(value, self.base.order - 2)

    def character(self, control: SmallField, character: int, value: int):
        assert control.prime == self.radix
        dot = sum(x * y for x, y in zip(control.coordinates(character), control.coordinates(value)))
        return self.power(self.omega, dot % self.radix)


def divisor_bit_bound(width: int):
    if width < 4:
        return width
    half_log = (width.bit_length() - 1) // 2
    return min(width, (width + half_log - 1) // half_log
               + (1 << half_log) * (width + 1).bit_length())


def divisors(integer: int):
    return [d for d in range(1, integer + 1) if integer % d == 0]


def make_control(radix: int, bound: int):
    degree = 1
    while radix**degree <= bound:
        degree += 1
    return SmallField(radix, degree)


def tag_value(divisor: int, integer: int, tag):
    if tag is None:
        return 1
    side, bit = tag
    value = divisor if side == 0 else integer // divisor
    return (value >> bit) & 1


def weighted_sum(integer, factors, control, algebra, point, character, tag):
    answer = 0
    width = integer.bit_length()
    for divisor in factors:
        if tag_value(divisor, integer, tag):
            h = control.hash_bits(divisor, width, point)
            answer = algebra.add(answer, algebra.character(control, character, h))
    return answer


def fourier_bucket(integer, factors, control, algebra, point, bucket, tag):
    answer = 0
    for character in range(control.order):
        coefficient = algebra.character(control, character, control.neg(bucket))
        value = weighted_sum(integer, factors, control, algebra, point, character, tag)
        answer = algebra.add(answer, algebra.mul(coefficient, value))
    return algebra.mul(answer, algebra.inverse_integer(control.order))


def verify_virtual_product(integer: int, bit_query):
    """Column multiplication recomputes every bit and retains only one carry."""
    width = integer.bit_length()
    positive = [False, False]
    for side in range(2):
        for bit in range(width):
            value = bit_query(side, bit)
            if value not in (0, 1):
                return False
            positive[side] |= value == 1
    if not all(positive):
        return False
    carry = 0
    for column in range(2 * width):
        total = carry
        for i in range(max(0, column - width + 1), min(width, column + 1)):
            total += bit_query(0, i) * bit_query(1, column - i)
        if total % 2 != ((integer >> column) & 1):
            return False
        carry = total // 2
    return carry == 0


def field_and_character_controls():
    inverses = orthogonality = 0
    zero_divisor_control = False
    for prime, degree in ((3, 1), (5, 1), (7, 1), (2, 1), (2, 2)):
        algebra = CharacterAlgebra(prime, degree)
        assert algebra.power(algebra.omega, algebra.radix) == 1
        assert sum(1 for x in range(algebra.order) if algebra.mul(x, 0) != 0) == 0
        if prime == 2 and degree == 2:
            zero_divisor_control = any(x and y and algebra.mul(x, y) == 0
                                      for x in range(algebra.order) for y in range(algebra.order))
            assert zero_divisor_control
        for control_degree in range(1, 4):
            control = SmallField(algebra.radix, control_degree)
            for value in range(1, control.order):
                assert control.mul(value, control.power(value, control.order - 2)) == 1
                inverses += 1
            normalization = algebra.inverse_integer(control.order)
            for value in range(control.order):
                total = 0
                for character in range(control.order):
                    total = algebra.add(total, algebra.character(control, character, value))
                assert algebra.mul(total, normalization) == (value == 0)
                orthogonality += 1
    # An absent normalization and a trivial characteristic-two character fail.
    odd = CharacterAlgebra(5)
    assert odd.scalar(8) != 1
    even = CharacterAlgebra(2)
    assert even.scalar(9) == 1  # trivial character would incorrectly detect every bucket
    return {"nonzero_control_inverses": inverses,
            "character_orthogonality": orthogonality,
            "split_quadratic_zero_divisor_control": zero_divisor_control,
            "missing_normalization_control": True,
            "trivial_character_control": True}


def recovery_controls():
    rng = Random(20261030)
    totals = {"integers": 0, "descriptors": 0, "accepted_descriptors": 0,
              "fourier_buckets": 0, "false_singletons_rejected": 0,
              "cancelled_nonempty_buckets": 0, "isolation_checks": 0}
    examples = {}
    for prime, base_degree in ((2, 1), (2, 2), (3, 1), (5, 1)):
        algebra = CharacterAlgebra(prime, base_degree)
        for integer in (2, 3, 4, 6, 8, 12, 18, 24, 30, 35):
            width = integer.bit_length()
            factors = divisors(integer)
            # The finite experiment supplies an exact support bound to keep it small.
            # Uniform bound (2) is independently checked below.
            control = make_control(algebra.radix, width * len(factors))
            accepted = set()
            isolated = set()
            samples = []
            for point in range(control.order):
                buckets = {}
                for divisor in factors:
                    h = control.hash_bits(divisor, width, point)
                    buckets.setdefault(h, []).append(divisor)
                for bucket in range(control.order):
                    members = buckets.get(bucket, [])
                    totals["descriptors"] += 1
                    if len(members) == 1:
                        isolated.add(members[0])
                    count = algebra.scalar(len(members))
                    if members and count == 0:
                        totals["cancelled_nonempty_buckets"] += 1
                        examples.setdefault("cancellation", {"N": integer, "characteristic": prime,
                                            "point": point, "bucket": bucket, "members": members})
                    if count != 1:
                        continue
                    def query(side, bit):
                        return algebra.scalar(sum(tag_value(d, integer, (side, bit)) for d in members))
                    good = verify_virtual_product(integer, query)
                    if good:
                        value = sum(query(0, bit) << bit for bit in range(width))
                        assert integer % value == 0
                        accepted.add(value)
                        totals["accepted_descriptors"] += 1
                    elif len(members) > 1:
                        totals["false_singletons_rejected"] += 1
                        examples.setdefault("false_singleton", {"N": integer, "characteristic": prime,
                                            "point": point, "bucket": bucket, "members": members})
                    if len(samples) < 3:
                        samples.append((point, bucket))
            assert isolated == set(factors)
            assert accepted == set(factors)
            totals["isolation_checks"] += len(factors)
            ordered = sorted(accepted)
            assert ordered == factors and min(x for x in ordered if x > 1) == factors[1]
            for low, high in ((0, integer), (2, integer // 2), (integer // 3, integer)):
                assert sum(low <= d <= high for d in ordered) == sum(low <= d <= high for d in factors)
            samples += [(rng.randrange(control.order), rng.randrange(control.order)) for _ in range(2)]
            for point, bucket in samples:
                members = [d for d in factors if control.hash_bits(d, width, point) == bucket]
                for tag in (None, (0, rng.randrange(width)), (1, rng.randrange(width))):
                    expected = algebra.scalar(sum(tag_value(d, integer, tag) for d in members))
                    actual = fourier_bucket(integer, factors, control, algebra, point, bucket, tag)
                    assert actual == expected
                    totals["fourier_buckets"] += 1
            totals["integers"] += 1
    assert totals["false_singletons_rejected"] > 0
    assert totals["cancelled_nonempty_buckets"] > 0
    return {**totals, "negative_control_examples": examples}


def elementary_bound_controls():
    # A sieve supplies independent divisor counts for every integer below this bound.
    limit = 10000
    counts = [0] * (limit + 1)
    for divisor in range(1, limit + 1):
        for multiple in range(divisor, limit + 1, divisor):
            counts[multiple] += 1
    for integer in range(2, limit + 1):
        assert counts[integer] <= 1 << divisor_bit_bound(integer.bit_length())
    large = []
    for exponent in (8, 16, 32, 64, 128, 256):
        width = 1 << exponent
        q = divisor_bit_bound(width)
        assert 0 < q <= width
        if exponent >= 32:
            assert q * exponent <= 3 * width
        large.append({"log2_bit_length": exponent,
                      "bound_bits": str(q), "three_n_over_log_n": q * exponent <= 3 * width})
    return {"sieved_integers": limit - 1, "large_integer_controls": large}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"status": "PASS", "scope": "Finite exact recovery identities; literal divisor oracle only.",
              "field_controls": field_and_character_controls(),
              "divisor_bound_controls": elementary_bound_controls(),
              "recovery_controls": recovery_controls()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
