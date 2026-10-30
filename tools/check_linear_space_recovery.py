#!/usr/bin/env python3
"""Exact controls for coordinatewise constant-ring witness recovery.

The finite evaluator reads an explicit witness fixture. It is an oracle
for checking substitutions, not an implementation of the huge clique solver.
"""

import argparse
from itertools import product
import json
from pathlib import Path

from check_bit_tagged_witness_recovery import (
    PrimeExtension, additive_bucket,
)


class TernaryField:
    """No field table; every polynomial coefficient is reduced immediately."""

    def __init__(self, degree):
        self.degree, self.size = degree, 3**degree
        self.polynomial = next(
            self.coordinates(candidate)+(1,)
            for candidate in range(self.size)
            if all(self.remainder(self.coordinates(candidate)+(1,),
                                  self.coordinates(divisor, order)+(1,))
                   for order in range(1, degree//2+1)
                   for divisor in range(3**order)))

    def coordinates(self, value, length=None):
        length = self.degree if length is None else length
        result = []
        for _ in range(length):
            value, digit = divmod(value, 3)
            result.append(digit)
        return tuple(result)

    @staticmethod
    def encode(coefficients):
        result = 0
        for value in reversed(coefficients):
            result = 3*result+value % 3
        return result

    @staticmethod
    def remainder(value, divisor):
        result = list(value)
        for i in range(len(result)-1, len(divisor)-2, -1):
            scalar = result[i]
            offset = i-len(divisor)+1
            for j, coefficient in enumerate(divisor):
                result[offset+j] = (result[offset+j]-scalar*coefficient) % 3
        return any(result[:len(divisor)-1])

    def add(self, first, second):
        return self.encode(tuple((a+b) % 3 for a, b in
                                 zip(self.coordinates(first), self.coordinates(second))))

    def multiply(self, first, second):
        a, b = self.coordinates(first), self.coordinates(second)
        result = [0]*(2*self.degree-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                result[i+j] = (result[i+j]+x*y) % 3
        for i in range(len(result)-1, self.degree-1, -1):
            scalar = result[i]
            for j, coefficient in enumerate(self.polynomial):
                position = i-self.degree+j
                result[position] = (result[position]-scalar*coefficient) % 3
        return self.encode(result[:self.degree])

    def power(self, value, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.multiply(result, value)
            value = self.multiply(value, value)
            exponent >>= 1
        return result

    def dot(self, first, second):
        result = 0
        for _ in range(self.degree):
            first, a = divmod(first, 3)
            second, b = divmod(second, 3)
            result = (result+a*b) % 3
        return result

    def hash(self, witness, point, width):
        result = 0
        for j in range(width-1, -1, -1):
            result = self.add(self.multiply(result, point), (witness >> j) & 1)
        return result


# F4=F2[w]/(w²+w+1); fixed-size arithmetic is the entire coefficient ring.
def f4_mul(a, b):
    raw = 0
    for j in range(2):
        if (b >> j) & 1:
            raw ^= a << j
    return raw ^ (7 if raw & 4 else 0)


PHASES = (1, 2, 3)


def ring_mul(first, second):
    a, b = first
    c, d = second
    return f4_mul(a, c), f4_mul(a, d) ^ f4_mul(b, c)


def query_fixture(witnesses, width, variable):
    """Read each variable on demand, keeping one ring term and one sum."""
    scalar = tag = 0
    for witness in witnesses:
        value = (1, 0)
        for j in range(width):
            if (witness >> j) & 1:
                value = ring_mul(value, variable(j))
        scalar ^= value[0]
        tag ^= value[1]
    return scalar, tag


def bucket_bit(witnesses, width, point, bucket, bit, control):
    scalar = tag = 0
    for character in range(control.size):
        def variable(j):
            # Recompute this power; there is no width-by-field-degree array.
            phase = PHASES[control.dot(character, control.power(point, j))]
            return phase, phase if j == bit else 0
        a, b = query_fixture(witnesses, width, variable)
        factor = PHASES[-control.dot(character, bucket) % 3]
        scalar ^= f4_mul(factor, a)
        tag ^= f4_mul(factor, b)
    return scalar, tag


def recover_minimum(witnesses, width, verifier):
    bound = 1 << max(0, (len(witnesses)-1).bit_length())
    degree = 1
    while 3**degree <= width*(bound-1):
        degree += 1
    control = TernaryField(degree)
    best, bit_queries, rejected = None, 0, 0
    isolated = set()  # Verification instrumentation only, not recovery state.
    for point in range(control.size):
        for bucket in range(control.size):
            candidate, valid = 0, True
            for bit in range(width):
                scalar, value = bucket_bit(witnesses, width, point, bucket, bit, control)
                bit_queries += 1
                valid &= scalar == 1 and value in (0, 1)
                if value == 1:
                    candidate |= 1 << bit
            if valid:
                if verifier(candidate):
                    best = candidate if best is None else min(best, candidate)
                else:
                    rejected += 1
            # Independent finite check that every fixture witness is isolated.
            members = [w for w in witnesses if control.hash(w, point, width) == bucket]
            if len(members) == 1:
                isolated.add(members[0])
    assert isolated == set(witnesses)
    return best, bit_queries, rejected


def elementary_factorization(number):
    answer, divisor = [], 2
    while divisor*divisor <= number:
        exponent = 0
        while number % divisor == 0:
            number //= divisor
            exponent += 1
        if exponent:
            answer.append((divisor, exponent))
        divisor += 1
    if number > 1:
        answer.append((number, 1))
    return answer


def divisors_from_factors(factors):
    """Mixed-radix enumeration: one exponent vector and one current product."""
    exponents = [0]*len(factors)
    while True:
        value = 1
        for (prime, _), exponent in zip(factors, exponents):
            value *= prime**exponent
        yield value
        for i, (_, limit) in enumerate(factors):
            if exponents[i] < limit:
                exponents[i] += 1
                break
            exponents[i] = 0
        else:
            return


def verify():
    fields = comparisons = projections = phase_queries = 0
    coefficient = PrimeExtension(2, 2)
    add = [[a ^ b for b in range(4)] for a in range(4)]
    mul = [[f4_mul(a, b) for b in range(4)] for a in range(4)]
    assert coefficient.polynomial == (1, 1, 1)
    for degree in range(1, 5):
        control, other = TernaryField(degree), PrimeExtension(3, degree)
        assert control.polynomial == other.polynomial
        fields += 1
        for a, b in product(range(control.size), repeat=2):
            assert control.add(a, b) == other.add(a, b)
            assert control.multiply(a, b) == other.multiply(a, b)
            comparisons += 2
        for a in range(control.size):
            assert control.power(a, control.size) == a
    examples = ((3, (1, 3, 6)), (4, (0, 5, 9, 14)), (5, (0, 17)))
    for width, witnesses in examples:
        degree = 1
        while 3**degree <= width*(len(witnesses)-1):
            degree += 1
        control, other = TernaryField(degree), PrimeExtension(3, degree)
        for point in range(control.size):
            basis = other.bit_basis(point, width)  # Independent comparison only.
            for j in range(width):
                assert control.power(point, j) == basis[j]
                phase_queries += 1
            for bucket in range(control.size):
                expected = additive_bucket(witnesses, basis, bucket, other,
                                           coefficient, 2, add, mul)
                direct = [w for w in witnesses if control.hash(w, point, width) == bucket]
                assert expected == (len(direct) % 2,)+tuple(
                    sum((w >> bit) & 1 for w in direct) % 2 for bit in range(width))
                for bit in range(width):
                    result = bucket_bit(witnesses, width, point, bucket, bit, control)
                    assert result == (expected[0], expected[bit+1])
                    projections += 1
    recovery_queries = rejected = 0
    finite_numbers = (2, 3, 4, 6, 9, 15, 17, 25, 49, 77)
    for number in finite_numbers:
        witnesses = tuple(d for d in range(2, number) if number % d == 0)
        result, calls, bad = recover_minimum(
            witnesses, number.bit_length(), lambda d: 1 < d < number and number % d == 0)
        result = number if result is None else result
        expected = next(d for d in range(2, number+1) if number % d == 0)
        assert result == expected
        recovery_queries += calls
        rejected += bad
    # An odd collision can decode a non-witness; candidate checks are necessary.
    result, calls, bad = recover_minimum((1, 2, 4), 3, lambda d: d in (1, 2, 4))
    assert result == 1 and bad > 0
    recovery_queries += calls
    rejected += bad
    divisor_checks = prefix_checks = 0
    for number in range(2, 1025):
        factors = elementary_factorization(number)
        generated = list(divisors_from_factors(factors))  # Independent test fixture.
        direct = [d for d in range(1, number+1) if number % d == 0]
        assert sorted(generated) == direct and len(generated) == len(set(generated))
        divisor_checks += 1
        width = number.bit_length()
        for free_bits in range(width+1):
            prefixes = {0, number >> free_bits, (number >> free_bits)//2}
            for prefix in prefixes:
                assert sum((d >> free_bits) == prefix for d in generated) == sum(
                    (d >> free_bits) == prefix for d in direct)
                prefix_checks += 1
        # A simple exact description budget for prime and exponent fields.
        assert sum(prime.bit_length()+exponent.bit_length()
                   for prime, exponent in factors) <= 3*width
    return {"status": "PASS", "control_fields": fields,
            "independent_field_operations": comparisons,
            "on_demand_power_checks": phase_queries,
            "coordinatewise_bucket_comparisons": projections,
            "least_factor_cases": len(finite_numbers),
            "recovery_bit_queries": recovery_queries,
            "invalid_collision_candidates_rejected": rejected,
            "factor_to_divisor_enumerations": divisor_checks,
            "exact_divisor_prefix_counts": prefix_checks,
            "coefficient_ring_cardinality": 16,
            "scope": "Finite recovery and enumeration checks; global bounds use the proof."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = verify()
    encoded = json.dumps(report, indent=2)+"\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
