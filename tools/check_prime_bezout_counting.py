#!/usr/bin/env python3
"""Exact small controls and a separately supplied count oracle.

Dense controls evaluate the full rational Lagrange sum. Wrapper controls
explicitly know the input's finite divisor list; they are not a provider.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb, factorial, lcm, prod
from pathlib import Path

from prime_bezout_counting import (
    is_prime, prime_rectangle_contract, smallest_prime_from_prime_coefficients,
    suffix_rectangle_contract,
)

ROOT = Path(__file__).resolve().parents[1]


class Quadratic:
    def __init__(self, a=0, b=0, d=2):
        self.a, self.b, self.d = Fraction(a), Fraction(b), d

    def of(self, other):
        if isinstance(other, Quadratic):
            assert other.d == self.d
            return other
        return Quadratic(other, d=self.d)

    def __add__(self, other):
        other = self.of(other)
        return Quadratic(self.a + other.a, self.b + other.b, self.d)

    __radd__ = __add__

    def __neg__(self):
        return Quadratic(-self.a, -self.b, self.d)

    def __sub__(self, other):
        return self + -self.of(other)

    def __rsub__(self, other):
        return self.of(other) - self

    def __mul__(self, other):
        other = self.of(other)
        return Quadratic(self.a * other.a + self.d * self.b * other.b,
                         self.a * other.b + self.b * other.a, self.d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.of(other)
        norm = other.a * other.a - self.d * other.b * other.b
        if not norm:
            raise ZeroDivisionError
        return self * Quadratic(other.a / norm, -other.b / norm, self.d)

    def __rtruediv__(self, other):
        return self.of(other) / self

    def residue(self, prime):
        def reduce(q):
            assert q.denominator % prime
            return q.numerator * pow(q.denominator, -1, prime) % prime
        return reduce(self.a), reduce(self.b)


def multiply(first, second):
    result = [Quadratic(d=first[0].d) for _ in range(len(first) + len(second) - 1)]
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    return result


def linear_product(roots, d):
    answer = [Quadratic(1, d=d)]
    for root in roots:
        answer = multiply(answer, [Quadratic(-root, d=d), Quadratic(1, d=d)])
    return answer


def rational_order(value, prime):
    if not value:
        return None
    numerator, denominator = abs(value.numerator), value.denominator
    order = 0
    while numerator % prime == 0:
        numerator //= prime
        order += 1
    while denominator % prime == 0:
        denominator //= prime
        order -= 1
    return order


def dense_leading_coefficient(contract, *, check_pair=False):
    prime, count = contract["prime"], contract["count"]
    if count > 25:
        raise ValueError("Dense controls are limited to 25 nodes per axis")
    aa = [contract["first_a"] + prime * i for i in range(count)]
    bb = [contract["first_b"] + prime * j for j in range(count)]
    omega = Quadratic(0, 1, contract["quadratic_nonresidue"])
    z = contract["n"] + prime ** contract["padding_exponent"] * omega
    scale = prime ** contract["clearing_exponent"] * omega
    answer = Quadratic(d=omega.d)
    if check_pair:
        vv = [Quadratic(d=omega.d) for _ in aa]
    for a in aa:
        denominator = prod((z - a * b for b in bb), start=Quadratic(1, d=omega.d))
        denominator *= prod(a - other for other in aa if other != a)
        weight = scale / denominator
        assert weight.residue(prime) == (int(any(a * b == contract["n"] for b in bb)), 0)
        answer += weight
        if check_pair:
            basis = linear_product([other for other in aa if other != a], omega.d)
            vv = [old + weight * coefficient for old, coefficient in zip(vv, basis)]
    real, imaginary = answer.residue(prime)
    assert imaginary == 0
    if check_pair:
        pp = linear_product(aa, omega.d)
        qq = [Quadratic(1, d=omega.d)]
        for b in bb:
            qq = multiply(qq, [z, Quadratic(-b, d=omega.d)])
        numerator = [-x for x in multiply(vv, qq)]
        numerator[0] += scale
        uu = [Quadratic(d=omega.d) for _ in aa]
        for i in range(len(uu) - 1, -1, -1):
            uu[i] = numerator[i + count]
            for j, coefficient in enumerate(pp):
                numerator[i + j] -= uu[i] * coefficient
        assert all(x.a == x.b == 0 for x in numerator)
        common = linear_product([contract["first_a"]] * (count - 1), omega.d)
        expected_v = [(real * x).residue(prime) for x in common]
        expected_u = [(contract["first_b"] * real * x).residue(prime) for x in common]
        assert [x.residue(prime) for x in vv] == expected_v
        assert [x.residue(prime) for x in uu] == expected_u
        orders = [order for x in vv + uu for coordinate in (x.a, x.b)
                  if (order := rational_order(coordinate, prime)) is not None]
        assert min(orders) >= 0
        primitive_scale = max(0, contract["clearing_exponent"] - min(orders))
        assert (primitive_scale == contract["clearing_exponent"]) == bool(real)
    return real, aa, bb


def dense_controls():
    cases = derivative_units = padding_controls = complete_pairs = 0
    examples = []
    for prime in (3, 5, 7):
        for exponent in (0, 1, 2):
            if prime ** exponent > 25:
                continue
            for n in (11, 35, 77, 143):
                if n % prime == 0:
                    continue
                for a0 in (1, 2, prime + 1):
                    b0 = n * pow(a0, -1, prime) % prime
                    contract = prime_rectangle_contract(n, prime, a0, b0, exponent)
                    check_pair = contract["count"] <= 7 and cases % 3 == 0
                    value, aa, bb = dense_leading_coefficient(contract, check_pair=check_pair)
                    complete_pairs += check_pair
                    expected = sum(a * b == n for a in aa for b in bb) % prime
                    assert value == expected
                    count, order = contract["count"], contract["derivative_order"]
                    for nodes in (aa, bb):
                        for a in nodes:
                            derivative = prod(a - other for other in nodes if other != a)
                            assert derivative % (prime ** order) == 0
                            unit = derivative // (prime ** order)
                            assert unit % prime == pow(-1, exponent, prime)
                            derivative_units += 1
                    if count <= 7 and cases % 7 == 0:
                        other = dict(contract)
                        other["padding_exponent"] += 1
                        other["clearing_exponent"] += 1
                        assert dense_leading_coefficient(other)[0] == value
                        padding_controls += 1
                    if expected and len(examples) < 8:
                        examples.append({"n": n, "prime": prime, "count": count,
                                         "first_a": a0, "first_b": b0,
                                         "leading_coefficient": value,
                                         "factor_pairs": [[a, b] for a in aa for b in bb
                                                          if a * b == n]})
                    cases += 1
    return {"dense_exact_quadratic_coefficients": cases,
            "complete_pairs_and_maximal_scales": complete_pairs,
            "literal_derivative_unit_controls": derivative_units,
            "independent_padding_controls": padding_controls, "examples": examples}


def wrapper_controls():
    total_queries = used_oracles = cases = 0
    examples = []
    inputs = list(range(2, 258)) + [17 * 19, 23 * 29, 101 * 103, 107 * 109, 104729]
    for n in inputs:
        divisors = [a for a in range(2, n // 2 + 1) if n % a == 0]

        def finite_known_factor_oracle(contract):
            nonlocal total_queries
            total_queries += 1
            prime = contract["prime"]
            first_a, first_b, count = contract["first_a"], contract["first_b"], contract["count"]
            last_a, last_b = first_a + prime * (count - 1), first_b + prime * (count - 1)
            return sum(first_a <= a <= last_a and (a - first_a) % prime == 0
                       and first_b <= n // a <= last_b
                       and (n // a - first_b) % prime == 0 for a in divisors) % prime

        result = smallest_prime_from_prime_coefficients(n, finite_known_factor_oracle)
        expected = min(divisors) if divisors else n
        assert result["smallest_prime_factor"] == expected
        used_oracles += bool(result["coefficient_queries"])
        if n in (169, 221, 257, 323, 667, 10403, 11663, 104729):
            examples.append({"n": n, "proper_divisors": divisors, **result})
        cases += 1
    return {"conditional_wrapper_inputs": cases, "inputs_requiring_supplied_oracle": used_oracles,
            "finite_known_factor_oracle_queries": total_queries, "selected_controls": examples}


def primorial_controls():
    checks = 0
    for middle in range(1, 81):
        bound = 2 * middle
        central = comb(bound, middle)
        common_multiple = lcm(*range(1, bound + 1))
        assert common_multiple % central == 0
        assert central >= 2 ** middle
        primorial = prod(p for p in range(2, bound + 1) if is_prime(p))
        assert common_multiple <= primorial ** bound.bit_length()
        checks += 1
    for prime in (3, 5, 7, 11):
        for exponent in range(1, 4):
            count = prime ** exponent
            alpha = factorial(count - 1)
            while alpha % prime == 0:
                alpha //= prime
            assert alpha % prime == pow(-1, exponent, prime)
            for index in range(count):
                assert comb(count - 1, index) % prime == pow(-1, index, prime)
                checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    dense, wrapper = dense_controls(), wrapper_controls()
    large = suffix_rectangle_contract((1 << 401) + 41, 17, 11, 3)
    receipt = {"theorem": "Prime-adic precision one returns factor counts",
               "scope": "Exact dense controls and explicitly supplied finite factor-count oracles; no global compressed coefficient provider or improved factoring score.",
               **dense, **wrapper, "primorial_lucas_and_factorial_controls": primorial_controls(),
               "large_encoded_contract": {"input_bits": large["n"].bit_length(),
                   "prime": large["prime"], "precision": large["precision"],
                   "count_bits": large["count"].bit_length(),
                   "padding_exponent": large["padding_exponent"],
                   "clearing_exponent_bits": large["clearing_exponent"].bit_length(),
                   "scope": "Contract construction only; no large coefficient computed"}}
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
