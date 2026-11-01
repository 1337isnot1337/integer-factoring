#!/usr/bin/env python3
"""Check the one-coefficient reduction using an explicitly factored oracle.

The finite oracle enumerates small integer divisors. It is evidence for
the reduction and normalization, not a fast global coefficient provider.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import json
from math import comb, factorial
from pathlib import Path

from balanced_bezout_modular import _unit_product
from interval_unit_jets import QuadraticJets
from leading_bezout_factoring import (
    coefficient_has_factor, leading_coefficient_contract,
    smallest_prime_from_leading_coefficient,
)


@lru_cache(maxsize=None)
def factorial_odd_unit(size_bits, precision):
    """Odd part of (2^size_bits-1)!, using complete odd progressions."""
    assert precision >= 3
    ring = QuadraticJets(precision, 0)
    result = 1
    for depth in range(1, min(size_bits, precision) + 1):
        value = _unit_product(ring, (0, 1 << (depth - 1)),
                              ring.constant(1), ring.constant(2))[0]
        assert value[1] == 0
        result = result * value[0] % ring.modulus
    return result


@lru_cache(maxsize=None)
def literal_binomial_mod(size_bits, index, precision):
    # Only the small known-divisor controls call this finite reference.
    if index > 4096:
        raise ValueError("Literal binomial control exceeds its declared finite cap")
    return comb((1 << size_bits) - 1, index) % (1 << precision)


def known_factor_coefficient(contract, proper_divisors):
    n, first = contract["n"], contract["first_a"]
    precision, size_bits = contract["precision"], contract["domain_size_bits"]
    count, modulus = contract["rows"], 1 << precision
    omega = factorial_odd_unit(size_bits, precision)
    common_inverse = pow(omega * omega % modulus, -1, modulus)
    coefficient = 0
    rows = []
    for a in proper_divisors:
        if a < first:
            continue
        b = n // a
        i, j = (a - first) // 2, (b - 3) // 2
        relative = (pow(a, -(count - 1), modulus) * (-1) ** (i + j)
                    * literal_binomial_mod(size_bits, i, precision)
                    * literal_binomial_mod(size_bits, j, precision)) % modulus
        assert relative == a % modulus
        coefficient = (coefficient + common_inverse * relative) % modulus
        rows.append(a)
    recovered_sum = omega * omega * coefficient % modulus
    assert recovered_sum == sum(rows) < n * n
    assert coefficient_has_factor(contract, coefficient) == bool(rows)
    return coefficient, len(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    factorial_controls = 0
    for size_bits in range(1, 7):
        value = factorial((1 << size_bits) - 1)
        value >>= (value & -value).bit_length() - 1
        for precision in range(3, 13):
            assert factorial_odd_unit(size_bits, precision) == value % (1 << precision)
            factorial_controls += 1
    records = []
    calls = factor_rows = 0
    for n in range(3, 258, 2):
        proper = [a for a in range(3, n, 2) if n % a == 0 and n // a >= 3]

        def oracle(contract):
            nonlocal calls, factor_rows
            value, support = known_factor_coefficient(contract, proper)
            calls += 1
            factor_rows += support
            return value

        result = smallest_prime_from_leading_coefficient(n, oracle)
        expected = next(d for d in range(2, n + 1) if n % d == 0)
        assert result["smallest_prime_factor"] == expected
        assert result["coefficient_queries"] <= n.bit_length() + 1
        if n in (3, 9, 15, 25, 45, 81, 105, 169, 225, 257):
            records.append({"n": n, "proper_divisors": proper, **result})
    for n in (2, 4, 16, 100):
        result = smallest_prime_from_leading_coefficient(n, lambda _: (_ for _ in ()).throw(AssertionError()))
        assert result["smallest_prime_factor"] == 2 and result["coefficient_queries"] == 0
    long_n = (1 << 401) + 45
    large = leading_coefficient_contract(long_n, 3)
    large_summary = {"input_bits": large["input_bits"], "precision": large["precision"],
                     "domain_log2_size": large["domain_size_bits"],
                     "padding_exponent": large["padding_exponent"],
                     "clearing_exponent_bits": large["clearing_exponent"].bit_length(),
                     "requested_coefficients": large["requested_coefficients"],
                     "scope": "Contract construction only; no large coefficient or factor computed"}
    result = {"theorem": "One global leading Bezout coefficient suffices for factoring",
              "literal_factorial_unit_controls": factorial_controls,
              "odd_integer_reduction_controls": 128,
              "finite_known_factor_oracle_queries": calls,
              "exact_relative_factor_weight_controls": factor_rows,
              "selected_controls": records, "large_encoded_contract": large_summary,
              "scope": "The oracle is explicitly supplied with finite divisor lists. No fast global coefficient provider or improved factoring score is claimed."}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
