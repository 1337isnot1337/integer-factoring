#!/usr/bin/env python3
"""Finite exact checks for the handwritten exponent-one-half composition.

No large clique is materialized. The prefix-search checks use an independent
small divisor-count oracle, not the Appendix B evaluator. The proof in
research/half-rate-factoring.md supplies the all-input and cost arguments.
"""

from fractions import Fraction
import json
from math import isqrt


def log_k(x0: int, level: int) -> int:
    if x0 < 18 or level < 0:
        raise ValueError("source schedule requires x0 >= 18 and level >= 0")
    return (3 * x0 + 6) * 2**level - 6 * level - 2 * x0 - 6


def top_parameters(n: int, p: int) -> dict[str, int]:
    if p < 2 or n < p**40:
        raise ValueError("top constructor requires p >= 2 and p**40 <= n")
    # Least s with n <= 2**(4*s+2), computed without real logarithms.
    s = max(0, ((n - 1).bit_length() + 1) // 4)
    scale = 1 << s
    a = 2 * scale
    f = 8 * scale**2 + 1
    g = 3 * p**8
    degree = 7 * (a // g + 2) * a * f + 3 * (a * a // g + 2) * f
    ordered_total = 10 * a * a * f
    width = max(1, (p * ordered_total + g * degree) // p**3)
    cluster_bits = 8 * (a // g + 2) * a * f + 2 * (a * a // g + 2) * f
    local_bits = cluster_bits * (2 * g // p)
    return {"s": s, "S": scale, "G": g, "lambda": width,
            "H": local_bits, "B": local_bits + 1}


def explicit_pad_exponent(p: int, level: int, width: int, bits: int, c0: int) -> int:
    if p < 2 or level < 0 or min(width, bits) < 1 or c0 < 0:
        raise ValueError("invalid pad parameters")
    r = 1 << (level // 2)
    return p * p * width // r + 1 + c0 * p * p * (p * (width + bits)).bit_length()


def factor_from_prefix_counts(n_integer: int, count) -> int:
    """Prefix driver only; `count(k,x)` must give the exact divisor count."""
    if n_integer < 0:
        raise ValueError("natural number required")
    if n_integer < 2:
        return 2 if n_integer == 0 else 1
    prefix = 0
    for k in reversed(range(n_integer.bit_length())):
        left = 2 * prefix
        nontrivial = count(k, left) - ((1 >> k) == left)
        if nontrivial < 0:
            raise ValueError("prefix oracle violates the divisor-one count")
        prefix = left if nontrivial else left + 1
    return prefix


def verify() -> dict:
    boundary_cases = pad_cases = 0
    max_length_bits = 0
    for x0 in (18, 19, 32):
        for level in range(5):
            exponent = log_k(x0, level)
            next_exponent = log_k(x0, level + 1)
            assert next_exponent == 2 * exponent + 6 * level + 2 * x0
            assert exponent < next_exponent <= 10 * exponent
            p = 1 << exponent
            r = 1 << (level // 2)
            assert r <= p
            lo = 1 << (40 * exponent)
            hi = 1 << (40 * next_exponent)
            for n in (lo, lo + 1, 2 * lo, hi - 1):
                assert p**40 <= n < p**400
                params = top_parameters(n, p)
                scale, width, bits = params["S"], params["lambda"], params["B"]
                assert n <= 4 * scale**4 <= 16 * n
                assert params["s"] == 0 or n > 4 * (scale // 2)**4
                assert p**9 <= scale
                assert p * params["G"] * scale**3 <= 3 * scale**4
                assert p * p * width <= 4717 * n
                assert p * params["H"] <= 6768 * n
                assert p * bits <= 6769 * n
                assert p * (width + bits) < (1 << 14) * n
                assert n.bit_length() <= 160 * (3 * x0 + 6) * r * r
                for c0 in (0, 1, 7, 41):
                    q = explicit_pad_exponent(p, level, width, bits, c0)
                    residual = q - c0 * p * p * (p * (width + bits)).bit_length()
                    assert residual * r > p * p * width
                    assert q * r <= (4718 + 414 * c0) * n
                    pad_cases += 1
                max_length_bits = max(max_length_bits, n.bit_length())
                boundary_cases += 1

    # The floor-plus-one is necessary at every nonintegral quotient.
    assert Fraction(2 * 2 * 1, 8) > (2 * 2 * 1 // 8)

    try:
        top_parameters(2**40 - 1, 2)
    except ValueError:
        pass
    else:
        raise AssertionError("cutoff negative control was accepted")

    prefix_queries = 0
    for value in range(2, 1025):
        divisors = [d for d in range(1, value + 1) if value % d == 0]
        queries = []

        def count(k, x):
            assert 0 <= k < value.bit_length()
            assert 0 <= x < 1 << (value.bit_length() - k)
            queries.append((k, x))
            return sum((d >> k) == x for d in divisors)

        result = factor_from_prefix_counts(value, count)
        expected = next((d for d in range(2, isqrt(value) + 1) if value % d == 0), value)
        assert result == expected
        assert len(queries) == value.bit_length()
        prefix_queries += len(queries)
    assert factor_from_prefix_counts(0, None) == 2
    assert factor_from_prefix_counts(1, None) == 1

    return {
        "scope": "finite arithmetic and prefix driver; no large evaluator or asymptotic test",
        "schedule_and_top_boundary_cases": boundary_cases,
        "explicit_pad_cases": pad_cases,
        "largest_tested_bit_length_parameter_bits": max_length_bits,
        "prefix_inputs": 1023,
        "prefix_queries": prefix_queries,
        "small_input_convention_cases": 2,
        "negative_controls": ["missing floor-plus-one", "input below top cutoff"],
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
