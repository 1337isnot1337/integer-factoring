#!/usr/bin/env python3
"""Exact finite checks for the rational crossing-to-factoring transfer.

These checks exercise the sufficient pad and cutoff inequalities. They do
not supply a crossing family or run an asymptotic factoring algorithm.
"""

from fractions import Fraction
import json


def schedule_log(x0: int, level: int) -> int:
    if x0 < 18 or level < 0:
        raise ValueError("x0 >= 18 and nonnegative level required")
    return (3 * x0 + 6) * 2**level - 6 * level - 2 * x0 - 6


def pad_exponent(p: int, level: int, width: int, bits: int,
                 numerator: int, denominator: int, c0: int) -> int:
    if not (0 < denominator and denominator <= 2 * numerator < 2 * denominator):
        raise ValueError("rho must lie in [1/2,1)")
    if p < 2 or level < 0 or width < 1 or bits < 1 or c0 < 0:
        raise ValueError("invalid clique or pad parameter")
    return (numerator**level * p * p * width // denominator**level + 1
            + c0 * p * p * (p * (width + bits)).bit_length())


def verify() -> dict:
    boundary_cases = pad_cases = 0
    for x0 in (18, 21):
        for level in range(5):
            exponent = schedule_log(x0, level)
            next_exponent = schedule_log(x0, level + 1)
            assert next_exponent == 2 * exponent + 6 * level + 2 * x0
            p = 1 << exponent
            lo, hi = 1 << (40 * exponent), 1 << (40 * next_exponent)
            for n in (lo, lo + 1, hi - 1):
                assert p**40 <= n < p**400
                assert 2**level <= p
                assert n.bit_length() <= 80 * (3 * x0 + 6) * 2**level
                # Extremal widths permitted by the constructor's inequalities.
                width, bits = 4717 * n // (p * p), 6769 * n // p
                size = (p * (width + bits)).bit_length()
                assert p * (width + bits) < (1 << 14) * n
                assert size <= n.bit_length() + 14 <= 414 * p
                assert p * p * size <= 414 * n // p
                boundary_cases += 1
                for a, b in ((1, 2), (3, 5), (2, 3), (7, 10), (3, 4)):
                    scale = Fraction(a, b)**level
                    assert scale >= Fraction(1, p)
                    for c0 in (0, 1, 17, 83):
                        q = pad_exponent(p, level, width, bits, a, b, c0)
                        # Independent rational comparison detects floor direction.
                        real_sufficient = scale * p * p * width + c0 * p * p * size
                        assert q > real_sufficient
                        assert q <= (4718 + 414 * c0) * n * scale
                        pad_cases += 1

    targets = []
    for a, b in ((7, 10), (2, 3), (3, 5), (1, 2)):
        margin = b*b - 2*a*a
        assert margin > 0
        targets.append({"rho_numerator": a, "rho_denominator": b,
                        "strict_half_exponent_margin": margin,
                        "family_supplied": False})
    assert 2 * 3**2 > 4**2  # rho=3/4 does not beat exponent one half.
    assert Fraction(1, 2) > 1 // 2  # A floor alone is not a sufficient pad.
    rejected = 0
    for a, b in ((0, 1), (1, 3), (1, 1), (2, 1), (1, 0), (-1, 2)):
        try:
            pad_exponent(2, 1, 1, 1, a, b, 1)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("invalid rho accepted")

    return {
        "scope": "finite pad and cutoff arithmetic; conditional transfer theorem is separate",
        "schedule_boundary_cases": boundary_cases,
        "pad_cases": pad_cases,
        "invalid_parameter_negative_controls": rejected,
        "other_negative_controls": ["rho=3/4 does not beat one half", "floor without plus one"],
        "conditional_targets": targets,
        "arithmetic": "integer and Fraction only",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
