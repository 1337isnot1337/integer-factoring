#!/usr/bin/env python3
"""Exact ceilings and volume controls for the adaptive prime schedule.

The finite replay does not estimate the effective constants in the
Thorner--Zaman theorem. Their asymptotic consequence is proved in the
companion note. No primality search is attempted at enormous bit widths.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import isqrt
from pathlib import Path


COUNTS = Counter()


def ceildiv(a, b):
    return (a + b - 1) // b


def block_size(epsilon):
    threshold = ceildiv(128 * epsilon.denominator, epsilon.numerator)
    return 1 << (threshold - 1).bit_length()


def initial_parameters(n, epsilon):
    u = (n + 1).bit_length()
    delta = epsilon / 48
    radius = isqrt((n * delta.numerator) // (u * delta.denominator))
    if radius < 1:
        return u, block_size(epsilon), 0
    initial = 1 << (radius.bit_length() - 1)
    assert Fraction(initial * initial) <= delta * n / u
    assert Fraction(4 * initial * initial) > delta * n / u
    return u, block_size(epsilon), initial


def stage_parameters(n, epsilon, stage):
    u, block, initial = initial_parameters(n, epsilon)
    A = initial >> stage
    if A < 2 * block:
        return None
    P, T, D = A // block, A * A, A * A // 2
    t = D.bit_length() - 1
    w = 2 * ceildiv(n + 1, 2 * (T - A - 1))
    residue, carry = 2 * w + t + 1, w + t
    width = 6 * block * block * residue + 2 * block * residue + 4 * block * carry
    return {"u": u, "K": block, "A0": initial, "A": A, "P": P,
            "T": T, "D": D, "t": t, "w": w, "f": residue,
            "c": carry, "lambda": width}


def parameter_controls():
    epsilons = [Fraction(1), Fraction(1, 2), Fraction(1, 4), Fraction(1, 16),
                Fraction(1, 64), Fraction(3, 7)]
    worst = Fraction(0)
    max_width_scaled = Fraction(0)
    for epsilon in epsilons:
        cutoff = ceildiv((1 << 40) * epsilon.denominator, epsilon.numerator)
        inputs = {cutoff, cutoff + 1, cutoff + 2}
        for exponent in [40, 41, 45, 50, 60, 80, 120, 200, 400, 800]:
            inputs.update((1 << exponent) + offset for offset in [-2, -1, 0, 1, 2]
                          if (1 << exponent) + offset >= cutoff)
        # Directly approach several A0 jump thresholds, with u recomputed
        # at every candidate rather than treated as a real logarithm.
        for a_exponent in [10, 15, 20, 40, 100]:
            for u_guess in range(2 * a_exponent, 2 * a_exponent + 20):
                candidate = ceildiv((1 << (2 * a_exponent)) * 48 * u_guess
                                    * epsilon.denominator, epsilon.numerator)
                for offset in [-1, 0, 1]:
                    if candidate + offset >= cutoff:
                        inputs.add(candidate + offset)
        for n in sorted(inputs):
            u, block, initial = initial_parameters(n, epsilon)
            if initial < 2 * block:
                continue
            final_stage = (initial // (2 * block)).bit_length() - 1
            stages = {0, 1, 2, 3, final_stage // 2, max(0, final_stage - 1), final_stage}
            for stage in sorted(stages):
                p = stage_parameters(n, epsilon, stage)
                if p is None:
                    continue
                A, P, K, t, w = p["A"], p["P"], p["K"], p["t"], p["w"]
                assert A >= 2 * P and A >= 2 * K and A == P * K
                assert not (w & 1)
                assert w >= Fraction(n + 1, A * A - A - 1)
                assert w < Fraction(n + 1, A * A - A - 1) + 2
                assert Fraction(A * A * w, n) < 1 + epsilon / 100
                assert Fraction(A * A * t, n) <= epsilon / 48
                assert Fraction(A * A, n) <= epsilon / (48 * u)
                exact_volume = Fraction(P * P * p["lambda"], n)
                expanded_volume = Fraction(A * A, n) * (
                    (12 + Fraction(8, K)) * w + (6 + Fraction(6, K)) * t
                    + 6 + Fraction(2, K))
                assert exact_volume == expanded_volume
                assert exact_volume < 12 + epsilon
                lower_log_ratio = Fraction(2 * w, 2 * (A.bit_length() - 1))
                assert lower_log_ratio >= Fraction(2 * (4 ** stage), 1) / (epsilon / 48)
                assert Fraction(2 * (A - 1).bit_length()) <= p["f"]
                assert P * P * p["lambda"] >= 6 * A * A * (2 * (A.bit_length() - 1))
                worst = max(worst, (exact_volume - 12) / epsilon)
                max_width_scaled = max(max_width_scaled, (Fraction(A * A * w, n) - 1) / epsilon)
                COUNTS["exact_parameter_stage_controls"] += 1
            COUNTS["exact_input_boundary_controls"] += 1
    universal_margin = (Fraction(12, 100) + Fraction(1, 16) + Fraction(1, 1600)
                        + (6 + Fraction(3, 64)) / 48
                        + (6 + Fraction(1, 64)) / (48 * 41))
    assert universal_margin < 1
    return {"epsilon_values": [str(x) for x in epsilons],
            "uniform_proved_epsilon_margin": str(universal_margin),
            "largest_observed_volume_excess_divided_by_epsilon": str(worst),
            "largest_observed_digit_capacity_excess_divided_by_epsilon": str(max_width_scaled)}


def is_prime(value):
    if value < 2:
        return False
    for divisor in range(2, isqrt(value) + 1):
        if value % divisor == 0:
            return False
    return True


def dyadic_prime(order, exponent):
    lower = 1 << exponent
    first = lower + ((1 - lower) % order)
    if first <= lower:
        first += order
    for candidate in range(first, 2 * lower, order):
        if is_prime(candidate):
            return candidate
    return None


def prime_controls():
    for A in [4, 8, 16]:
        T, D = A * A, A * A // 2
        t = D.bit_length() - 1
        for w in [2, 4]:
            exponent = 2 * w + t
            prime = dyadic_prime(T, exponent)
            assert prime is not None
            assert (1 << exponent) < prime < (1 << (exponent + 1))
            assert (prime - 1) % T == 0
            assert prime > D * (1 << (2 * w))
            nonresidue = next(g for g in range(2, prime)
                              if pow(g, (prime - 1) // 2, prime) == prime - 1)
            root = pow(nonresidue, (prime - 1) // T, prime)
            assert pow(root, T, prime) == 1
            assert pow(root, T // 2, prime) == prime - 1
            COUNTS["small_exact_dyadic_prime_and_root_controls"] += 1
    # This is outside the even-w production schedule. It is a substantive
    # control against assuming every dyadic interval contains the needed
    # progression prime without an analytic size hypothesis.
    assert dyadic_prime(16, 5) is None
    COUNTS["below_hypothesis_empty_prime_interval_controls"] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    parameters = parameter_controls()
    prime_controls()
    result = {"status": "passed", "claim": "exact adaptive prime-constructor ceilings and volume",
              "scope": "Analytic prime-distribution constants are proved symbolically, not numerically estimated.",
              "parameters": parameters, "counts": dict(sorted(COUNTS.items()))}
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
