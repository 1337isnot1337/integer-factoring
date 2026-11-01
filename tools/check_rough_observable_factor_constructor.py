#!/usr/bin/env python3
"""Exact controls for pretesting into directly observable operand lengths."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import isqrt
from pathlib import Path
import random

from check_observable_boundary_constructor import BoundaryNetwork, fa_message, interpolate_input
from check_prime_transform_boundary_constructor import PrimeTransform, PrimeLengthTransform
from check_width_sensitive_prime_constructor import initial_parameters, stage_parameters


COUNTS = Counter()


def divisors(number):
    result = set()
    for candidate in range(1, isqrt(number) + 1):
        if number % candidate == 0:
            result.add(candidate)
            result.add(number // candidate)
    return sorted(result)


def parity_free_support_controls():
    for number in range(2, 4097):
        n = number.bit_length()
        all_divisors = divisors(number)
        least = all_divisors[1]
        for u in range(1, 7):
            small = next((d for d in all_divisors[1:] if d <= min(1 << u, isqrt(number))), None)
            if small is not None:
                assert small == least
                COUNTS["small_factor_pretest_returns"] += 1
                continue
            if (1 << u) >= isqrt(number):
                assert least == number
                COUNTS["complete_trial_prime_returns"] += 1
                continue
            retained = [d for d in all_divisors if u + 1 <= d.bit_length() <= (n + 1) // 2]
            assert bool(retained) == (least != number)
            assert len(all_divisors) <= 1 << ((n - 1) // u)
            COUNTS["rough_support_and_completeness_controls"] += 1
            if not retained:
                continue
            prefix = 0
            for remaining in range(n - 1, -1, -1):
                left = prefix << 1
                if any(d >> remaining == left for d in retained):
                    prefix = left
                else:
                    prefix = left + 1
                    assert any(d >> remaining == prefix for d in retained)
                COUNTS["zero_padded_prefix_decisions"] += 1
            assert prefix == least
            COUNTS["least_prime_prefix_recoveries"] += 1


def height_controls():
    for A in [4, 8, 16, 32, 64]:
        for P in range(2, A // 2 + 1):
            K = (A + P - 1) // P
            if K > A // 2:
                continue
            for w in range(1, 7):
                cutoff = A * (K - 1) * w
                maximum_n = w * (A * A - A - 1) - 1
                for n in sorted({2 * cutoff + 1, maximum_n - 1, maximum_n}):
                    if not (2 * cutoff < n <= maximum_n):
                        continue
                    low, high = cutoff + 1, (n + 1) // 2
                    for ell in sorted({low, high, (low + high) // 2}):
                        if ell > high:
                            continue
                        dp, dq = (ell + w - 1) // w, (n - ell + w) // w
                        assert dp <= A * (A - K) and dq <= A * (A - K)
                        assert (dp + A - 1) // A <= A - K
                        assert (dq + A - 1) // A <= A - K
                        assert Fraction(n, cutoff) < 2 * P
                        COUNTS["exact_small_height_boundaries"] += 1
    for epsilon in [Fraction(1), Fraction(1, 2), Fraction(1, 16), Fraction(3, 7)]:
        for exponent in [48, 64, 128, 256, 512, 1024]:
            n = 1 << exponent
            _, K, initial = initial_parameters(n, epsilon)
            if initial < 2 * K:
                continue
            last = (initial // (2 * K)).bit_length() - 1
            for stage in sorted({0, 1, 2, last // 2, last}):
                parameter = stage_parameters(n, epsilon, stage)
                if parameter is None:
                    continue
                A, P, w = parameter["A"], parameter["P"], parameter["w"]
                cutoff = A * (K - 1) * w
                assert cutoff * P < A * A * w
                assert 12 * cutoff < P * parameter["lambda"]
                assert Fraction(n, cutoff) < 2 * P
                if 2 * cutoff < n:
                    for ell in [cutoff + 1, (n + 1) // 2]:
                        dp, dq = (ell + w - 1) // w, (n - ell + w) // w
                        assert max(dp, dq) <= A * (A - K)
                        COUNTS["exact_large_height_boundaries"] += 1
                COUNTS["exact_pretest_budget_controls"] += 1


def direct_interpolation_controls():
    randomizer = random.Random(0x0B5E8E)
    for model, block in [(PrimeTransform(2, 1), (0, 1)),
                         (PrimeTransform(4, 1), (2, 3)),
                         (PrimeTransform(4, 2), (0, 1, 2, 3))]:
        recipients = tuple(j for j in range(model.a) if j not in block)
        for height in range(len(recipients) + 1):
            for residue in block:
                for _ in range(4):
                    vector = tuple(randomizer.randrange(model.base) for _ in range(height))
                    observations = tuple(fa_message(model, residue, vector, j) for j in recipients)
                    decoded = tuple(interpolate_input(model, residue, height, recipients, observations))
                    assert decoded == (vector,)
                    assert max(0, height - len(recipients)) == 0
                    COUNTS["zero_free_digit_interpolation_controls"] += 1
    for left, right in [(17, 19), (19, 23), (29, 31), (257, 263), (719, 1451)]:
        number, ell = left * right, left.bit_length()
        model = PrimeLengthTransform(2, number.bit_length(), ell)
        cutoff = model.a * model.digit_bits
        assert ell >= cutoff + 1 and ell <= (number.bit_length() + 1) // 2
        assert not any(number % d == 0 for d in range(2, (1 << cutoff) + 1))
        for cut in [0, number.bit_length() // 2, number.bit_length()]:
            network = BoundaryNetwork(model, 2, number, cut, left >> cut)
            trace = network.trace(left, right)
            for group, block in enumerate(network.blocks):
                assert all(network.heights[family, residue] <= model.a - len(block)
                           for family in range(2) for residue in block)
                boundary = network.boundary(group, trace)
                actual, trials = network.entry(group, boundary, True)
                expected, reference_trials = network.entry(group, boundary, False)
                assert actual == expected == 1 and trials == 1
                COUNTS["rough_prime_network_entry_controls"] += 1
                COUNTS["reference_local_assignments"] += reference_trials
                assert network.entry(group, boundary, True, padding=True)[0] == 0
                COUNTS["invalid_padding_rejections"] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    parity_free_support_controls()
    height_controls()
    direct_interpolation_controls()
    result = {"status": "passed", "claim": "rough factor search with directly observable Boolean entries",
              "counts": dict(sorted(COUNTS.items())),
              "scope": "No free digit, norm bound, or lattice algorithm is used in the direct decoder."}
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
