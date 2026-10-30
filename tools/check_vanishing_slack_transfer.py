#!/usr/bin/env python3
"""Exact finite checks for the vanishing-slack score transfer.

This checks recurrence and arithmetic statements, not new crossing tensors
and not the infinite-input factoring theorem. No floating-point arithmetic
or third-party packages are used.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path


def ceil_fraction(value: Fraction) -> int:
    return -(-value.numerator // value.denominator)


def schedule_log(j: int, x0: int = 18) -> int:
    return (3 * x0 + 6) * (1 << j) - 6 * j - 2 * x0 - 6


def slow_error(label_length: int) -> Fraction:
    """A monotone error of order 1/sqrt(log log lambda)."""
    assert label_length >= 1
    log_lambda = max(1, label_length.bit_length() - 1)
    log_log_lambda = max(1, log_lambda.bit_length() - 1)
    return Fraction(1, 16 + isqrt(log_log_lambda))


def check_schedule() -> dict[str, int]:
    cases = 0
    for x0 in (18, 19, 25, 100):
        previous = x0
        for j in range(257):
            current = schedule_log(j, x0)
            assert current == previous
            next_value = 2 * (3 * j + x0 + current)
            assert current < next_value <= 10 * current
            assert current >= 18 * (j + 1)
            previous = next_value
            cases += 1
    return {"exact_schedule_cases": cases}


def check_rounded_path(
    ell: int, log_n: int, rho: Fraction, use_cutoff_volume: bool
) -> dict[str, int]:
    """Check one rounded worst-envelope path at the exact source schedule."""
    schedule = [1 << schedule_log(j) for j in range(ell + 1)]
    p = schedule[-1]
    n = 1 << log_n
    sqrt_n = isqrt(n)
    ceil_sqrt_n = sqrt_n + (sqrt_n * sqrt_n != n)
    assert p**40 <= n < (1 << (40 * schedule_log(ell + 1)))
    assert p**400 > n
    assert p**40 <= n

    # This power of two is at most n^(9/20), so it is a safe tail start.
    lower_label = 1 << (9 * log_n // 20)
    assert lower_label**20 <= n**9
    rounding_envelope = Fraction(p, sqrt_n)
    beta = rho + slow_error(lower_label) + rounding_envelope
    assert Fraction(1, 2) <= rho <= beta < 1

    # A rational replacement for exp(3/(16rho)).  Since the sum of the
    # nonnegative relative schedule errors is <= 3/(16rho) < 1, every
    # finite product is at most the corresponding geometric series.
    d_bound = 1 / (1 - Fraction(3, 16) / rho)

    top_label = (
        ceil_sqrt_n // (p * p) + 1 if use_cutoff_volume else n // (p * p)
    )
    top_volume = p * p * top_label
    assert top_volume <= n
    current_volume = top_volume
    visited: list[tuple[int, int]] = []
    query_sum = Fraction(0)
    entropy_sum = Fraction(0)
    grouping_log_sum = 0
    expansions = 0

    for j in range(ell, -1, -1):
        visited.append((j, current_volume))
        assert current_volume <= d_bound * top_volume * beta ** (ell - j)
        if j == 0 or current_volume <= ceil_sqrt_n:
            break

        label_length, remainder = divmod(current_volume, schedule[j] ** 2)
        assert remainder == 0 and label_length >= 1
        assert label_length**20 >= n**9
        assert label_length >= lower_label
        local_error = slow_error(label_length)
        assert local_error <= slow_error(lower_label)
        k = schedule[j - 1]
        assert k * k <= p
        assert Fraction(k * k, current_volume) <= rounding_envelope

        zeta = Fraction(1, 1 << (j - 1 + 5))
        proposed = (rho + local_error + 3 * zeta) * current_volume
        child_label = ceil_fraction(proposed / (k * k))
        child_volume = child_label * k * k
        assert proposed <= child_volume < proposed + k * k
        assert child_volume <= (beta + 3 * zeta) * current_volume

        query_sum += Fraction(current_volume, k)
        # An integer entropy envelope with delta=1.  The +1 rounding
        # term is retained rather than silently discarded.
        denominator = label_length.bit_length()
        entropy = -(-label_length // denominator)
        entropy_sum += Fraction(schedule[j] ** 2 * entropy, 4)
        assert Fraction(entropy, label_length) <= (
            Fraction(1, lower_label.bit_length()) + Fraction(1, lower_label)
        )

        group_count_log = (1 << (3 * (j - 1) + 18 + schedule_log(j - 1)))
        group_count_log *= schedule_log(j - 1)
        assert group_count_log <= (isqrt(p) + 1) * schedule_log(ell)
        grouping_log_sum += group_count_log
        expansions += 1
        current_volume = child_volume

    volume_sum = sum(v for _, v in visited)
    geometric_bound = d_bound * top_volume / (1 - beta)
    assert volume_sum <= geometric_bound
    assert query_sum <= (
        d_bound * top_volume * beta**ell / ((1 << 17) - 1)
    )
    entropy_tail = Fraction(1, lower_label.bit_length()) + Fraction(1, lower_label)
    assert entropy_sum <= entropy_tail * geometric_bound / 4
    assert grouping_log_sum <= ell * (isqrt(p) + 1) * schedule_log(ell)

    stopped_level = visited[-1][0]
    if stopped_level > 0:
        assert current_volume <= ceil_sqrt_n
    else:
        assert current_volume <= d_bound * top_volume * beta**ell
    return {
        "expanded_calls": expansions,
        "early_cutoff": int(stopped_level > 0),
        "log2_input_length": log_n,
    }


def check_grid_paths() -> dict[str, int]:
    paths = expansions = early = max_log_n = 0
    for ell in range(5):
        endpoints = (40 * schedule_log(ell), 40 * schedule_log(ell + 1) - 1)
        for log_n in endpoints:
            for rho in (Fraction(1, 2), Fraction(2, 3), Fraction(7, 10)):
                for cutoff in (False, True):
                    record = check_rounded_path(ell, log_n, rho, cutoff)
                    paths += 1
                    expansions += record["expanded_calls"]
                    early += record["early_cutoff"]
                    max_log_n = max(max_log_n, record["log2_input_length"])
    return {
        "rounded_grid_paths": paths,
        "expanded_grid_calls": expansions,
        "paths_reaching_early_cutoff": early,
        "largest_log2_of_input_bit_length": max_log_n,
    }


def check_series() -> dict[str, int]:
    cases = 0
    for rho in (Fraction(1, 2), Fraction(2, 3), Fraction(7, 10), Fraction(7, 8)):
        for length in range(1, 257):
            # K_(j-1) >= 2^(18j) reduces the infinite query series to
            # the geometric series checked here.  Check finite sums
            # with their exact tails, not only finite partial sums.
            ratio = Fraction(1, (1 << 18)) / rho
            partial = ratio * (1 - ratio**length) / (1 - ratio)
            tail = ratio ** (length + 1) / (1 - ratio)
            assert partial + tail == ratio / (1 - ratio)
            assert partial + tail <= Fraction(1, (1 << 17) - 1)
            cases += 1
    return {"query_series_with_exact_tail_cases": cases}


def check_score_not_constant() -> dict[str, int]:
    divergence_cases = score_cases = 0
    for rho in (Fraction(1, 2), Fraction(2, 3), Fraction(7, 10)):
        for m in (8, 16, 32, 64, 128):
            beta = rho + Fraction(1, m)
            assert beta < 1
            ratio = (beta / rho) ** (m * m)
            # Bernoulli proves an unbounded ratio, while beta -> rho.
            assert ratio >= 1 + Fraction(m, 1) / rho
            divergence_cases += 1
        for q in (1, 2, 4, 8, 16, 32, 64, 128):
            m = 8 * q
            beta = rho + Fraction(1, m)
            # log2(beta/rho) < 1/q, certified with integer arithmetic.
            assert beta**q < 2 * rho**q
            score_cases += 1
    return {
        "unbounded_constant_ratio_cases": divergence_cases,
        "vanishing_score_loss_certificates": score_cases,
    }


def check_rough_support() -> dict[str, int]:
    # The identity sum(a_i+1) is not the divisor count; verify the
    # required product inequality directly on independent exponent
    # vectors.  The prime-size premise itself is a theorem assumption.
    from itertools import product

    cases = 0
    for length in range(7):
        for exponents in product(range(1, 6), repeat=length):
            omega = sum(exponents)
            divisors = 1
            for exponent in exponents:
                divisors *= exponent + 1
            assert divisors <= 1 << omega
            cases += 1
    for n in range(2, 8193):
        b = isqrt(n)
        b += b * b != n
        assert -(-n // b) <= b
        cases += 1
    return {"rough_support_and_threshold_cases": cases}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results: dict[str, object] = {
        "status": "passed",
        "scope": "Exact finite recurrence checks; no new crossing family or improved factoring exponent.",
        "arithmetic": "Python integers and fractions only",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    for routine in (
        check_schedule,
        check_grid_paths,
        check_series,
        check_score_not_constant,
        check_rough_support,
    ):
        results.update(routine())
    encoded = json.dumps(results, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
