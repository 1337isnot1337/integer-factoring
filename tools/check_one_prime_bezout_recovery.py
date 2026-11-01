#!/usr/bin/env python3
"""Exact finite controls for the conditional one-field recovery wrapper.

The oracle below explicitly uses supplied finite divisor lists. These controls
do not implement or time a compressed normalized Bezout coefficient provider.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import prod
from pathlib import Path

from one_prime_bezout_recovery import (
    one_prime_recovery_plan, prime_above_bits, smallest_prime_from_one_field,
)
from prime_bezout_counting import is_prime

ROOT = Path(__file__).resolve().parents[1]


def literal_divisors(n):
    return [d for d in range(2, n // 2 + 1) if n % d == 0]


def oracle_from_divisors(n, divisors, requested_primes):
    def oracle(contract):
        prime = contract["prime"]
        requested_primes.add(prime)
        assert contract["n"] == n and contract["precision"] == 1
        a0, b0, count = (contract[key] for key in ("first_a", "first_b", "count"))
        return sum(a0 <= a <= a0 + prime * (count - 1)
                   and (a - a0) % prime == 0
                   and b0 <= n // a <= b0 + prime * (count - 1)
                   and (n // a - b0) % prime == 0 for a in divisors) % prime
    return oracle


def run_case(n, prime, divisors):
    requested = set()
    answer = smallest_prime_from_one_field(n, oracle_from_divisors(n, divisors, requested), prime)
    expected = min(divisors, default=n)
    assert answer["smallest_prime_factor"] == expected
    assert len(requested) <= 1
    if requested:
        assert requested == {answer["prime"]}
    assert answer["coefficient_queries"] == (answer["prime"] - 1) * answer["suffix_queries"]
    for interval in answer["prime_intervals"]:
        selected = [d for d in divisors if interval["lower"] <= d < interval["upper"]]
        assert all(is_prime(d) for d in selected)
        assert len(selected) < answer["prime"]
        assert len(selected) == interval["exact_prime_divisor_count"]
        assert all(d >= interval["lower"] for d in divisors if is_prime(d))
    return answer


def verify():
    cases = queries = used_provider = interval_checks = 0
    for n in range(2, 258):
        divisors = literal_divisors(n)
        for prime in (3, 5, 7, 11, 17, None):
            answer = run_case(n, prime, divisors)
            cases += 1
            queries += answer["coefficient_queries"]
            used_provider += answer["coefficient_queries"] > 0
            interval_checks += len(answer["prime_intervals"])

    controls = []
    for factors, prime in [
        ((11, 13, 17, 19), 7),
        ((101, 103), 17),
        ((101, 101), 17),
        ((97, 101, 103, 107, 109, 113), 13),
        ((1009,), 11),
    ]:
        assert all(is_prime(p) for p in factors)
        n = prod(factors)
        divisors = sorted({prod(factors[i] for i in range(len(factors)) if mask >> i & 1)
                           for mask in range(1 << len(factors))} - {1, n})
        answer = run_case(n, prime, divisors)
        controls.append({"n": n, "prime": prime,
                         "ordinary_proper_divisor_count_mod_prime": len(divisors) % prime,
                         "prime_intervals": answer["prime_intervals"],
                         "smallest_prime_factor": answer["smallest_prime_factor"],
                         "coefficient_queries": answer["coefficient_queries"]})
        cases += 1
        queries += answer["coefficient_queries"]
        used_provider += answer["coefficient_queries"] > 0
        interval_checks += len(answer["prime_intervals"])
    assert controls[0]["ordinary_proper_divisor_count_mod_prime"] == 0
    assert controls[0]["smallest_prime_factor"] == 11
    assert len(controls[1]["prime_intervals"]) >= 3
    assert controls[2]["smallest_prime_factor"] == 101

    for bits in range(1, 129):
        prime = prime_above_bits(bits)
        assert prime > bits and prime <= 16 * (bits + 1) ** 4 and is_prime(prime)

    large = one_prime_recovery_plan((1 << 1000) + 129, 101)
    assert large["input_bits"] == 1001
    assert large["trial_bits"] == 10
    assert large["factor_multiplicity_bound"] == 100 < large["prime"]
    default_large = one_prime_recovery_plan((1 << 400) + 87)
    assert default_large["prime"] > default_large["input_bits"]
    assert default_large["trial_bits"] == 1

    rejected = 0
    for args in [(1, 3), (15, 2), (15, 9), (15, -3)]:
        try:
            one_prime_recovery_plan(*args)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("Invalid contract was accepted")

    source_paths = ["tools/one_prime_bezout_recovery.py", "tools/check_one_prime_bezout_recovery.py",
                    "tools/prime_bezout_counting.py"]
    return {
        "theorem": "One prime field and certified prime-only intervals recover the smallest prime",
        "finite_wrapper_cases": cases, "cases_using_supplied_provider": used_provider,
        "finite_known_divisor_coefficient_queries": queries,
        "certified_prime_interval_controls": interval_checks,
        "prime_selection_controls": 128, "invalid_contracts_rejected": rejected,
        "cancellation_and_square_controls": controls,
        "large_contract_only": {
            "input_bits": large["input_bits"], "prime": large["prime"],
            "trial_bits": large["trial_bits"],
            "factor_multiplicity_bound": large["factor_multiplicity_bound"],
            "default_input_bits": default_large["input_bits"],
            "default_prime": default_large["prime"],
            "scope": "Binary parameter construction only; no large coefficient supplied or evaluated",
        },
        "source_sha256": {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
                          for path in source_paths},
        "scope": "The global normalized coefficient provider remains an explicit hypothesis; no improved factoring score is claimed.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    args = parser.parse_args()
    answer = json.dumps(verify(), indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(answer)
    print(answer, end="")


if __name__ == "__main__":
    main()
