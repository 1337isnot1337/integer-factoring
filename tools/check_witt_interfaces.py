#!/usr/bin/env python3
"""Fast exact replay of the positive Witt duality interface and finite jets."""

import argparse
from fractions import Fraction
import hashlib
import json
from math import factorial, prod
from pathlib import Path
from time import monotonic

from witt_artin_hasse import artin_hasse_coefficients, hall_gram, vacuum_map
from witt_exchange_jets import calculate
from weighted_energy_pencil import rank_two


def mobius(n):
    value, prime = 1, 2
    while prime * prime <= n:
        if n % prime == 0:
            n //= prime
            value = -value
            if n % prime == 0:
                return 0
        prime += 1
    return -value if n > 1 else value


def euler_coefficients(N):
    """Independent rational Euler product over odd indices."""
    answer = [Fraction(1)] + [Fraction(0)] * N
    for n in range(1, N + 1, 2):
        exponent = Fraction(mobius(n), n)
        local = [Fraction(1)]
        for k in range(1, N // n + 1):
            local.append(local[-1] * (exponent + k - 1) / k)
        updated = [Fraction(0)] * (N + 1)
        for i, a in enumerate(answer):
            for k, b in enumerate(local[:(N-i)//n+1]):
                updated[i + n*k] += a * b
        answer = updated
    return answer


def determinant(matrix):
    """Fraction-free Bareiss elimination, independent of parity ranks."""
    A = [row[:] for row in matrix]
    sign, denominator = 1, 1
    for k in range(len(A) - 1):
        if not A[k][k]:
            pivot = next((i for i in range(k+1, len(A)) if A[i][k]), None)
            if pivot is None:
                return 0
            A[k], A[pivot] = A[pivot], A[k]
            sign = -sign
        pivot = A[k][k]
        for i in range(k+1, len(A)):
            for j in range(k+1, len(A)):
                numerator = A[i][j] * pivot - A[i][k] * A[k][j]
                assert numerator % denominator == 0
                A[i][j] = numerator // denominator
            A[i][k] = 0
        denominator = pivot
    return sign * A[-1][-1]


def oddpart(n):
    return n // (n & -n)


def check_length_unitriangular(basis, rows):
    checks = 0
    for i, row in enumerate(rows):
        assert row >> i & 1
        row ^= 1 << i
        while row:
            bit = row & -row
            j = bit.bit_length() - 1
            assert sum(basis[j]) > sum(basis[i])
            row ^= bit
            checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-h", type=int, default=16)
    parser.add_argument("--determinant-max", type=int, default=12)
    parser.add_argument("--coefficient-limit", type=int, default=32)
    parser.add_argument("--output")
    args = parser.parse_args()
    started = monotonic()
    grades, triangular_entries, determinant_checks = [], 0, 0
    for h in range(1, args.max_h + 1):
        basis, gram, rows = hall_gram(h)
        assert rank_two(rows) == len(basis)
        vacuum_ranks = []
        for neutral in (False, True):
            other, vacuum = vacuum_map(h, neutral=neutral)
            assert other == basis
            actual = rank_two(vacuum)
            assert actual == len(basis)
            triangular_entries += check_length_unitriangular(basis, vacuum)
            vacuum_ranks.append(actual)
        if h <= args.determinant_max:
            expected = prod(oddpart(factorial(a)) for m in basis for a in m)
            assert determinant(gram) == expected
            determinant_checks += 1
        grades.append({"h": h, "dimension": len(basis),
                       "hall_gram_rank": len(basis),
                       "charged_vacuum_rank": vacuum_ranks[0],
                       "neutral_vacuum_rank": vacuum_ranks[1]})

    recurrence = artin_hasse_coefficients(args.coefficient_limit)
    assert recurrence == euler_coefficients(args.coefficient_limit)
    assert all(value.denominator & 1 for value in recurrence)

    jets = []
    for neutral in (False, True):
        result = calculate(16, neutral=neutral)
        assert result["h"] == 2
        assert rank_two(result["vacuum_gram_rows"]) == 2
        systems = result["finite_exchange_systems"]
        assert [item["angle_degree"] for item in systems] == list(range(13))
        assert all(item["solution"] is None for item in systems)
        jets.append({"neutral": neutral, "energy_window": 16,
                     "state_count": result["total_states"],
                     "vacuum_rank": 2,
                     "inconsistent_angle_degrees": list(range(13)),
                     "system_ranks": [item["rank"] for item in systems]})

    source_names = ["check_witt_interfaces.py", "witt_artin_hasse.py",
                    "witt_exchange_jets.py", "weighted_energy_pencil.py",
                    "check_modular_crossing.py"]
    directory = Path(__file__).resolve().parent
    report = {
        "source_sha256": {name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
                          for name in source_names},
        "checked_grades": grades,
        "full_rank_gram_checks": len(grades),
        "full_rank_vacuum_checks": 2 * len(grades),
        "length_unitriangular_maps": 2 * len(grades),
        "strict_length_increasing_entries": triangular_entries,
        "independent_integer_determinants": determinant_checks,
        "artin_hasse_euler_product_coefficient_checks": len(recurrence),
        "charged_and_neutral_exchange": jets,
        "scope": "Exact finite verification of the all-weight duality proofs; finite exchange failure only at h=2 and listed windows/degrees.",
        "elapsed_seconds": round(monotonic() - started, 3),
    }
    if args.output:
        Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items()
                      if key not in {"source_sha256", "checked_grades",
                                     "charged_and_neutral_exchange"}}, indent=2))


if __name__ == "__main__":
    main()
