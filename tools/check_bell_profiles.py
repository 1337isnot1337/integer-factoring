#!/usr/bin/env python3
"""Independent affine-solver controls and exact all-length bound receipts."""

import argparse
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import random

from bell_crossing import rank_mod, transform
from search_bell_profiles import profiles
from search_bell_synthesis import AffineTernary


def solver_controls() -> int:
    rng = random.Random(338922)
    checks = 0
    for variables in range(1, 13):
        for _ in range(12):
            system = AffineTernary(variables)
            equations = []
            for _ in range(2 * variables + 4):
                equation = [rng.randrange(3) for _ in range(variables + 1)]
                trial = equations + [equation]
                feasible = rank_mod([r[:-1] for r in trial], 3) == rank_mod(trial, 3)
                one = sum(1 << j for j, x in enumerate(equation) if x == 1)
                two = sum(1 << j for j, x in enumerate(equation) if x == 2)
                before = system.basis.copy()
                assert system.add(one, two) == feasible
                if feasible:
                    equations = trial
                else:
                    assert system.basis == before
                solution = system.solution()
                assert all(sum(a * x for a, x in zip(row[:-1], solution)) % 3 == row[-1]
                           for row in equations)
                checks += 1
    return checks


def literal_ranks(coefficients: tuple[tuple[int, ...], ...], length: int) -> tuple[int, ...]:
    local = tuple(zip(*coefficients))
    degree = len(coefficients) - 1
    ranks = [0] * (degree * length + 1)
    for word in product(local, repeat=length):
        polynomial = [1]
        for factor in word:
            new = [0] * (len(polynomial) + degree)
            for i, a in enumerate(polynomial):
                for j, b in enumerate(factor):
                    new[i + j] = (new[i + j] + a * b) % 3
            polynomial = new
        for i, a in enumerate(polynomial):
            ranks[i] += bool(a)
    return tuple(ranks)


def check_record(data: dict) -> int:
    bits = data["bits"]
    saved = data["best"]
    coefficients = tuple(map(tuple, saved["coefficients_by_degree"]))
    # Do not run the optimizer in this readback: validate supplied rational
    # data and exact physical coefficients independently of search choices.
    other = tuple(transform(row, bits, 3) for row in coefficients)
    assert list(map(list, other)) == saved["transformed_coefficients_by_degree"]
    assert all(sum(col) % 3 == 1 for col in zip(*coefficients))
    assert all(sum(col) % 3 == 1 for col in zip(*other))
    rl, rh = profiles(coefficients)
    cl, ch = profiles(other)
    for key, value in (("row_low", rl), ("row_high", rh),
                       ("column_low", cl), ("column_high", ch)):
        assert list(value) == saved[key]
    row, col = literal_ranks(coefficients, 2), literal_ranks(other, 2)
    for bound in saved["bounds"]:
        if bound["orientation"] == "row-low/column-high":
            low, high, lower_ranks, upper_ranks = rl, ch, row, col
        else:
            low, high, lower_ranks, upper_ranks = cl, rh, col, row
        assert list(low) == bound["low_counts"]
        assert list(high) == bound["high_counts"]
        x = Fraction(*bound["scale"])
        assert x >= 1
        a = sum(Fraction(n) / x ** i for i, n in enumerate(low))
        b = sum(Fraction(n) * x ** i for i, n in enumerate(high))
        exact = a * b
        assert exact == Fraction(*bound["pressure"])
        for i, (r, c) in enumerate(zip(lower_ranks, upper_ranks)):
            assert r <= x ** i * a ** 2
            assert c <= b ** 2 / x ** i
            assert r * c <= exact ** 2
        numerator, denominator = bound["rho_upper_exact"]
        assert exact ** denominator <= 2 ** (4 * bits * numerator)
        assert bound["strict_below_critical_certified"] == (exact ** 577 < 2 ** (4 * bits * 408))
    return 2 * len(row)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("records", nargs="*")
    args = parser.parse_args()
    checks = solver_controls()
    ranks = 0
    for name in args.records:
        ranks += check_record(json.loads(Path(name).read_text()))
    assert 2 * 408 ** 2 < 577 ** 2
    print(f"Passed {checks} dense affine controls and {ranks} literal coefficient-rank bounds.")


if __name__ == "__main__":
    main()
