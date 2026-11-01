#!/usr/bin/env python3
"""Independent exact controls for the Fermat-torus binomial circuit."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from check_frobenius_rational_trace import field_of_degree
from fermat_torus_binomial_contraction import torus_binomial_sum, dyadic_box_star_sum


def binary_power(value, exponent, multiply, one=1):
    result = one
    while exponent:
        if exponent & 1:
            result = multiply(result, value)
        value = multiply(value, value)
        exponent >>= 1
    return result


def reference_sum(field, m, factors, monomial=0, scalar=1):
    modulus = (1 << m) + 1
    torus = [value for value in range(1, field.size)
             if field.power(value, modulus) == 1]
    assert len(torus) == modulus
    result = 0
    for value in torus:
        term = field.mul(scalar, field.power(value, monomial % modulus))
        for a, b, exponent in factors:
            term = field.mul(term, a ^ field.mul(
                b, field.power(value, pow(2, exponent, modulus))))
        result ^= term
    return result, len(torus)


def subset_reference(field, m, factors, monomial, scalar):
    """Full monomial expansion; no binary carry algorithm is used."""
    modulus = (1 << m) + 1
    result = 0
    for choices in product((0, 1), repeat=len(factors)):
        exponent, coefficient = monomial, scalar
        for choice, (a, b, power) in zip(choices, factors):
            coefficient = field.mul(coefficient, b if choice else a)
            if choice:
                exponent += pow(2, power, modulus)
        if exponent % modulus == 0:
            result ^= coefficient
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(0xFE2A72026)
    counts = Counter()
    models = []
    for m in range(1, 6):
        field = field_of_degree(2 * m)
        modulus = (1 << m) + 1
        models.append({"m": m, "torus_order": modulus,
                       "field_degree": field.degree,
                       "irreducible_polynomial": field.polynomial})
        for trial in range(60):
            factors = [(rng.randrange(field.size), rng.randrange(field.size),
                        rng.randrange(12 * m + 1)) for _ in range(trial % 15)]
            monomial = rng.randrange(-8 * modulus, 8 * modulus + 1)
            scalar = rng.randrange(field.size)
            actual = torus_binomial_sum(
                m, factors, multiply=field.mul, monomial=monomial, scalar=scalar)
            expected, points = reference_sum(field, m, factors, monomial, scalar)
            assert actual == expected
            counts["complete_torus_sums"] += 1
            counts["enumerated_reference_points"] += points
            counts["zero_coefficient_factors"] += sum(not a or not b for a, b, _ in factors)
            counts["negative_reduced_exponents"] += sum(e % (2*m) >= m for _, _, e in factors)
        torus = [value for value in range(1, field.size)
                 if field.power(value, modulus) == 1]
        for trial in range(24):
            boxes = [(rng.randrange(field.size), rng.randrange(6 * m + 1),
                      rng.randrange(6)) for _ in range(trial % 4)]
            monomial = rng.randrange(-modulus, modulus + 1)
            actual = dyadic_box_star_sum(m, boxes, multiply=field.mul, monomial=monomial)
            expected = 0
            for value in torus:
                term = field.power(value, monomial % modulus)
                for gamma, exponent, width in boxes:
                    ratio = field.mul(gamma, field.power(value, pow(2, exponent, modulus)))
                    subtotal, power = 0, 1
                    for _ in range(1 << width):
                        subtotal ^= power
                        power = field.mul(power, ratio)
                    term = field.mul(term, subtotal)
                expected ^= term
            assert actual == expected
            counts["complete_dyadic_box_stars"] += 1
            counts["dyadic_reference_points"] += len(torus)

        # Coefficients in E[epsilon]/epsilon^2 include genuine zero divisors.
        bits, mask = field.degree, field.size - 1
        def dual_multiply(first, second):
            a, b = first & mask, first >> bits
            c, d = second & mask, second >> bits
            return field.mul(a, c) | ((field.mul(a, d) ^ field.mul(b, c)) << bits)
        for trial in range(16):
            factors = [(rng.randrange(field.size**2), rng.randrange(field.size**2),
                        rng.randrange(6 * m + 1)) for _ in range(1 + trial % 6)]
            if trial % 4 == 0:
                a, _, exponent = factors[0]
                factors[0] = (a, 1 << bits, exponent)
            monomial = rng.randrange(-modulus, modulus + 1)
            actual = torus_binomial_sum(m, factors, multiply=dual_multiply, monomial=monomial)
            expected = 0
            for value in torus:
                term = field.power(value, monomial % modulus)
                for a, b, exponent in factors:
                    term = dual_multiply(term, a ^ dual_multiply(
                        b, field.power(value, pow(2, exponent, modulus))))
                expected ^= term
            assert actual == expected
            counts["zero_divisor_algebra_sums"] += 1
            counts["zero_divisor_reference_points"] += len(torus)

    # Large m controls use complete monomial expansion with small r, so
    # neither the verifier nor the tested circuit enumerates its torus.
    field = field_of_degree(2)
    for m in (20, 37, 64, 127):
        for trial in range(8):
            factors = [(rng.randrange(4), rng.randrange(4), rng.randrange(10 * m))
                       for _ in range(6 + trial % 7)]
            monomial = rng.randrange(-1 << (m + 4), 1 << (m + 4))
            scalar = rng.randrange(4)
            actual = torus_binomial_sum(
                m, factors, multiply=field.mul, monomial=monomial, scalar=scalar)
            assert actual == subset_reference(field, m, factors, monomial, scalar)
            counts["large_bitlength_complete_expansions"] += 1
            counts["large_bitlength_expanded_monomials"] += 1 << len(factors)

    # The arithmetic gate schedule depends only on exponent data.
    scheduled_counts = []
    for coefficient_seed in range(4):
        calls = Counter()
        def counted_multiply(a, b):
            calls["multiplications"] += 1
            return field.mul(a, b)
        def counted_add(a, b):
            calls["additions"] += 1
            return a ^ b
        factors = [(coefficient_seed, (coefficient_seed + i) % 4, exponent)
                   for i, exponent in enumerate((0, 1, 3, 4, 7, 11))]
        torus_binomial_sum(4, factors, multiply=counted_multiply, add=counted_add, monomial=-19)
        scheduled_counts.append(dict(calls))
    assert all(row == scheduled_counts[0] for row in scheduled_counts)
    counts["coefficient_independent_gate_schedules"] = len(scheduled_counts)

    names = (Path(__file__).name, "fermat_torus_binomial_contraction.py",
             "check_frobenius_rational_trace.py", "check_carry_fourier_elimination.py")
    report = {
        "status": "passed",
        "scope": "complete finite torus sums, literal dyadic boxes, coefficient algebras with zero divisors, and large binary-exponent monomial controls",
        "counts": dict(counts),
        "models": models,
        "fixed_gate_counts": scheduled_counts[0],
        "source_sha256": {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names},
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "source_sha256"}))


if __name__ == "__main__":
    main()
