#!/usr/bin/env python3
"""Independent exact controls for quadratic-parity Walsh Bell paths."""

import argparse
from fractions import Fraction
from itertools import product
import json
from math import comb
from pathlib import Path
import random

from bell_crossing import physical_matrix, rank_mod
from bell_walsh_constructions import (
    coefficients, energy_polynomials, folded_energies, folded_retiming,
    krawtchouk_matrix,
    norm_energies, parity, pressure, quadratic_two_level, radial_energies,
    radial_histogram, receipt, reciprocal_retime, shifted_dual, validate_energies,
    valuation_histogram,
)


def convolve(a, b, prime):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % prime
    return out


def covector(dimension, prime=None):
    a, remainder = divmod(dimension, 4)
    result = [0] * (dimension + 1)
    for j in range(a + 1):
        for s in range(remainder + 1):
            result[a + 2 * j + s] += comb(a, j) * comb(remainder, s)
    if prime:
        result = [x % prime for x in result]
    return result


def radial_transform(values, dimension, prime):
    k = krawtchouk_matrix(dimension, prime)
    first = [sum(x * y for x, y in zip(row, values)) % prime for row in k]
    signed = [(-1) ** ((w + 1) // 2) * value % prime
              for w, value in enumerate(first)]
    scale = pow(1 << dimension, -1, prime)
    return [scale * sum(x * y for x, y in zip(row, signed)) % prime for row in k]


def polynomial_remainder(value, divisor):
    while value and value.bit_length() >= divisor.bit_length():
        value ^= divisor << (value.bit_length() - divisor.bit_length())
    return value


def irreducible_binary(polynomial, degree):
    if polynomial.bit_length() != degree + 1:
        return False
    return all(polynomial_remainder(polynomial, divisor)
               for d in range(1, degree // 2 + 1)
               for divisor in range((1 << d) | 1, 1 << (d + 1), 2))


def check(saved_paths=()):
    rng = random.Random(2026110107)
    counts = {
        "independent_energy_polynomials": 0,
        "literal_walsh_coefficient_comparisons": 0,
        "exact_endpoint_and_reciprocity_coefficients": 0,
        "physical_cut_rank_controls": 0,
        "radial_histogram_controls": 0,
        "invariant_covector_columns": 0,
        "detected_piece_rank_bounds": 0,
        "literal_tensor_coefficient_bounds": 0,
        "saved_search_receipts": 0,
        "norm_constructors": 0,
        "retimed_reciprocity_coefficients": 0,
    }
    for prime in (3, 5, 7):
        for degree in range(1, 13):
            actual = energy_polynomials(degree, prime)
            for energy in range(degree + 1):
                independent = [1]
                for _ in range(degree - energy):
                    independent = convolve(independent, [1, 1], prime)
                for _ in range(energy):
                    independent = convolve(independent, [1, -1], prime)
                assert actual[energy] == tuple(independent)
                counts["independent_energy_polynomials"] += 1

    for bits in (1, 2, 3):
        base, shift = folded_energies(bits)
        q = len(base)
        sources = [base, tuple(3 * e for e in base)]
        sources.append(tuple(e + (2 if parity(v & 3) else 0)
                             for v, e in enumerate(base)))
        sources.append(tuple(0 if v == 0 else (e & 1) + 2 * rng.randint(1, 3)
                             for v, e in enumerate(base)))
        for prime in (3, 5, 7):
            for energies in sources:
                rows = coefficients(energies, bits, prime, shift)
                degree = len(rows) - 1
                polys = energy_polynomials(degree, prime)
                scale = pow(pow(2, degree, prime), -1, prime)
                for j, row in enumerate(rows):
                    literal = tuple(scale * sum((-1) ** parity(v & u) * polys[e][j]
                                                for u, e in enumerate(energies)) % prime
                                    for v in range(q))
                    assert row == literal
                    counts["literal_walsh_coefficient_comparisons"] += q
                    assert shifted_dual(row, bits, prime, shift) == rows[degree - j]
                    counts["exact_endpoint_and_reciprocity_coefficients"] += q
                assert all(sum(row[v] for row in rows) % prime == 1 for v in range(q))
                if bits <= 2:
                    for j in range(degree + 1):
                        physical = tuple(rows[j][v ^ shift] for v in range(q))
                        for column, expected in ((False, rows[j]), (True, rows[degree - j])):
                            assert rank_mod(physical_matrix(physical, bits, prime, column), prime) == sum(bool(x) for x in expected)
                            counts["physical_cut_rank_controls"] += 1

    # The independent radial computation uses only m+1 Krawtchouk classes.
    for bits in (1, 2, 3, 4):
        for prime in (3, 5, 7):
            for _ in range(8):
                levels = (0,) + tuple(((w + 1) // 2 & 1) + 2 * rng.randint(1, 4)
                                      for w in range(1, 2 * bits + 1))
                energies, shift = radial_energies(levels, bits)
                assert radial_histogram(levels, bits, prime) == valuation_histogram(coefficients(energies, bits, prime, shift))
                counts["radial_histogram_controls"] += 1

    # Full dual-action identities, including characteristics that divide
    # many of the orbit multiplicities. The multiplicities below stay integers.
    for dimension in range(2, 34, 2):
        a = dimension // 4
        minimum_rank = comb(dimension, a)
        for prime in (3, 5, 7):
            p = covector(dimension, prime)
            assert sum(p) % prime == pow(2, dimension - 3 * a, prime)
            assert all(p[j] == 0 for j in range(a))
            assert all(p[j] == 0 for j in range(dimension - a + 1, dimension + 1))
            for j in range(dimension + 1):
                unit = [int(i == j) for i in range(dimension + 1)]
                other = radial_transform(unit, dimension, prime)
                assert sum(x * y for x, y in zip(p, other)) % prime == p[j]
                counts["invariant_covector_columns"] += 1
            for _ in range(30):
                values = [rng.randrange(prime) for _ in range(dimension + 1)]
                if sum(x * y for x, y in zip(p, values)) % prime == 0:
                    continue
                other = radial_transform(values, dimension, prime)
                row_rank = sum(comb(dimension, j) for j, x in enumerate(values) if x)
                col_rank = sum(comb(dimension, j) for j, x in enumerate(other) if x)
                assert row_rank >= minimum_rank and col_rank >= minimum_rank
                counts["detected_piece_rank_bounds"] += 1

    # A two-level example cannot be a positive-code quotient at bits >= 3:
    # its degree is two, below the code quotient lower bound degree >= bits.
    two_level = []
    for bits in range(1, 7):
        result = receipt(quadratic_two_level(bits), bits)
        q, d = 4 ** bits, 2 ** bits
        if bits % 2 == 0:
            expected = [1, (q - d) // 2, (q + d) // 2 - 1]
        else:
            expected = [1, (q + d) // 2 - 1, (q - d) // 2]
        assert result["valuation_histogram"] == expected
        two_level.append({key: result[key] for key in
                          ("physical_bits", "degree", "valuation_histogram", "rho_pressure")})

    # Independently verify the moduli of the genuinely coupled norm examples.
    norm_profiles = []
    for bits, modulus in ((1, 3), (2, 7), (3, 11), (4, 19), (5, 37), (6, 67)):
        assert irreducible_binary(modulus, bits)
        energies, shift = norm_energies(bits, modulus)
        assert max(energies) == bits
        for prime in (3, 5, 7):
            result = receipt(energies, bits, prime, shift)
            counts["norm_constructors"] += 1
            norm_profiles.append({key: result[key] for key in
                                  ("physical_bits", "prime", "degree", "valuation_histogram", "rho_pressure")})

    # Retiming retains each old coefficient and preserves the full endpoint
    # and reciprocal identities. The new pressure uses the actual new orders.
    retimings = []
    for bits in range(1, 5):
        energies, shift = folded_energies(bits)
        rows = coefficients(energies, bits, 3, shift)
        result = folded_retiming(bits)
        retimed = reciprocal_retime(rows, result["retimed_degrees"])
        for j, row in enumerate(retimed):
            assert shifted_dual(row, bits, 3, shift) == retimed[-1 - j]
            counts["retimed_reciprocity_coefficients"] += len(row)
        assert all(sum(row[v] for row in retimed) % 3 == 1
                   for v in range(len(retimed[0])))
        histogram = valuation_histogram(retimed)
        bound = pressure(histogram, Fraction(*result["scale"]))
        assert [bound.numerator, bound.denominator] == result["pressure"]
        retimings.append(result)

    # A strict all-length improvement of the fixed physical-dimension-16
    # folded seed. Both separating inequalities are integer comparisons.
    retimed16 = retimings[-1]
    assert retimed16["retimed_degrees"] == [0, 22, 40, 58, 80]
    new_bound = Fraction(*retimed16["pressure"])
    assert new_bound ** 75 < 2 ** 1024
    scale = Fraction(13, 4)
    old_root = scale ** 2 + 9 * scale + 36 + 84 / scale + 126 / scale ** 2
    derivative = 2 * scale + 9 - 84 / scale ** 2 - 252 / scale ** 3
    old_lower = old_root - derivative ** 2 / 4
    assert old_lower == Fraction(2192857863, 19307236)
    assert old_lower > Fraction(2271, 20)
    assert Fraction(2271, 20) ** 75 > 2 ** 512

    # Literal length-two scalar products give the exact diagonal ranks.
    energies, shift = radial_energies((0, 1, 3, 2, 4), 2)
    rows = coefficients(energies, 2, 3, shift)
    degree = len(rows) - 1
    polynomials = [list(column) for column in zip(*rows)]
    ranks = [0] * (2 * degree + 1)
    for left, right in product(polynomials, repeat=2):
        for j, value in enumerate(convolve(left, right, 3)):
            ranks[j] += bool(value)
    histogram = valuation_histogram(rows)
    for scale in (Fraction(1), Fraction(3, 2), Fraction(3), Fraction(5)):
        bound = pressure(histogram, scale) ** 2
        for j in range(len(ranks)):
            assert ranks[j] * ranks[-1 - j] <= bound
            counts["literal_tensor_coefficient_bounds"] += 1

    for filename in saved_paths:
        saved = json.loads(Path(filename).read_text())
        rebuilt = receipt(tuple(saved["energies"]), saved["physical_bits"],
                          saved["prime"], saved["quadratic_shift"])
        for key in ("degree", "valuation_histogram", "coefficients_by_degree"):
            assert rebuilt[key] == saved[key]
        scale = Fraction(*saved["scale"])
        bound = pressure(saved["valuation_histogram"], scale)
        assert [bound.numerator, bound.denominator] == saved["pressure"]
        assert saved["strict_sufficient_target"] == (bound ** 577 < 2 ** (4 * saved["physical_bits"] * 408))
        assert saved["search"]["steps_completed"] > 0
        counts["saved_search_receipts"] += 1

    # Exact negative controls for the two independent construction hypotheses.
    base, shift = folded_energies(2)
    for index, change in ((1, 1), (1, -base[1])):
        invalid = list(base)
        invalid[index] += change
        try:
            validate_energies(invalid, 2, shift)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid spectral energy accepted")

    assert 3 ** 5 < 2 ** 8  # H_2(1/4) > 4/5 > 1/sqrt(2).
    return {
        "status": "passed", "counts": counts,
        "radial_rank_bound": "both ranks >= binom(m,floor(m/4)) in one common summand",
        "radial_optimal_leading_parameter": "H_2(1/4), including arbitrary recombinations in the full radial Bell space",
        "two_level_profiles": two_level,
        "norm_profiles": norm_profiles,
        "folded_retimings": retimings,
        "dimension_16_retiming": {
            "strict_separation": "rho_retimed < 64/75 < rho_original",
            "original_pressure_root_lower_bound": [old_lower.numerator, old_lower.denominator],
            "proof": "positive second derivative >= 2 and exact rational power comparisons",
        },
        "factoring_score": "1/2; unchanged",
        "scope": "ordinary proof plus exact finite controls; no new Lean theorem or global nonradial search optimum",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipts", nargs="*")
    parser.add_argument("--output")
    args = parser.parse_args()
    result = check(args.receipts)
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("two_level_profiles", "norm_profiles", "folded_retimings")}, indent=2))


if __name__ == "__main__":
    main()
