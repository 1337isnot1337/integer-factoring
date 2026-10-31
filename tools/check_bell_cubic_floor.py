#!/usr/bin/env python3
"""Exact controls for the outer-rank-one cubic Bell pressure restriction."""

import argparse
from fractions import Fraction
import json
from math import comb
import random

from bell_crossing import transform


def check():
    rng = random.Random(20261101)
    endpoint_cases = pressure_comparisons = 0
    quartic_cases = quartic_comparisons = 0
    for bits in (1, 2, 3):
        d, q = 1 << bits, 1 << (2 * bits)
        minimum = (q - d) // 2 - 1
        for prime in (3, 5, 7):
            for _ in range(24):
                a0, b3 = [0] * q, [0] * q
                a0[rng.randrange(q)] = rng.randrange(1, prime)
                b3[rng.randrange(q)] = rng.randrange(1, prime)
                a3 = transform(tuple(b3), bits, prime)
                a1 = tuple(rng.randrange(prime) for _ in range(q))
                a2 = tuple((1 - a - b - c) % prime for a, b, c in zip(a0, a1, a3))
                row = (tuple(a0), a1, a2, a3)
                col = tuple(transform(values, bits, prime) for values in row)
                assert col[3] == tuple(b3)
                assert all(col[0]) and all(row[3])
                assert all(sum(p) % prime == 1 for p in zip(*row))
                assert all(sum(p) % prime == 1 for p in zip(*col))
                low = [next(i for i in range(4) if row[i][j]) for j in range(q)]
                high = [next(i for i in range(3, -1, -1) if col[i][j]) for j in range(q)]
                assert low.count(0) == 1 and high.count(3) == 1
                assert sum(x <= 2 for x in low) - 1 >= minimum
                assert sum(x >= 1 for x in high) - 1 >= minimum
                for x in (Fraction(1), Fraction(3, 2), Fraction(3), Fraction(8)):
                    a = sum(x ** (-i) for i in low)
                    b = sum(x ** i for i in high)
                    assert a >= 1 + minimum / x ** 2
                    assert b >= x ** 3 + minimum * x
                    assert a * b >= x ** 3 * (1 + minimum / x ** 2) ** 2
                    assert sum(x ** 0 for _ in col[0]) * sum(x ** 3 for _ in row[3]) >= q ** 2
                    pressure_comparisons += 1
                endpoint_cases += 1
                middle = tuple(rng.randrange(prime) for _ in range(q))
                a3q = tuple((1 - a - b - c - e) % prime
                            for a, b, c, e in zip(a0, a1, middle, a3))
                rowq = (tuple(a0), a1, middle, a3q, a3)
                colq = tuple(transform(values, bits, prime) for values in rowq)
                assert colq[4] == tuple(b3)
                lowq = [next(i for i in range(5) if rowq[i][j]) for j in range(q)]
                highq = [next(i for i in range(4, -1, -1) if colq[i][j]) for j in range(q)]
                for x in (Fraction(1), Fraction(3, 2), Fraction(3), Fraction(8)):
                    a = sum(x ** (-i) for i in lowq)
                    b = sum(x ** i for i in highq)
                    bound = 1 + minimum / x ** 3 + (q - minimum - 1) / x ** 4
                    assert a >= bound and b >= x ** 4 * bound
                    assert a * b >= x ** 4 * bound ** 2
                    quartic_comparisons += 1
                quartic_cases += 1
    for bits in range(3, 33):
        d = 1 << bits
        minimum = (d * d - d) // 2 - 1
        assert 64 * minimum >= 27 * d * d
        assert 64 * minimum - 27 * d * d == (d - 8) * (5 * d + 8)
    assert 3 ** 7 > 2 ** 11
    assert 3 ** 21 > 2 ** 32
    assert 2 * 5 ** 2 > 7 ** 2
    m = 27
    for x in (Fraction(1), Fraction(3, 2), Fraction(3), Fraction(8), Fraction(32)):
        lower = x ** 3 * (1 + m / x ** 2) ** 2
        assert lower ** 2 >= Fraction(256 * m ** 3, 27)
    assert Fraction(3) ** 3 * (1 + Fraction(27, 9)) ** 2 == 432
    quartic_certificates = []
    for d, threshold, left, right in ((8, 20, 2, 4), (16, 53, 4, 5)):
        m = (d * d - d) // 2 - 1
        n = d * d - m - 1
        power = [n, m, -threshold, 0, 1]
        shifted = [sum(Fraction(power[i] * comb(i, j) * left ** (i - j)
                                * (right - left) ** j) for i in range(j, 5))
                   for j in range(5)]
        bernstein = [sum(shifted[j] * Fraction(comb(k, j), comb(4, j))
                         for j in range(k + 1)) for k in range(5)]
        assert min(bernstein) > 0
        # Left/right derivative signs for x^2+m/x+n/x^2.
        assert max(0, 2 * left ** 3 - m) * left < 2 * n
        assert 2 * right ** 3 - m > 0
        assert (2 * right ** 3 - m) * right > 2 * n
        for x in (left, right):
            assert Fraction(x * x) + Fraction(m, x) + Fraction(n, x * x) > threshold
        assert threshold ** 7 > d ** 10
        quartic_certificates.append({"physical_dimension": d, "M": m, "remaining": n,
                                     "sqrt_pressure_floor": threshold,
                                     "compact_interval": [left, right],
                                     "positive_bernstein_coefficients":
                                         [[b.numerator, b.denominator] for b in bernstein]})
    relaxed_at_32 = Fraction(6394, 49)
    assert relaxed_at_32 == 7 ** 2 + Fraction(495, 7) + Fraction(528, 49)
    assert relaxed_at_32 ** 17 < 32 ** 24
    assert 2 * 12 ** 2 < 17 ** 2
    return {"status": "passed", "random_seed": 20261101,
            "exact_endpoint_cases": endpoint_cases,
            "exact_pressure_comparisons": pressure_comparisons,
            "exact_quartic_endpoint_cases": quartic_cases,
            "exact_quartic_pressure_comparisons": quartic_comparisons,
            "quartic_certificates": quartic_certificates,
            "dimension_32_relaxation": "At x=7, necessary sqrt-pressure is 6394/49 with rho<12/17<1/sqrt(2); no path supplied",
            "dimension_inequality_powers": [3, 32],
            "dimension_eight_pressure_floor": 432,
            "uniform_rho_floor": "strictly greater than 5/7 > 1/sqrt(2)",
            "scope": "Outer-rank-one cubic pressure at all d>=8, and quartic pressure at d=8,16; not actual-rank exclusions",
            "factoring_score": "1/2; unchanged"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = check()
    if args.output:
        with open(args.output, "w") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
