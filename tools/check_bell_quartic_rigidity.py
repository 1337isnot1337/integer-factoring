#!/usr/bin/env python3
"""Exact coefficient-to-Hadamard-rigidity identities for quartic Bell paths."""

import argparse
from fractions import Fraction
import json
import random

from bell_crossing import matmul, rank_mod, sign, transform


def reshape(values, bits, prime):
    d = 1 << bits
    return [[sign(a & b) * values[a * d + b] % prime for b in range(d)]
            for a in range(d)]


def check():
    rng = random.Random(440993)
    identities = comparisons = 0
    for bits in (1, 2, 3, 4):
        d, q = 1 << bits, 1 << (2 * bits)
        for prime in (3, 5, 7):
            h = [[sign(a & b) % prime for b in range(d)] for a in range(d)]
            for support in (0, 1, d // 2, d, q // 3):
                shift = rng.randrange(q)
                def dual(values):
                    return transform(tuple(values[j ^ shift] for j in range(q)), bits, prime)
                a0 = [0] * q
                a0[rng.randrange(q)] = rng.randrange(1, prime)
                a1 = [0] * q
                for j in rng.sample(range(q), support):
                    a1[j] = rng.randrange(1, prime)
                a4, a3 = dual(a0), dual(a1)
                a2 = tuple((1 - a - b - c - e) % prime
                           for a, b, c, e in zip(a0, a1, a3, a4))
                assert dual(a2) == a2
                rows = (tuple(a0), tuple(a1), a2, a3, a4)
                orders = [next(i for i in range(5) if rows[i][j]) for j in range(q)]
                r1, r2 = sum(i <= 1 for i in orders), sum(i <= 2 for i in orders)
                center = reshape(a2, bits, prime)
                low = tuple((a + b) % prime for a, b in zip(a0, a1))
                low_matrix, dual_low = reshape(low, bits, prime), reshape(dual(low), bits, prime)
                residual = [[(h[i][j] - center[i][j]) % prime for j in range(d)] for i in range(d)]
                assert residual == [[(low_matrix[i][j] + dual_low[i][j]) % prime
                                     for j in range(d)] for i in range(d)]
                assert rank_mod(low_matrix, prime) == rank_mod(dual_low, prime)
                translated = reshape(tuple(low[j ^ shift] for j in range(q)), bits, prime)
                expected = matmul(matmul(h, list(map(list, zip(*translated))), prime), h, prime)
                expected = [[x * pow(d, -1, prime) % prime for x in row] for row in expected]
                assert expected == dual_low
                rank = rank_mod(residual, prime)
                changed = sum(bool(x) for row in center for x in row)
                assert rank <= 2 * r1 and changed <= r2
                minimum = (q - d) // 2 - 1
                for x in (Fraction(1), Fraction(3, 2), Fraction(3), Fraction(9)):
                    f = sum(x ** (2 - order) for order in orders)
                    pressure = f ** 2
                    assert r1 <= f / x and r2 <= f
                    assert changed ** 2 <= pressure
                    if minimum > 0:
                        assert f >= minimum / x
                        assert rank * minimum <= 2 * pressure
                    comparisons += 1
                identities += 1
    return {"status": "passed", "random_seed": 440993,
            "exact_quartic_matrix_identities": identities,
            "rational_pressure_comparisons": comparisons,
            "fields": [3, 5, 7], "physical_dimensions": [2, 4, 8, 16],
            "necessary_rigidity_bound": "R_H(floor(2Q/M)) <= floor(sqrt(Q)); M=(d^2-d)/2-1",
            "critical_rank_exponent": "2sqrt(2)-2",
            "critical_changed_entry_exponent": "sqrt(2)",
            "scope": "Necessary condition for self-reciprocal quartic Bell paths with rank-one constant coefficient; not a converse or a factoring improvement",
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
