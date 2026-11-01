#!/usr/bin/env python3
"""Replay common functionals, compatible block counts, and small rank oracles."""

import argparse
from hashlib import sha256
from fractions import Fraction
from math import comb
import json
from pathlib import Path

from check_modular_crossing import rank_two
from check_modular_recombination import crossing
from check_pairing_word_barrier import multiply, power, R, H, C, D
from check_word_symmetry_barrier import log_bounds
from pairing_type_blocks import monomials, optimize, table_ranks


def boolean_rows(m, coefficients):
    rows = [0] * (1 << m)
    for column in range(1 << m):
        for support in coefficients:
            if not support & column:
                rows[support | column] |= 1 << column
    return rows


def boolean_checks():
    cases = 0
    for m in range(4):
        for mask in range(1, 1 << (1 << m)):
            supports = [s for s in range(1 << m) if mask >> s & 1]
            rows = boolean_rows(m, supports)
            actual = rank_two(rows)
            for support in supports:
                free = [u for u in range(1 << m) if not u & support]
                minor = [
                    sum(((rows[support | u] >> v) & 1) << j for j, v in enumerate(free))
                    for u in free
                ]
                assert rank_two(minor) == len(free)
                assert actual >= 1 << (m - support.bit_count())
                cases += 1
    return cases


def even_count(t, d, k, l):
    total = sum(
        comb(t, a) * comb(t - a, j) * 2**a * (d * (d - 1)) ** j * (d - 2) ** (t - a - j)
        for a in range(k, t + 1)
        for j in range(l, t - a + 1)
    )
    assert total % (1 << (k + l)) == 0
    return total >> (k + l)


def odd_count(t, d, a, l):
    total = comb(t, a) * sum(
        comb(t - a, j) * (d - 1) ** (t - a - j) * (d * (d - 1)) ** j
        for j in range(l, t - a + 1)
    )
    assert total % (1 << l) == 0
    return total >> l


def count_checks():
    cases = 0
    for d in (2, 4, 6, 8):
        for s in (1, 2, 4, 8, 16):
            t, bound = 2 * s, comb(2 * s, s) * (d * d - 2) ** s
            k_tail = (
                sum(comb(t, a) * 2**a * (d * d - 2) ** (t - a) for a in range(s, t + 1))
                >> s
            )
            j_tail = (
                sum(
                    comb(t, j) * (d * (d - 1)) ** j * d ** (t - j)
                    for j in range(s, t + 1)
                )
                >> s
            )
            assert even_count(t, d, 0, 0) == d ** (2 * t)
            assert even_count(t, d, s, 0) == k_tail >= bound
            assert even_count(t, d, 0, s) == j_tail >= k_tail
            assert bound * (t + 1) >= (4 * (d * d - 2)) ** s
            cases += 1
    for s in (1, 2, 4, 8, 16, 32):
        t, d = 2 * s, 3
        first = odd_count(t, d, s, 0)
        second = odd_count(t, d, 0, s)
        assert first == comb(t, s) * 8**s
        assert second * (1 << (s + 1)) >= 8**t
        assert min(first, second) * 2 * (t + 1) >= 1 << (5 * s)
        cases += 1
    return cases


def polynomial_checks():
    cases = 0
    for s in (1, 2, 4, 8, 16, 32):
        target = power(multiply(R, C), s)
        for row in (True, False):
            base = C if row else R
            alpha, beta = (R | H, H | C) if row else (C | H, H | R)
            even = power(base, 2 * s) ^ multiply(power(base, s), power(alpha, s))
            even ^= multiply(power(base, s), power(beta, s))
            odd = multiply(power(base, s), power(D, s) ^ power(beta, s))
            assert even == odd == target
            assert len(target) % 2 == 1
            cases += 1
    # Frobenius's two/three-term expansion is not valid for an arbitrary s.
    s = 3
    wrong = power(C, 2 * s) ^ multiply(power(C, s), power(R | H, s))
    wrong ^= multiply(power(C, s), power(H | C, s))
    assert wrong != power(multiply(R, C), s)
    return cases


def kron(left, right):
    q, result = len(right), []
    for row in left:
        indices = [j for j in range(row.bit_length()) if row >> j & 1]
        result.extend(sum(v << (j * q) for j in indices) for v in right)
    return result


def direct_coefficients(t, d):
    q = d * d
    diagonal = sum(1 << (i * d + i) for i in range(d))
    r = [diagonal if i // d == i % d else 0 for i in range(q)]
    swap = crossing(d)
    h = [a ^ b ^ (1 << i) for i, (a, b) in enumerate(zip(r, swap))]
    coefficients = {(0, 0, 0): [1]}
    for _ in range(t):
        following = {}
        for monomial, matrix in coefficients.items():
            for j, local in enumerate((r, h, swap)):
                target = tuple(v + (i == j) for i, v in enumerate(monomial))
                rows = kron(matrix, local)
                if target in following:
                    following[target] = [a ^ b for a, b in zip(following[target], rows)]
                else:
                    following[target] = rows
        coefficients = following
    return [coefficients[m] for m in monomials(t)]


def direct_checks():
    cases = 0
    for d, t in ((2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (4, 1), (4, 2)):
        matrices, ranks = direct_coefficients(t, d), table_ranks(t, d)
        for mask in range(1, len(ranks)):
            rows = [0] * len(matrices[0])
            for j, matrix in enumerate(matrices):
                if mask >> j & 1:
                    rows = [a ^ b for a, b in zip(rows, matrix)]
            assert rank_two(rows) == ranks[mask]
            cases += 1
    return cases


def binary_checks():
    cases = 0
    for t in range(2, 8193):
        even = t - (t % 2)
        blocks = [1 << j for j in range(1, even.bit_length()) if even >> j & 1]
        assert sum(blocks) == even
        assert all((m // 2) & (m // 2 - 1) == 0 for m in blocks)
        # Each denominator is at most 2(t+1); there are at most log2(t)+1.
        loss = sum((2 * (m + 1)).bit_length() for m in blocks)
        assert loss <= t.bit_length() * (t.bit_length() + 2)
        cases += 1
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-word-length", type=int, default=4, choices=range(1, 6))
    args = parser.parse_args()
    result = {
        "status": "all fixed pairing seeds obstructed; no factoring lower bound",
        "boolean_triangular_minors": boolean_checks(),
        "compatible_block_estimates": count_checks(),
        "common_polynomial_identities": polynomial_checks(),
        "independent_direct_rank_comparisons": direct_checks(),
        "binary_class_schedules": binary_checks(),
    }
    assert 14**5 > 2**19 and 2 * 29**2 > 40**2
    assert 3**3 < 2**5 and 2 * 3**2 > 4**2
    l2, u2 = log_bounds(2)
    l3, u3 = log_bounds(3)
    l14, u14 = log_bounds(14)
    rates = [
        (2, Fraction(3, 4), Fraction(3, 4)),
        (3, 5 * l2 / (4 * u3), 5 * u2 / (4 * l3)),
        (4, (2 + l14 / u2) / 8, (2 + u14 / l2) / 8),
    ]
    l11, _ = log_bounds(11)
    _, u5 = log_bounds(5)
    assert l11 / (2 * u5) > rates[-1][2]
    assert rates[-1][2] < Fraction(3, 4)
    intervals = []
    for d, lo, hi in rates:
        assert 2 * lo * lo > 1 and lo <= hi < 1
        scale = 10**12
        a = lo.numerator * scale // lo.denominator
        b = -((-hi.numerator * scale) // hi.denominator)
        intervals.append(
            {
                "seed": d,
                "lower": f"{a//scale}.{a%scale:012d}",
                "upper": f"{b//scale}.{b%scale:012d}",
            }
        )
    result["rho_intervals"] = intervals
    optima = []
    for t in range(1, args.max_word_length + 1):
        for d in (2, 3, 4, 5, 6):
            optimum = optimize(t, d)
            optimum.pop("seconds")
            optima.append(optimum)
    result["finite_full_type_optima"] = optima
    paths = [Path(__file__), Path(__file__).with_name("pairing_type_blocks.py")]
    result["source_sha256"] = {
        p.name: sha256(p.read_bytes()).hexdigest() for p in paths
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in result.items() if k != "finite_full_type_optima"})
    )


if __name__ == "__main__":
    main()
