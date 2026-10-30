#!/usr/bin/env python3
"""Exact structural checks for the all-energy weighted Weyl filtration."""

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
from time import monotonic


def labels_at(h, variables):
    if variables == 0:
        return [()] if h == 0 else []
    w = 1 << (variables - 1)
    return [low + (a,) for a in range(h // w + 1)
            for low in labels_at(h - a*w, variables-1)]


def binomial_odd(n, k):
    return 0 <= k <= n and k & ~n == 0


def binomial_polynomial(n):
    value, subset = 0, n
    while True:
        value |= 1 << subset
        if subset == 0:
            return value
        subset = (subset - 1) & n


def polynomial_multiply(a, b):
    value = 0
    while b:
        bit = b & -b
        value ^= a << (bit.bit_length()-1)
        b ^= bit
    return value


def entry(alpha, yy, beta):
    value = 1
    for j, (a, y, b) in enumerate(zip(alpha, yy, beta)):
        local = 0
        for k in range(min(a, y) + 1):
            if binomial_odd(y, k) and binomial_odd(b, a-k):
                local ^= binomial_polynomial(k << j)
        value = polynomial_multiply(value, local)
    return value


def binary_basis(vectors):
    basis = {}
    for value in vectors:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return list(basis.values())


def matrix_vector(rows, vector):
    return sum(((row & vector).bit_count() & 1) << i
               for i, row in enumerate(rows))


def finite_filtration(h):
    labels = labels_at(h, h.bit_length())
    blocks = defaultdict(list)
    for alpha in labels:
        for yy in labels:
            blocks[tuple(a-y for a, y in zip(alpha, yy))].append((alpha, yy))
    dimensions = []
    for pairs in blocks.values():
        n = len(pairs)
        coefficients = [[0] * n for _ in range(h+1)]
        for a, (alpha, yy) in enumerate(pairs):
            for b, (beta, _) in enumerate(pairs):
                polynomial = entry(alpha, yy, beta)
                assert polynomial.bit_count() & 1 == (a == b)
                while polynomial:
                    bit = polynomial & -polynomial
                    coefficients[bit.bit_length()-1][a] |= 1 << b
                    polynomial ^= bit
        nilpotent_generators = [[row ^ (1 << i) for i, row in enumerate(coefficients[0])]]
        nilpotent_generators.extend(A for A in coefficients[1:] if any(A))
        space = [1 << i for i in range(n)]
        profile = [n]
        while space:
            following = binary_basis(matrix_vector(A, v)
                                     for A in nilpotent_generators for v in space)
            if len(following) >= len(space):
                return False, dimensions, len(blocks)
            profile.append(len(following))
            space = following
        dimensions.append(profile)
    return True, dimensions, len(blocks)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    args = parser.parse_args()
    started = monotonic()
    quotient_checks = 0
    coefficient_changes = 0
    for B in range(1, 7):
        n, K = 1 << (B-1), (1 << B)-1
        for delta in range(1, n+1):
            N = K + delta
            outside = list(range(delta)) + list(range(K+1, N+1))
            assert all(N-a in outside for a in outside)
            for a in outside:
                for b in range(delta, K+1):
                    assert not binomial_odd(N-a, b-a)
                    quotient_checks += 1
        powers = [binomial_polynomial(i) for i in range(n+1)]
        for alpha in range(n+1):
            for y in range(n+1):
                for beta in range(n+1):
                    u = beta+y-alpha
                    if not 0 <= u <= n:
                        continue
                    original = converted = 0
                    for k in range(alpha+1):
                        if binomial_odd(y, k) and binomial_odd(beta, alpha-k):
                            original ^= powers[k]
                        if binomial_odd(beta, k) and binomial_odd(K-u, alpha-k):
                            converted ^= powers[alpha-k] << k
                    assert original == converted
                    coefficient_changes += 1

    simple_weight_exclusions = 0
    for R in range(1, 13):
        h, N = 1 << R, 2*R*(1 << R)+1
        mu = N-2*h
        for m in range(mu, N+1, 4):
            e = (m-mu)//2
            assert e == 0 or (e & ~m) != 0
            simple_weight_exclusions += 1

    column_witnesses = 0
    for R in range(1, 6):
        h = 1 << R
        sums = labels_at(2*h, R+1)
        exceptions = 0
        for sigma in sums:
            bits = [(1 << (j+k), j, 1 << k)
                    for j, value in enumerate(sigma)
                    for k in range(value.bit_length()) if value >> k & 1]
            if len(bits) == 1:
                exceptions += 1
                continue
            remaining, alpha = h, [0] * (R+1)
            for weight, j, value in sorted(bits, reverse=True):
                if weight <= remaining:
                    alpha[j] += value
                    remaining -= weight
            assert remaining == 0
            uu = [s-a for s, a in zip(sigma, alpha)]
            assert all(a >= 0 and u >= 0 and (a & u) == 0
                       for a, u in zip(alpha, uu))
            assert sum(a << j for j, a in enumerate(alpha)) == h
            column_witnesses += 1
        assert exceptions == R+1

    finite = []
    for h in (2, 4, 8):
        passed, profiles, block_count = finite_filtration(h)
        assert passed
        finite.append({"energy": h, "blocks": block_count,
                       "maximum_filtration_length": max(len(p)-1 for p in profiles),
                       "all_common_nilpotent_filtrations": True})
    negative, _, _ = finite_filtration(6)
    assert not negative

    report = {
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "local_quotient_invariance_coefficients": quotient_checks,
        "local_source_to_sl2_coefficient_identities": coefficient_changes,
        "simple_weight_bit_exclusions": simple_weight_exclusions,
        "column_block_no_carry_witnesses": column_witnesses,
        "finite_row_filtrations": finite,
        "non_prime_power_energy_six_negative_control": True,
        "scope": "Structural finite checks supplement the written all-energy quotient and composition-series proof; no asymptotic claim is inferred from a finite matrix table.",
        "elapsed_seconds": round(monotonic()-started, 3),
    }
    if args.output:
        Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
