#!/usr/bin/env python3
"""Exact checks of the integral Gauss coefficients and weighted Weyl algebra."""

import argparse
from collections import defaultdict
import hashlib
import json
from math import comb
from pathlib import Path
from time import monotonic

from check_weighted_weyl_filtration import labels_at, entry, binomial_polynomial


def choose(n, k):
    if k < 0:
        return 0
    if n >= 0:
        return comb(n, k) if k <= n else 0
    return (-1)**k * comb(k-n-1, k)


def multiply(A, B):
    result = []
    for row in A:
        value = 0
        while row:
            bit = row & -row
            value ^= B[bit.bit_length()-1]
            row ^= bit
        result.append(value)
    return result


def shifted_polynomial(value):
    answer = 0
    while value:
        bit = value & -value
        answer ^= binomial_polynomial(bit.bit_length()-1)
        value ^= bit
    return answer


def check_shell(h):
    labels = labels_at(h, h.bit_length())
    blocks = defaultdict(list)
    for alpha in labels:
        for yy in labels:
            blocks[tuple(a-y for a, y in zip(alpha, yy))].append((alpha, yy))
    H = 2*h + (1 << h.bit_length())-1
    products = 0
    additive_failures = 0
    for pairs in blocks.values():
        n = len(pairs)
        B = [[0]*n for _ in range(h+1)]
        for r, (alpha, yy) in enumerate(pairs):
            for c, (beta, _) in enumerate(pairs):
                value = shifted_polynomial(entry(alpha, yy, beta))
                while value:
                    bit = value & -value
                    B[bit.bit_length()-1][r] |= 1 << c
                    value ^= bit
        assert B[0] == [1 << i for i in range(n)]
        for i in range(h+1):
            for j in range(h+1):
                actual = multiply(B[i], B[j])
                expected = [0]*n
                for k in range(min(i, j)+1):
                    index = i+j-k
                    if index > h:
                        continue
                    coefficient = choose(index, i)*choose(index, j)*choose(H-index-1, k)
                    if coefficient & 1:
                        expected = [a ^ b for a, b in zip(expected, B[index])]
                assert actual == expected
                simple = B[i+j] if i+j <= h and choose(i+j, i) & 1 else [0]*n
                additive_failures += actual != simple
                products += 1
    if h >= 2 and h & (h-1) == 0:
        assert additive_failures == 0
    if h == 6:
        assert additive_failures > 0
    return {"energy": h, "labels": len(labels), "row_blocks": len(blocks),
            "exact_coefficient_products": products,
            "additive_law_failures": additive_failures}


def check_opposite_shell(h):
    labels = labels_at(h, h.bit_length())
    blocks = defaultdict(list)
    for alpha in labels:
        for yy in labels:
            blocks[tuple(a+y for a, y in zip(alpha, yy))].append((alpha, yy))
    products = 0
    for pairs in blocks.values():
        n = len(pairs)
        ids = {pair: i for i, pair in enumerate(pairs)}
        Q = [[0]*n for _ in range(h+1)]
        for r, (alpha, _) in enumerate(pairs):
            for c, (beta, uu) in enumerate(pairs):
                value = entry(alpha, uu, beta)
                while value:
                    bit = value & -value
                    Q[bit.bit_length()-1][r] |= 1 << c
                    value ^= bit
        K = [[Q[h-i][ids[(yy, alpha)]] for alpha, yy in pairs] for i in range(h+1)]
        I, zero = [1 << i for i in range(n)], [0]*n
        assert K[0] == I
        for i in range(h+1):
            assert multiply(Q[i], Q[i]) == (I if i == h else zero)
            for j in range(h+1):
                expected = K[i+j] if i+j <= h and choose(i+j, i) & 1 else zero
                assert multiply(K[i], K[j]) == expected
                assert multiply(Q[i], Q[j]) == multiply(Q[j], Q[i])
                products += 1
    return {"energy": h, "labels": len(labels), "sum_blocks": len(blocks),
            "exact_normalized_products": products,
            "coefficient_squares_and_commutation": True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    args = parser.parse_args()
    started = monotonic()
    gauss = 0
    ordinary_gauss = 0
    for a in range(7):
        for b in range(7):
            H0 = a+b+1
            for i in range(7):
                for j in range(7):
                    left = choose(b+j, j)*choose(a+j, i) if b+j >= i else 0
                    right = 0
                    for k in range(min(i, j)+1):
                        if i-k > min(a, b):
                            continue
                        right += (choose(H0+j-i+k-1, k)*choose(a, i-k)
                                  *choose(b+j-i, j-k))
                    assert left == right
                    gauss += 1
                    ordinary_left = choose(a, j)*choose(b+j, i)
                    ordinary_right = sum(choose(a-b+i-j, k)*choose(b, i-k)
                                         *choose(a+i-k, j-k)
                                         for k in range(min(i, j)+1))
                    assert ordinary_left == ordinary_right
                    ordinary_gauss += 1

    carry = 0
    for N in range(256):
        subsets = []
        value = N
        while True:
            subsets.append(value)
            if value == 0:
                break
            value = (value-1) & N
        for I in subsets:
            for J in subsets:
                K = I+J-N
                if K >= 0 and K & N == 0:
                    assert K == 0 and I & J == 0
                carry += 1
    shells = [check_shell(h) for h in (1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16)]
    opposite = [check_opposite_shell(h) for h in (1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16)]
    source = Path(__file__)
    dependency = source.with_name("check_weighted_weyl_filtration.py")
    result = {"source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "coefficient_helper_sha256": hashlib.sha256(dependency.read_bytes()).hexdigest(),
              "integral_gauss_coefficient_checks": gauss,
              "integral_sl2_gauss_coefficient_checks": ordinary_gauss,
              "binary_carry_checks": carry,
              "shells": shells,
              "opposite_shells": opposite,
              "scope": "The complete finite checks replay the written all-energy coefficient algebra and power-of-two additive specialization; mixed row orientations are not covered.",
              "elapsed_seconds": round(monotonic()-started, 3)}
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
