#!/usr/bin/env python3
"""Exact all-characteristic Jacobian certificates for subspace compression.

Standard library only. The chart/properness proof is in the companion note.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from random import Random
from time import monotonic


A_TYPES = {
    1: [[[0], [1]], [[1], [1]], [[1], [1]]],
    2: [[[0, 0], [1, 1], [0, 1]],
        [[0, 0], [1, 0], [1, 1]],
        [[0, 1], [0, 0], [0, 1]]],
}
D_TYPES = {
    1: [[[1, 0]], [[1, 1]], [[1, 0]]],
    2: [[[0, 1, 0], [1, 0, 0]],
        [[0, 1, 0], [0, 0, 1]],
        [[0, 0, 1], [1, 1, 1]]],
}


def zero(rows, columns):
    return [[0] * columns for _ in range(rows)]


def identity(d):
    return [[int(i == j) for j in range(d)] for i in range(d)]


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, column)) for column in zip(*b)]
            for row in a]


def add(a, b, factor=1):
    return [[x + factor * y for x, y in zip(arow, brow)]
            for arow, brow in zip(a, b)]


def determinant(matrix):
    """Fraction-free Bareiss elimination, including exact division checks."""
    a = [row[:] for row in matrix]
    n = len(a)
    assert all(len(row) == n for row in a)
    if not n:
        return 1
    sign, previous = 1, 1
    for col in range(n - 1):
        pivot_row = next((i for i in range(col, n) if a[i][col]), None)
        if pivot_row is None:
            return 0
        if pivot_row != col:
            a[col], a[pivot_row] = a[pivot_row], a[col]
            sign = -sign
        pivot = a[col][col]
        for i in range(col + 1, n):
            for j in range(col + 1, n):
                numerator = pivot * a[i][j] - a[i][col] * a[col][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
            a[i][col] = 0
        previous = pivot
    return sign * a[-1][-1]


def integer_inverse(matrix):
    n = len(matrix)
    rows = [[Fraction(x) for x in row + unit]
            for row, unit in zip(matrix, identity(n))]
    for col in range(n):
        pivot = next(i for i in range(col, n) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        scale = rows[col][col]
        rows[col] = [x / scale for x in rows[col]]
        for i in range(n):
            if i != col:
                scale = rows[i][col]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[col])]
    answer = []
    for row in rows:
        assert all(x.denominator == 1 for x in row[n:])
        answer.append([int(x) for x in row[n:]])
    assert multiply(matrix, answer) == identity(n)
    assert multiply(answer, matrix) == identity(n)
    return answer


def rank_mod(matrix, prime):
    rows = [[x % prime for x in row] for row in matrix]
    if not rows:
        return 0
    answer = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(answer, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[answer], rows[pivot] = rows[pivot], rows[answer]
        scale = pow(rows[answer][col], -1, prime)
        rows[answer] = [x * scale % prime for x in rows[answer]]
        for i in range(answer + 1, len(rows)):
            scale = rows[i][col]
            if scale:
                rows[i] = [(x - scale * y) % prime
                           for x, y in zip(rows[i], rows[answer])]
        answer += 1
        if answer == len(rows):
            break
    return answer


def jacobian(As, Ds):
    assert len(As) == len(Ds)
    ell, s = len(As[0]), len(As[0][0])
    q, p = len(Ds[0]), len(Ds[0][0])
    rows = []
    for A, D in zip(As, Ds):
        assert len(A) == ell and all(len(row) == s for row in A)
        assert len(D) == q and all(len(row) == p for row in D)
        for u in range(q):
            for v in range(s):
                row = [0] * (p * s + q * ell)
                for h in range(p):
                    row[h * s + v] = D[u][h]
                for h in range(ell):
                    row[p * s + u * ell + h] = -A[h][v]
                rows.append(row)
    return rows


def basic_records():
    records = []
    for a_type, d_type in ((1, 1), (2, 1), (1, 2), (2, 2)):
        J = jacobian(A_TYPES[a_type], D_TYPES[d_type])
        order = 3 * a_type * d_type
        minor = [row[:order] for row in J]
        det = determinant(minor)
        assert det == (-1 if a_type == d_type == 1 else 1)
        inverse = integer_inverse(minor)
        records.append({"A_type": a_type, "D_type": d_type,
                        "rows": len(J), "columns": len(J[0]),
                        "selected_columns": list(range(order)),
                        "determinant": det, "minor": minor,
                        "integer_inverse": inverse})
    return records


def split_system(d, s, q, a, b):
    assert a >= 1 and b >= 1 and s + a * q <= d and b * s + q <= d
    p, ell, m = d - s, d - q, a + b
    As, Ds = [zero(ell, s) for _ in range(m)], [zero(q, p) for _ in range(m)]
    blocks = []
    for index in range(a):
        for u in range(q):
            Ds[index][u][index * q + u] = 1
        rows = list(range(index * q * s, (index + 1) * q * s))
        cols = [(index * q + u) * s + v for u in range(q) for v in range(s)]
        blocks.append((rows, cols, identity(q * s)))
    for index in range(b):
        for v in range(s):
            As[a + index][index * s + v][v] = 1
        rows = list(range((a + index) * q * s, (a + index + 1) * q * s))
        cols = [p * s + u * ell + index * s + v for u in range(q) for v in range(s)]
        blocks.append((rows, cols, [[-x for x in row] for row in identity(q * s)]))
    return As, Ds, blocks


def typed_blocks(total):
    parts = [2] * (total // 2) + ([1] if total % 2 else [])
    result, small_offset, large_offset = [], 0, 0
    for part in parts:
        result.append((part, small_offset, large_offset))
        small_offset += part
        large_offset += part + 1
    return result, large_offset


def three_system(d, s, q, a=0, b=0):
    assert a >= 0 and b >= 0
    assert s + (a + 1) * q + (q + 1) // 2 <= d
    assert (b + 1) * s + q + (s + 1) // 2 <= d
    p, ell, m = d - s, d - q, 3 + a + b
    As, Ds = [zero(ell, s) for _ in range(m)], [zero(q, p) for _ in range(m)]
    a_parts, ell0 = typed_blocks(s)
    d_parts, p0 = typed_blocks(q)
    for part, small, large in a_parts:
        for index in range(3):
            for u in range(part + 1):
                for v in range(part):
                    As[index][large + u][small + v] = A_TYPES[part][index][u][v]
    for part, small, large in d_parts:
        for index in range(3):
            for u in range(part):
                for v in range(part + 1):
                    Ds[index][small + u][large + v] = D_TYPES[part][index][u][v]
    blocks = []
    for ap, ai, ao in a_parts:
        for dp, do, di in d_parts:
            base = jacobian(A_TYPES[ap], D_TYPES[dp])
            size = 3 * ap * dp
            rows = [index * q * s + (do + u) * s + ai + v
                    for index in range(3) for u in range(dp) for v in range(ap)]
            cols = []
            for col in range(size):
                if col < (dp + 1) * ap:
                    u, v = divmod(col, ap)
                    cols.append((di + u) * s + ai + v)
                else:
                    u, v = divmod(col - (dp + 1) * ap, ap + 1)
                    cols.append(p * s + (do + u) * ell + ao + v)
            blocks.append((rows, cols, [row[:size] for row in base]))
    for extra in range(a):
        index = 3 + extra
        for u in range(q):
            Ds[index][u][p0 + extra * q + u] = 1
        rows = list(range(index * q * s, (index + 1) * q * s))
        cols = [(p0 + extra * q + u) * s + v for u in range(q) for v in range(s)]
        blocks.append((rows, cols, identity(q * s)))
    for extra in range(b):
        index = 3 + a + extra
        for v in range(s):
            As[index][ell0 + extra * s + v][v] = 1
        rows = list(range(index * q * s, (index + 1) * q * s))
        cols = [p * s + u * ell + ell0 + extra * s + v
                for u in range(q) for v in range(s)]
        blocks.append((rows, cols, [[-x for x in row] for row in identity(q * s)]))
    return As, Ds, blocks


def verify_system(As, Ds, blocks, counters):
    J = jacobian(As, Ds)
    rows = [row for block_rows, _, _ in blocks for row in block_rows]
    cols = [col for _, block_cols, _ in blocks for col in block_cols]
    assert len(rows) == len(set(rows)) == len(J)
    assert len(cols) == len(set(cols)) == len(J)
    assert sorted(rows) == list(range(len(J)))
    cursor = 0
    for block_rows, block_cols, expected in blocks:
        assert [[J[row][col] for col in block_cols] for row in block_rows] == expected
        outside = cols[:cursor] + cols[cursor + len(block_cols):]
        assert all(J[row][col] == 0 for row in block_rows for col in outside)
        assert abs(determinant(expected)) == 1
        cursor += len(block_cols)
        counters["unimodular_blocks"] += 1
    selected = [[J[row][col] for col in cols] for row in rows]
    if len(J) <= 24:
        assert abs(determinant(selected)) == 1
        counters["whole_integer_determinants"] += 1
    for prime in (2, 3, 5):
        assert rank_mod(selected, prime) == len(J)
        counters["finite_field_jacobian_ranks"] += 1
    counters["assembled_systems"] += 1


def compression_bound(d, m):
    best, witness = 0, None
    for a in range(1, m):
        b = m - a
        for s in range(1, d):
            q = min((d - s) // a, d - b * s)
            if q > 0 and s * q > best:
                best, witness = s * q, {"family": "split", "s": s, "q": q, "a": a, "b": b}
    for a in range(max(0, m - 2)):
        b = m - 3 - a
        for s in range(1, d):
            q = min(2 * (d - s) // (2 * a + 3), d - (b + 1) * s - (s + 1) // 2)
            if q > 0 and s * q > best:
                best, witness = s * q, {"family": "three", "s": s, "q": q, "a": a, "b": b}
    return best, witness


@lru_cache(None)
def geometric_complement(d, m):
    """Geometric complement bound only; exact small-count bounds are separate.

    A and B are twice the allocation coefficients. Even values use the
    split-map system; odd values use the compatible three-map system.
    Symmetry permits A<=B. For m=0,1 this helper returns zero.
    """
    assert d >= 1 and m >= 0
    if m < 2 or m >= 2 * d - 1:
        return 0
    best = 0
    for A in range(2, m + 1):
        B = 2 * m - A
        for s in range(1, (2*d-2)//B+1):
            q = min(2 * (d - s) // A, (2 * d - B * s) // 2)
            if q > 0:
                best = max(best, s * q)
    return best


def brute_bound(d, m):
    best = 0
    for s in range(1, d):
        for q in range(1, d):
            split = any(s + a * q <= d and (m - a) * s + q <= d
                        for a in range(1, m))
            three = any(s + (a + 1) * q + (q + 1) // 2 <= d
                        and (m - 2 - a) * s + q + (s + 1) // 2 <= d
                        for a in range(max(0, m - 2)))
            if split or three:
                best = max(best, s * q)
    return best


def random_matrix(rows, columns, rng, modulus=None):
    return [[rng.randrange(modulus) if modulus else rng.randrange(-2, 3)
             for _ in range(columns)] for _ in range(rows)]


def chart_controls(rng, counters):
    for d in range(2, 8):
        for s in range(1, d):
            q = d - s
            ell, p = d - q, d - s
            A = random_matrix(ell, s, rng)
            B = random_matrix(ell, p, rng)
            D = random_matrix(q, p, rng)
            X = random_matrix(p, s, rng)
            Y = random_matrix(q, ell, rng)
            C = add(add(multiply(Y, A), multiply(multiply(Y, B), X)), multiply(D, X), -1)
            P = [arow + brow for arow, brow in zip(A, B)] + [crow + drow for crow, drow in zip(C, D)]
            inclusion = identity(s) + X
            image = multiply(P, inclusion)
            upper = add(A, multiply(B, X))
            assert image == upper + multiply(Y, upper)
            for t in (-2, -1, 0, 1, 2):
                Xt = [[t * x for x in row] for row in X]
                Yt = [[t * x for x in row] for row in Y]
                Ct = add(add(multiply(Yt, A), multiply(multiply(Yt, B), Xt)), multiply(D, Xt), -1)
                linear = add(multiply(Y, A), multiply(D, X), -1)
                quadratic = multiply(multiply(Y, B), X)
                expected = [[t * x + t * t * y for x, y in zip(row, qrow)]
                            for row, qrow in zip(linear, quadratic)]
                assert Ct == expected
            counters["integer_chart_identities"] += 1


def complement_controls(rng, counters):
    for prime in (2, 3, 5):
        for d in range(2, 7):
            n = d * d
            E = [[int(i == v and u == j) for j in range(d) for v in range(d)]
                 for i in range(d) for u in range(d)]
            for m in range(2, min(2 * d, 8) + 1):
                bound, _ = compression_bound(d, m)
                for _ in range(3):
                    residual = [row[:] for row in E]
                    for _ in range(m):
                        P = random_matrix(d, d, rng, prime)
                        Q = random_matrix(d, d, rng, prime)
                        for i in range(d):
                            for u in range(d):
                                for j in range(d):
                                    for v in range(d):
                                        residual[i * d + u][j * d + v] -= P[i][j] * Q[u][v]
                    assert rank_mod(residual, prime) >= bound
                    counters["actual_product_complements"] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = monotonic()
    counters = Counter()
    rng = Random(202610303)
    records = basic_records()
    counters["basic_unimodular_minors"] = len(records)
    counters["basic_integer_inverses"] = len(records)

    for s in range(1, 5):
        for q in range(1, 5):
            for a, b in ((0, 0), (1, 0), (0, 1), (1, 1), (2, 0)):
                d = max(s + (a + 1) * q + (q + 1) // 2,
                        (b + 1) * s + q + (s + 1) // 2)
                for padding in (0, 1):
                    verify_system(*three_system(d + padding, s, q, a, b), counters)
    for a in range(1, 5):
        for b in range(1, 5):
            for s, q in ((1, 1), (1, 3), (2, 2), (3, 1)):
                d = max(s + a * q, b * s + q)
                verify_system(*split_system(d, s, q, a, b), counters)

    selected = []
    for d in range(2, 21):
        for m in range(2, 2 * d + 2):
            value, witness = compression_bound(d, m)
            assert value == geometric_complement(d, m)
            counters["doubled_allocation_comparisons"] += 1
            if d <= 12:
                assert value == brute_bound(d, m)
                counters["optimized_vs_brute_bounds"] += 1
            old = max(s * (d - (m - 1) * s) for s in range(d + 1))
            assert value >= old
            if m >= 2 * d - 1:
                assert value == 0
                counters["finite_closure_zero_tails"] += 1
            if m == 2 * d - 2:
                assert value >= 1
                counters["positive_tail_boundary_controls"] += 1
            if d in (5, 8, 9, 10, 16, 20) and m in (3, 4, 5, 6, 7, 8):
                selected.append({"d": d, "m": m, "bound": value,
                                 "old_flag_bound": old, "witness": witness})
            counters["old_flag_comparisons"] += 1
    assert compression_bound(5, 3)[0] == 4
    assert compression_bound(9, 3)[0] == 12
    assert compression_bound(8, 4)[0] == 6
    chart_controls(rng, counters)
    complement_controls(rng, counters)
    result = {"status": "PASS", "scope": "Universal all-characteristic subspace-compression certificates; no factoring improvement",
              "counts": dict(sorted(counters.items())), "basic_certificates": records,
              "selected_bounds": selected, "random_seed": 202610303,
              "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "seconds": monotonic() - started}
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("basic_certificates", "selected_bounds")}, indent=2))
    print(json.dumps({"selected_bounds": selected}, indent=2))


if __name__ == "__main__":
    main()
