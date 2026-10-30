#!/usr/bin/env python3
"""Exact paired-chain roots of a rank-two modified transpose over F2."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from check_nilpotent_flag_crossings import (
    apply_columns, compose_columns, det_poly, inverse, nilpotent_family,
    phi, rank, reverse_at_degree, smith_series, vector_mask,
)


def power(columns, exponent):
    result = [1 << i for i in range(len(columns))]
    while exponent:
        if exponent & 1:
            result = compose_columns(columns, result)
        columns = compose_columns(columns, columns)
        exponent >>= 1
    return result


def row_rank(columns, q):
    return rank([vector_mask(row) for row in phi(columns, q)])


def construct(q):
    assert q >= 2 and not q & (q - 1)
    n, half = q * q, q * q // 2
    identity = sum(1 << (q * i + i) for i in range(q))
    # N is block diagonal with two-by-two all-one blocks.
    V = []
    for i in range(q):
        for j in range(q):
            correction = sum(
                (int(x // 2 == j // 2) & int(i // 2 == y // 2)) << (q * x + y)
                for x in range(q) for y in range(q)
            )
            V.append((1 << (q * j + i)) ^ correction)
    T = [v ^ (1 << i) for i, v in enumerate(V)]
    assert not any(compose_columns(T, T))
    assert rank(T) == half
    assert row_rank(V, q) == 2
    # An explicit trace-zero preimage of I under T.
    b0 = sum(1 << (q * (2 * i) + 2 * i + 1) for i in range(q // 2))
    assert apply_columns(T, b0) == identity
    assert (b0 & identity).bit_count() % 2 == 0
    image_basis, preimages = [identity], [b0]
    # Standard basis for the trace-zero matrix hyperplane.
    candidates = [(1 << j) ^ (1 if (1 << j) & identity else 0)
                  for j in range(1, n)]
    for b in candidates:
        a = apply_columns(T, b)
        if rank(image_basis + [a]) > len(image_basis):
            image_basis.append(a)
            preimages.append(b)
    assert len(image_basis) == half - 1
    assert rank(image_basis + [T[0]]) == half
    image_basis.append(T[0])
    preimages.append(1)  # E_00 has trace one.
    basis = image_basis + preimages
    assert rank(basis) == n
    inverse_rows = inverse([[basis[j] >> i & 1 for j in range(n)] for i in range(n)])
    inverse_rows = [vector_mask(row) for row in inverse_rows]
    S = []
    for j in range(n):
        coordinates = sum(((row >> j) & 1) << i for i, row in enumerate(inverse_rows))
        S.append(apply_columns(basis, coordinates >> 1))
    U = [v ^ (1 << i) for i, v in enumerate(S)]
    assert rank(S) == n - 1
    assert not any(power(S, n))
    assert power(S, half) == T
    assert power(U, half) == V
    assert apply_columns(S, identity) == 0
    assert all((v & identity).bit_count() % 2 == 0 for v in S)
    expected_last = [identity if (identity >> j) & 1 else 0 for j in range(n)]
    assert power(S, n - 1) == expected_last
    powers, current = [], U
    for _ in range((n - 1).bit_length()):
        powers.append(row_rank(current, q))
        current = compose_columns(current, current)
    result = dict(q=q, nilpotent_index=n, column_orders=[0] + [n] * (n - 1),
                  dyadic_power_row_ranks=powers, row_endpoint_at_zero="identity",
                  crossing_endpoint_at_one=True, terminal_row_rank=2)
    if q <= 4:
        result.update(S_columns=S, U_columns=U, paired_basis=basis)
    if q == 2:
        K = nilpotent_family(S, q, n)
        result["canonical_q2_row_orders"] = smith_series(reverse_at_degree(K, n - 1), 64)
        result["canonical_q2_row_determinant_mask"] = det_poly(K)
    return result


def small_control():
    # A separate rank-two root, with unimodular row polynomial.
    U = [3, 13, 2, 10]
    S = [v ^ (1 << i) for i, v in enumerate(U)]
    V = [14, 11, 13, 7]
    assert compose_columns(U, U) == V
    assert row_rank(U, 2) == 2
    K = nilpotent_family(S, 2, 4)
    assert det_poly(K) == 1
    mu = smith_series(reverse_at_degree(K, 3), 64)
    assert mu == [0, 0, 1, 11]
    return dict(q=2, U_columns=U, squared_columns=V,
                row_determinant_mask=1, row_orders=mu, column_orders=[0, 4, 4, 4])


def feistel_controls():
    f = [1, 1, 2, 3]
    q, n = len(f), len(f) ** 2
    U = [1 << (q * j + ((i + f[j]) % q)) for i in range(q) for j in range(q)]
    state, seen = 0, set()
    while state not in seen:
        seen.add(state)
        state = U[state].bit_length() - 1
    assert state == 0 and len(seen) == n
    r = len(set(f))
    image = {((y + f[x]) % q, (x - f[y]) % q)
             for x in range(q) for y in range(q)}
    r2 = len(image)
    assert row_rank(U, q) == r
    assert row_rank(compose_columns(U, U), q) == r2
    assert r * r2 >= n
    assert 17 * 17 > 288  # Exact strict comparison at q=8.
    return dict(function=f, full_cycle_length=n, row_rank=r, squared_row_rank=r2,
                exact_q8_threshold_comparison=[289, 288])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", default="2,4,8,16")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = dict(status="exact construction; no improved factoring exponent",
                  paired_chain_controls=[construct(int(q)) for q in args.sizes.split(",")],
                  rank_two_q2_control=small_control(), feistel_control=feistel_controls())
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
