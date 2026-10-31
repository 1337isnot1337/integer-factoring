#!/usr/bin/env python3
"""Exact sparse Schur kernels, Boolean down-set holonomy, and cut profiles.

This verifies a construction and compact sufficient profiles. It does not
claim an improved factoring exponent or transfer lower bounds through energy
projection.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
import random

from check_block_smith_crossings import smith_profile
from check_schur_gaussian_reduction import (
    evaluate, matrix_multiply, partition_from_frobenius,
    polynomial_multiply, transpose,
)
from schur_fermion import BinaryField, crossing_value_binary


ROOT = Path(__file__).resolve().parent.parent


def binomial_polynomial(exponent, step=1):
    value = 1
    shift = step
    while exponent:
        if exponent & 1:
            value ^= value << shift
        exponent >>= 1
        shift <<= 1
    return value


def poly_matrix_multiply(a, b):
    result = [[0] * len(b[0]) for _ in a]
    for i, row in enumerate(a):
        for k, value in enumerate(row):
            if value:
                for j, other in enumerate(b[k]):
                    result[i][j] ^= polynomial_multiply(value, other)
    return result


def downset(bits, radius):
    return tuple(a for a in range(1 << bits) if a.bit_count() <= radius)


def holonomy_matrices(indices, m):
    n = len(indices)
    p = [[int(b & ~a == 0) for b in indices] for a in indices]
    g = [[1 << (a & b) for b in indices] for a in indices]
    h = [[1 << (m - 1 - (a | b)) for b in indices] for a in indices]
    zz = [[1 << (a & ~b) for b in indices] for a in indices]
    zz2 = poly_matrix_multiply(zz, zz)
    f = [[zz2[i][j] ^ (zz2[i][j] << 1) ^ int(i == j)
          for j in range(n)] for i in range(n)]
    gh = poly_matrix_multiply(g, h)
    linf = [[gh[i][j] ^ (gh[i][j] << 1) ^ (int(i == j) << m)
             for j in range(n)] for i in range(n)]
    return p, g, h, zz, f, linf


def sparse_kernel(indices, m):
    # Polynomial bits use w, with z=w^2 and sqrt(z+1)=w+1.
    n = len(indices)
    result = [[0] * (4 * n) for _ in range(4 * n)]
    q = {a: binomial_polynomial(a, 2) for a in range(m)}
    for i, a in enumerate(indices):
        for j, b in enumerate(indices):
            gv = 1 << (2 * (a & b))
            hv = 1 << (2 * (m - 1 - (a | b)))
            blocks = {(0, 1): (hv ^ (hv << 1)) << (2 * m),
                      (1, 0): (gv ^ (gv << 1)) << (2 * m),
                      (1, 3): gv ^ (gv << 1),
                      (3, 1): hv ^ (hv << 1)}
            if a & ~b == 0:
                low = q[a] ^ (q[a] << 1)
                high = q[m - 1 - b] ^ (q[m - 1 - b] << 1)
                blocks.update({(0, 2): high, (2, 0): low,
                               (2, 3): low, (3, 2): high})
            if b & ~a == 0:
                blocks.update({(1, 2): 1 << (2 * m), (2, 1): 1 << (2 * m)})
            if a == b:
                blocks.update({(0, 3): 1, (3, 0): 1})
            for (row, column), value in blocks.items():
                result[row * n + i][column * n + j] = value
    return result


def cross_blocks(matrix, n, cut):
    first, second = ((1, 3), (2, 0)) if cut == 0 else ((1, 0), (2, 3))
    result = []
    for rows, columns in ((first, second), (second, first)):
        result.append([[matrix[a * n + i][b * n + j]
                        for b in columns for j in range(n)]
                       for a in rows for i in range(n)])
    return result


def infinity_gauge(matrix, indices, m):
    n = len(indices)
    row_heights = ([0] * n + [-1 - 2 * a for a in indices]
                   + [-1 - 2 * a for a in indices] + [0] * n)
    column_heights = ([-2 * m] * n + [-4 * m + 1 + 2 * a for a in indices]
                      + [-2 * m + 1 + 2 * a for a in indices] + [0] * n)
    result = []
    for i, row in enumerate(matrix):
        out = []
        for j, value in enumerate(row):
            replacement = 0
            while value:
                bit = value & -value
                degree = bit.bit_length() - 1
                exponent = -row_heights[i] - column_heights[j] - degree
                assert exponent >= 0
                replacement ^= 1 << exponent
                value ^= bit
            out.append(replacement)
        result.append(out)
    return result


def exact_orders(matrix):
    degree = max(value.bit_length() - 1 for row in matrix for value in row)
    return list(smith_profile(matrix, degree)[0])


def check_holonomy(indices, m):
    n = len(indices)
    p, g, h, zz, f, linf = holonomy_matrices(indices, m)
    ident = [[int(i == j) for j in range(n)] for i in range(n)]
    assert poly_matrix_multiply(p, p) == ident
    qp = [[binomial_polynomial(indices[i]) * int(i == j) for j in range(n)] for i in range(n)]
    assert poly_matrix_multiply(poly_matrix_multiply(p, qp), transpose(p)) == g
    pz = [[(1 << (a - b)) if b & ~a == 0 else 0 for b in indices] for a in indices]
    assert poly_matrix_multiply(poly_matrix_multiply(pz, qp), transpose(p)) == zz
    for i in range(n):
        for j in range(n):
            assert f[i][j].bit_length() <= m + 1
            reverse = sum(((f[i][j] >> e) & 1) << (m - e) for e in range(m + 1))
            assert reverse == linf[j][i]
    hs = exact_orders(h)
    fs = exact_orders(f)
    assert exact_orders(linf) == hs
    assert min(hs) >= 0 and max(hs) <= m - 1
    assert sum(hs) == n * (m - 1) - 2 * sum(indices)
    assert min(fs) >= 1
    assert sum(fs) <= n + 2 * sum(indices)
    return hs, fs


def check_profiles(bits, radius):
    assert 2 * radius < bits
    m = 1 << bits
    indices = downset(bits, radius)
    n = len(indices)
    hs, fs = check_holonomy(indices, m)
    matrix = sparse_kernel(indices, m)
    infinity = infinity_gauge(matrix, indices, m)
    zero_r = sorted(sum((exact_orders(a) for a in cross_blocks(matrix, n, 0)), []))
    zero_c = sorted(sum((exact_orders(a) for a in cross_blocks(matrix, n, 1)), []))
    inf_r = sorted(sum((exact_orders(a) for a in cross_blocks(infinity, n, 0)), []))
    inf_c = sorted(sum((exact_orders(a) for a in cross_blocks(infinity, n, 1)), []))
    assert zero_r == [0] * (2 * n) + [2 * m] * n + [4 * m] * n
    assert zero_c == sorted([0] * (3 * n) + [2 * v for v in hs])
    assert inf_r == [0] * (4 * n)
    assert inf_c == sorted([0] * n + [2 * m] * (2 * n) + [2 * m + 2 * v for v in fs])
    return {"bits": bits, "radius": radius, "M": m, "mode_count": n,
            "H_zero_orders_in_z": hs, "F_zero_orders_in_z": fs,
            "H_determinant_order": sum(hs),
            "F_determinant_degree": n + 2 * sum(indices),
            "one_particle_orders_in_w": {"zero_row": dict(Counter(zero_r)),
                                           "zero_column": dict(Counter(zero_c)),
                                           "gauged_infinity_row": dict(Counter(inf_r)),
                                           "gauged_infinity_column": dict(Counter(inf_c))}}


def check_schur_entries(bits, radius, r, field):
    m = 1 << bits
    indices = downset(bits, radius)
    n = len(indices)
    subsets = tuple(combinations(indices, r))
    labels = []
    for a, b in product(subsets, repeat=2):
        if sum(a) == sum(b):
            legs = tuple(sorted(m - 1 - v for v in b))
            labels.append((a, b, partition_from_frobenius(a, legs)))
    dimension = len(labels)
    assert dimension * (r * (m - 1) + 1) >= comb(n, r) ** 2
    lookup = {a: i for i, a in enumerate(indices)}
    p = [[int(b & ~a == 0) for b in indices] for a in indices]
    endpoint_checks = 0
    for x, v in product(labels, repeat=2):
        left = field.determinant([[p[lookup[a]][lookup[b]] for b in v[1]] for a in x[1]])
        right = field.determinant([[p[lookup[a]][lookup[b]] for b in x[0]] for a in v[0]])
        assert field.multiply(left, right) == int(x == v)
        endpoint_checks += 1
    polynomial = sparse_kernel(indices, m)
    w = 2
    z = field.power(w, 2)
    matrix = [[evaluate(value, w, field) for value in row] for row in polynomial]
    if dimension <= 6:
        queries = list(product(range(dimension), repeat=4))
    else:
        rng = random.Random(1000 * bits + 10 * radius + r)
        queries = [tuple(rng.randrange(dimension) for _ in range(4)) for _ in range(256)]
    for xi, yi, ui, vi in queries:
        x, y, u, v = [labels[j] for j in (xi, yi, ui, vi)]
        plus = (u[1], x[1], v[0], y[0])
        minus = (u[0], x[0], v[1], y[1])
        rows = [block * n + lookup[a] for block in range(4) for a in plus[block]]
        columns = [block * n + lookup[a] for block in range(4) for a in minus[block]]
        actual = field.determinant([[matrix[i][j] for j in columns] for i in rows])
        expected = crossing_value_binary(x[2], y[2], u[2], v[2], z, field)
        assert actual == expected, (bits, radius, r, xi, yi, ui, vi)
    return {"bits": bits, "radius": radius, "r": r, "h": m * r,
            "mode_count": n, "selected_labels": dimension,
            "full_occupation_dimension": comb(n, r) ** 2,
            "endpoint_pairing_entries": endpoint_checks,
            "complete_schur_entry_checks": len(queries)}


def negative_control():
    indices = (0, 3)
    p, g, _, _, _, _ = holonomy_matrices(indices, 4)
    q = [[binomial_polynomial(indices[i]) * int(i == j) for j in range(2)] for i in range(2)]
    product_ = poly_matrix_multiply(poly_matrix_multiply(p, q), transpose(p))
    assert product_ != g
    _, proper_g, proper_h, _, _, _ = holonomy_matrices((0, 1, 2), 4)
    off_diagonal = poly_matrix_multiply(proper_h, proper_g)[0][1]
    assert off_diagonal == 2
    return {"not_a_downset": list(indices),
            "factorization_difference_bits": product_[1][1] ^ g[1][1],
            "proper_downset": [0, 1, 2],
            "proper_downset_HG_off_diagonal_bits": off_diagonal}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    field = BinaryField(7)
    result = {
        "status": "exact sparse Schur construction and compact profiles; no improved factoring exponent",
        "profiles": [check_profiles(b, k) for b, k in ((3, 1), (4, 1), (5, 2))],
        "schur_checks": [check_schur_entries(b, k, r, field)
                         for b, k, r in ((3, 1, 1), (3, 1, 2), (3, 1, 3), (4, 2, 2))],
        "negative_control": negative_control(),
        "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ("tools/check_schur_downset_crossing.py",
                                       "tools/check_schur_gaussian_reduction.py",
                                       "tools/check_block_smith_crossings.py",
                                       "tools/schur_fermion.py")},
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
