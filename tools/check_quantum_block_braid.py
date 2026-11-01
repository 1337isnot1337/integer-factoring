#!/usr/bin/env python3
"""Exact local identities for the uniform quantum block-braid endpoint."""

import argparse
import json
from hashlib import sha256
from pathlib import Path
from itertools import product
from paired_path_crossing import MULTIPLY as M, power, splitting_blocks, path_entry
from bounded_spin_crossing import paths
from quantum_block_braid import polynomial_entry, braid_word
from check_fixed_field_descent import Field


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def scale(a, c):
    return [[M[x][c] for x in row] for row in a]


def add(a, b):
    return [[x ^ y for x, y in zip(u, v)] for u, v in zip(a, b)]


def mm(a, b):
    r = [[0] * len(b[0]) for _ in a]
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            if x:
                for j, y in enumerate(b[k]):
                    r[i][j] ^= M[x][y]
    return r


def transpose(a):
    return list(map(list, zip(*a)))


def lift(a):
    n = len(a)
    return [
        [a[i % n][j % n] if i // n == j // n else 0 for j in range(2 * n)]
        for i in range(2 * n)
    ]


def e(n, site, q):
    a = [[0] * (1 << n) for _ in range(1 << n)]
    for w in range(1 << n):
        bits = (w >> site) & 3
        if bits in (1, 2):
            a[w][w] = power(q, -1) if bits == 1 else q
            a[w ^ (3 << site)][w] = 1
    return a


def B(n, site, q):
    return add(scale(eye(1 << n), q), e(n, site, q))


def sym(h, q):
    i = [[0] * (h + 1) for _ in range(1 << h)]
    for w in range(1 << h):
        inv = sum(j - (w & ((1 << j) - 1)).bit_count() for j in range(h) if w >> j & 1)
        i[w][w.bit_count()] = power(q, -inv)
    g = mm(transpose(i), i)
    p = scale(eye(h + 1), 0)
    for r in range(h + 1):
        assert g[r][r]
        p[r][r] = power(g[r][r], -1)
    return i, mm(mm(i, p), transpose(i))


def direct_value(four_paths, parameter, q=8):
    n = len(four_paths[0])
    frame = [
        [path_entry(p, w, q, dual=i >= 2) for w in range(1 << n)]
        for i, p in enumerate(four_paths)
    ]
    values = [
        M[frame[0][w & ((1 << n) - 1)]][frame[1][w >> n]] for w in range(1 << (2 * n))
    ]
    local = add(scale(eye(4), q), scale(e(2, 0, q), parameter))
    for site in braid_word(n, n):
        for base in range(len(values)):
            if base & (3 << site):
                continue
            indices = [base | (b << site) for b in range(4)]
            result = mm(local, [[values[i]] for i in indices])
            for i, v in zip(indices, result):
                values[i] = v[0]
    result = 0
    for word, value in enumerate(values):
        dual = M[frame[2][word & ((1 << n) - 1)]][frame[3][word >> n]]
        result ^= M[value][dual]
    normalization = power(q, -n * n // 2)
    affine = normalization ^ M[1 ^ normalization][parameter]
    return M[M[result][normalization]][affine]


models = []
checks = 0
for q in range(2, 16):
    try:
        splitting_blocks(q)
    except ValueError:
        continue
    den = q ^ power(q, -1)

    def qi(k):
        return M[power(q, k) ^ power(q, -k)][power(den, -1)]

    syms = [sym(h, q) for h in range(5)]
    for h in range(4):
        ih, _ = syms[h]
        embedding = [
            [
                ih[w % (1 << h)][r % (h + 1)] if w // (1 << h) == r // (h + 1) else 0
                for r in range(2 * (h + 1))
            ]
            for w in range(1 << (h + 1))
        ]
        up = [[0] * (h + 2) for _ in range(2 * (h + 1))]
        for r in range(h + 1):
            up[r][r] = 1
            up[h + 1 + r][r + 1] = power(q, -h + r)
        assert mm(embedding, up) == syms[h + 1][0], ("up", q, h)
        checks += 1
    for h in range(1, 4):
        f = lift(syms[h][1])
        rec = add(
            f, scale(mm(mm(f, e(h + 1, h - 1, q)), f), M[qi(h)][power(qi(h + 1), -1)])
        )
        assert rec == syms[h + 1][1], ("JW", q, h)
        checks += 1
    for h in range(1, 4):
        ih, fh = syms[h]
        ip, fp = syms[h - 1]
        cup = [[0] * (1 << (h - 1)) for _ in range(1 << (h + 1))]
        for w in range(1 << (h - 1)):
            cup[w + (1 << (h - 1))][w] = 1
            cup[w + (1 << h)][w] = q
        u = [[0] * h for _ in range(2 * (h + 1))]
        for m in range(1, h + 1):
            u[m][m - 1] = qi(m)
            u[h + m][m - 1] = M[power(q, m)][qi(h - m + 1)]
        embedding = [
            [
                ih[w % (1 << h)][r % (h + 1)] if w // (1 << h) == r // (h + 1) else 0
                for r in range(2 * (h + 1))
            ]
            for w in range(1 << (h + 1))
        ]
        assert mm(mm(lift(fh), cup), ip) == scale(mm(embedding, u), power(qi(h), -1)), (
            "down",
            q,
            h,
        )
        checks += 1
    b0, b1 = B(3, 0, q), B(3, 1, q)
    assert mm(mm(b0, b1), b0) == mm(mm(b1, b0), b1)
    e0, e1 = e(3, 0, q), e(3, 1, q)
    assert mm(e0, e0) == scale(e0, qi(2))
    assert mm(mm(e0, e1), e0) == e0
    checks += 2
    c0, c1 = [[0] * 2 for _ in range(8)], [[0] * 2 for _ in range(8)]
    for bit in range(2):
        c0[1 + (bit << 2)][bit] = 1
        c0[2 + (bit << 2)][bit] = q
        c1[bit + 2][bit] = 1
        c1[bit + 4][bit] = q
    assert mm(mm(b0, b1), c0) == scale(c1, q), ("cupslide", q)
    checks += 2
    for h in range(1, 4):
        braid = eye(1 << (h + 1))
        for site in range(h - 1, -1, -1):
            braid = mm(B(h + 1, site, q), braid)
        f = syms[h][1]
        right = [
            [f[i >> 1][j >> 1] if i % 2 == j % 2 else 0 for j in range(1 << (h + 1))]
            for i in range(1 << (h + 1))
        ]
        assert mm(braid, lift(f)) == mm(right, braid), ("naturality", q, h)
        checks += 1
    models.append(q)
entries = 0
values_checked = 0
for n in (2, 4, 6):
    ps = paths(n)
    choices = {
        (0, 0, 0, 0),
        (0, len(ps) - 1, len(ps) - 1, 0),
        (0, len(ps) - 1, 0, len(ps) - 1),
    }
    if len(ps) >= 5:
        choices |= {(1, 3, 2, 4), (4, 3, 1, 2)}
    for labels in choices:
        four = tuple(ps[i] for i in labels)
        polynomial = polynomial_entry(four)
        total = 0
        for value in polynomial:
            total ^= value
        assert polynomial[0] == int(labels[0] == labels[2] and labels[1] == labels[3])
        assert total == int(labels[0] == labels[3] and labels[1] == labels[2])
        assert len(polynomial) == n * n + 2
        for parameter in (2, 8, 11):
            value = 0
            for coefficient in polynomial[::-1]:
                value = M[value][parameter] ^ coefficient
            assert value == direct_value(four, parameter)
            values_checked += 1
        entries += 1
    assert len(braid_word(n, n)) == n * n
ps = paths(4)
tensor = {
    ids: polynomial_entry(tuple(ps[i] for i in ids))
    for ids in product(range(2), repeat=4)
}
rank = Field(0b10011).rank
coefficient_rows = [
    [tensor[ids][j] for ids in product(range(2), repeat=4)] for j in range(18)
]
assert rank(coefficient_rows) == 5
profile = []
for j in range(18):
    row = [
        [tensor[p, q, r, s][j] for r, s in product(range(2), repeat=2)]
        for p, q in product(range(2), repeat=2)
    ]
    col = [
        [tensor[p, q, r, s][j] for q, s in product(range(2), repeat=2)]
        for p, r in product(range(2), repeat=2)
    ]
    profile.append([rank(row), rank(col)])
root = Path(__file__).resolve().parents[1]
files = (
    "tools/quantum_block_braid.py",
    "tools/check_quantum_block_braid.py",
    "tools/paired_path_crossing.py",
    "tools/bounded_spin_crossing.py",
    "tools/bounded_spin_support.py",
    "tools/check_fourier_multiplication_gates.py",
    "tools/check_fixed_field_descent.py",
)
result = {
    "status": "pass",
    "models": models,
    "local_identities": checks,
    "complete_polynomial_endpoint_entries": entries,
    "independent_direct_gate_value_controls": values_checked,
    "complete_four_site_profile": {
        "labels": 2,
        "span_dimension": 5,
        "rank_pairs": profile,
    },
    "scope": "Constant local identities support the all-length endpoint proof; no favorable two-cut rank theorem.",
    "source_sha256": {f: sha256((root / f).read_bytes()).hexdigest() for f in files},
}
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
if args.output:
    args.output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
