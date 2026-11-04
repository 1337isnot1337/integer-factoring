#!/usr/bin/env python3
"""Exact unweighted cut-codec image for the interval-4 Gaussian seed.

The native weighted oracle exposes attaining nested minors but not the
Smith transformations.  This checker computes the exact F4[x] Smith
decompositions of the two interval-4 cut matrices, verifies the
reconstructions, and counts the joint image of the two row/column codec
bases over the tensor support.  It also compares the exact image with the
nested-minor prefix image and records the cheaper prefix proxy at
interval 8.  Finite diagnostic only: the weighted codec and any transfer
or score change remain open.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from reflection_dual_rail_gaussian import tensor  # noqa: E402
from native_rectangular_smith_f4 import RectangularOracle  # noqa: E402

MUL = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))


def trim(p):
    while p and p[-1] == 0:
        p.pop()
    return p


def pmul(a, b):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] ^= MUL[x][y]
    return trim(out)


def padd(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) ^ (b[i] if i < len(b) else 0)
    return trim(out)


def pdivmod(a, b):
    a = a[:]
    q = [0] * max(0, len(a) - len(b) + 1)
    inv = {1: 1, 2: 3, 3: 2}[b[-1]]
    while len(a) >= len(b) and any(a):
        trim(a)
        if len(a) < len(b):
            break
        shift = len(a) - len(b)
        c = MUL[a[-1]][inv]
        q[shift] ^= c
        for i, x in enumerate(b):
            if x:
                a[shift + i] ^= MUL[c][x]
        trim(a)
    return trim(q), a


def deg(p):
    return len(p) - 1


def ident(n):
    return [[[1] if i == j else [] for j in range(n)] for i in range(n)]


def smith_with_transforms(matrix):
    n = len(matrix)
    A = [[entry[:] for entry in row] for row in matrix]
    U = ident(n)
    V = ident(n)

    def row_swap(i, j):
        A[i], A[j] = A[j], A[i]
        U[i], U[j] = U[j], U[i]

    def col_swap(i, j):
        for row in A:
            row[i], row[j] = row[j], row[i]
        for row in V:
            row[i], row[j] = row[j], row[i]

    def row_add(i, k, factor):
        if not any(factor):
            return
        A[i] = [padd(A[i][j], pmul(factor, A[k][j])) for j in range(n)]
        U[i] = [padd(U[i][j], pmul(factor, U[k][j])) for j in range(n)]

    def col_add(j, k, factor):
        if not any(factor):
            return
        for i in range(n):
            A[i][j] = padd(A[i][j], pmul(A[i][k], factor))
        for i in range(n):
            V[i][j] = padd(V[i][j], pmul(V[i][k], factor))

    for k in range(n):
        guard = 0
        while True:
            guard += 1
            if guard > 100000:
                raise ArithmeticError('Smith reduction guard exceeded.')
            best = None
            for i in range(k, n):
                for j in range(k, n):
                    if A[i][j]:
                        d = deg(A[i][j])
                        if best is None or d < best[0]:
                            best = (d, i, j)
            if best is None:
                return U, V, A
            _, i, j = best
            if i != k:
                row_swap(i, k)
            if j != k:
                col_swap(j, k)
            pivot = A[k][k]
            if not pivot:
                raise ArithmeticError('A nonzero pivot vanished.')
            changed = False
            for i in range(k + 1, n):
                if A[i][k]:
                    quotient, remainder = pdivmod(A[i][k], pivot)
                    row_add(i, k, quotient)
                    if remainder:
                        row_swap(i, k)
                        changed = True
                        break
            if changed:
                continue
            for j in range(k + 1, n):
                if A[k][j]:
                    quotient, remainder = pdivmod(A[k][j], pivot)
                    col_add(j, k, quotient)
                    if remainder:
                        col_swap(j, k)
                        changed = True
                        break
            if changed:
                continue
            bad = None
            for i in range(k + 1, n):
                for j in range(k + 1, n):
                    if A[i][j] and pdivmod(A[i][j], A[k][k])[1]:
                        bad = i
                        break
                if bad is not None:
                    break
            if bad is not None:
                row_add(k, bad, [1])
                changed = True
                continue
            break
    return U, V, A


def matmul(X, Y):
    n = len(X)
    out = [[[] for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            acc = []
            for k in range(n):
                acc = padd(acc, pmul(X[i][k], Y[k][j]))
            out[i][j] = acc
    return out


def poly(bits):
    coefficients = []
    while bits:
        coefficients.append(bits & 3)
        bits >>= 2
    return trim(coefficients)


def cut_matrices(interval):
    data = tensor(interval)
    q = data['q_x']
    n = q * q
    values = data['values']
    cut0 = [[[] for _ in range(n)] for _ in range(n)]
    cut1 = [[[] for _ in range(n)] for _ in range(n)]
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    entry = poly(values[x + q * (y + q * (u + q * v))][0])
                    cut0[x * q + y][u * q + v] = entry[:]
                    cut1[x * q + u][y * q + v] = entry[:]
    return data, q, cut0, cut1


def gauge_patterns(q):
    stride = 4 * (q - 1)
    length = 2 * stride
    yield 'zero', [0.0] * length + [0.1]
    yield 'infinity_positive', [0.05] * stride + [0.0] * stride + [0.1]
    yield 'zero_positive', [0.0] * stride + [0.05] * stride + [0.1]
    yield 'alternating', [0.05 if i % 2 else -0.05 for i in range(length)] + [0.1]
    for seed in (1, 2, 3):
        generator = random.Random(seed)
        yield f'random_{seed}', [generator.uniform(-0.1, 0.1) for _ in range(length)] + [0.1]
    for seed in (4, 5):
        generator = random.Random(seed)
        yield f'wide_{seed}', [generator.uniform(-0.3, 0.3) for _ in range(length)] + [0.2]


def prefix_image(data, q, params=None):
    oracle = RectangularOracle(data['values'], q, q, 0, True)
    if params is None:
        params = [0.0] * (oracle.parameters - 1) + [0.1]
    maps = {}
    for cut in (0, 1):
        prefixes = oracle.prefixes(params, cut)
        row_index = {}
        col_index = {}
        for i in range(q * q):
            row_index.setdefault(prefixes['rows'][i], i)
            col_index.setdefault(prefixes['cols'][i], i)
        maps[cut] = (row_index, col_index)
    row_pairs = set()
    col_pairs = set()
    support = 0
    values = data['values']
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    if not values[x + q * (y + q * (u + q * v))][0]:
                        continue
                    support += 1
                    row0 = maps[0][0][x * q + y]
                    col0 = maps[0][1][u * q + v]
                    row1 = maps[1][0][x * q + u]
                    col1 = maps[1][1][y * q + v]
                    row_pairs.add((row0, row1))
                    col_pairs.add((col0, col1))
    return support, len(row_pairs), len(col_pairs)


def exact_image(data, q, cut0, cut1):
    U0, V0, S0 = smith_with_transforms(cut0)
    U1, V1, S1 = smith_with_transforms(cut1)
    n = q * q
    reconstruction0 = matmul(matmul(U0, cut0), V0)
    reconstruction1 = matmul(matmul(U1, cut1), V1)
    verified = all(
        reconstruction[i][j] == S[i][j]
        for reconstruction, S in ((reconstruction0, S0), (reconstruction1, S1))
        for i in range(n) for j in range(n)
    )
    diagonal = all(
        not S[i][j]
        for S in (S0, S1)
        for i in range(n) for j in range(n) if i != j
    )
    values = data['values']
    row_pairs = set()
    col_pairs = set()
    raw_pairs = set()
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    if not values[x + q * (y + q * (u + q * v))][0]:
                        continue
                    i0, j0 = x * q + y, u * q + v
                    i1, j1 = x * q + u, y * q + v
                    row_pairs.add((tuple(map(tuple, U0[i0])), tuple(map(tuple, U1[i1]))))
                    col_pairs.add((tuple(map(tuple, (V0[k][j0] for k in range(n)))),
                                   tuple(map(tuple, (V1[k][j1] for k in range(n))))))
                    raw_pairs.add((tuple(cut0[i0][j0]), tuple(cut1[i1][j1])))
    return {
        'verified': verified,
        'diagonal': diagonal,
        'invariant_degrees_cut0': [deg(S0[i][i]) for i in range(n)],
        'invariant_degrees_cut1': [deg(S1[i][i]) for i in range(n)],
        'transform_max_degree_cut0': max(
            deg(U0[i][j]) if U0[i][j] else -1 for i in range(n) for j in range(n)),
        'transform_max_degree_cut1': max(
            deg(U1[i][j]) if U1[i][j] else -1 for i in range(n) for j in range(n)),
        'row_codec_pairs': len(row_pairs),
        'column_codec_pairs': len(col_pairs),
        'raw_value_pairs': len(raw_pairs),
    }


def partition_identity(data, q, cut0, cut1):
    oracle = RectangularOracle(data['values'], q, q, 0, True)
    params = [0.0] * (oracle.parameters - 1) + [0.1]
    maps = {}
    for cut in (0, 1):
        prefixes = oracle.prefixes(params, cut)
        row_index = {}
        col_index = {}
        for i in range(q * q):
            row_index.setdefault(prefixes['rows'][i], i)
            col_index.setdefault(prefixes['cols'][i], i)
        maps[cut] = (row_index, col_index)
    U0, V0, _ = smith_with_transforms(cut0)
    U1, V1, _ = smith_with_transforms(cut1)
    n = q * q
    prefix_row_classes = {}
    codec_row_classes = {}
    prefix_column_classes = {}
    codec_column_classes = {}
    values = data['values']
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    if not values[x + q * (y + q * (u + q * v))][0]:
                        continue
                    key = (x, y, u, v)
                    i0, i1 = x * q + y, x * q + u
                    j0, j1 = u * q + v, y * q + v
                    prefix_row_key = (maps[0][0][i0], maps[1][0][i1])
                    codec_row_key = (tuple(map(tuple, U0[i0])), tuple(map(tuple, U1[i1])))
                    prefix_column_key = (maps[0][1][j0], maps[1][1][j1])
                    codec_column_key = (
                        tuple(map(tuple, (V0[k][j0] for k in range(n)))),
                        tuple(map(tuple, (V1[k][j1] for k in range(n)))))
                    prefix_row_classes.setdefault(prefix_row_key, set()).add(key)
                    codec_row_classes.setdefault(codec_row_key, set()).add(key)
                    prefix_column_classes.setdefault(prefix_column_key, set()).add(key)
                    codec_column_classes.setdefault(codec_column_key, set()).add(key)
    prefix_row_partition = {frozenset(items) for items in prefix_row_classes.values()}
    codec_row_partition = {frozenset(items) for items in codec_row_classes.values()}
    prefix_column_partition = {frozenset(items) for items in prefix_column_classes.values()}
    codec_column_partition = {frozenset(items) for items in codec_column_classes.values()}
    return {
        'row_identical': prefix_row_partition == codec_row_partition,
        'column_identical': prefix_column_partition == codec_column_partition,
        'identical': (prefix_row_partition == codec_row_partition
                      and prefix_column_partition == codec_column_partition),
        'row_class_sizes': sorted(len(items) for items in prefix_row_partition),
        'column_class_sizes': sorted(len(items) for items in prefix_column_partition),
        'row_common_classes': len(prefix_row_partition & codec_row_partition),
        'column_common_classes': len(prefix_column_partition & codec_column_partition),
        'total_classes': len(prefix_row_partition),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    data4, q4, cut04, cut14 = cut_matrices(4)
    exact = exact_image(data4, q4, cut04, cut14)
    partition = partition_identity(data4, q4, cut04, cut14)
    support4, prefix_rows4, prefix_cols4 = prefix_image(data4, q4)

    data8 = tensor(8)
    support8, prefix_rows8, prefix_cols8 = prefix_image(data8, data8['q_x'])
    gauge_counts = []
    for name, params in gauge_patterns(data8['q_x']):
        _, rows, cols = prefix_image(data8, data8['q_x'], params)
        gauge_counts.append({'gauge': name, 'row_pairs': rows, 'column_pairs': cols})
    gauge_stable = all(
        entry['row_pairs'] == prefix_rows8 and entry['column_pairs'] == prefix_cols8
        for entry in gauge_counts)
    distinct8 = len({entry[0] for entry in data8['values'] if entry[0]})
    n4 = q4 * q4
    n8 = data8['q_x'] ** 2

    result = {
        'status': 'pass',
        'scope': (
            'Exact unweighted F4[x] Smith transformations for the interval-4 '
            'two-cut Gaussian seed, and the nested-minor prefix image at '
            'interval 8. The interval-4 codec image is exact for this '
            'specialization; the interval-8 entry is only the exposed '
            'prefix proxy, not the Smith codec. No weighted codec, transfer, '
            'or score change is claimed.'
        ),
        'interval4': {
            'support': support4,
            'codec': exact,
            'prefix_row_pairs': prefix_rows4,
            'prefix_col_pairs': prefix_cols4,
            'codec_fraction': exact['row_codec_pairs'] / (n4 * n4),
            'prefix_fraction': prefix_rows4 / (n4 * n4),
            'partition_identity': partition,
        },
        'interval8': {
            'support': support8,
            'prefix_row_pairs': prefix_rows8,
            'prefix_col_pairs': prefix_cols8,
            'prefix_fraction': prefix_rows8 / (n8 * n8),
            'distinct_tensor_values': distinct8,
            'gauge_counts': gauge_counts,
            'gauge_stable': gauge_stable,
        },
        'source_sha256': {
            name: sha256((HERE / name).read_bytes()).hexdigest()
            for name in ('check_two_cut_codec_joint_image.py',
                         'reflection_dual_rail_gaussian.py',
                         'native_rectangular_smith_f4.py',
                         'native_rectangular_smith_f4.cpp')
        },
    }
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(json.dumps({
        'status': 'pass',
        'interval4_support': support4,
        'interval4_codec_row_pairs': exact['row_codec_pairs'],
        'interval4_prefix_row_pairs': prefix_rows4,
        'interval4_partition_identical': partition['identical'],
        'interval8_support': support8,
        'interval8_prefix_row_pairs': prefix_rows8,
        'interval8_prefix_fraction': result['interval8']['prefix_fraction'],
        'interval8_gauge_stable': gauge_stable,
    }, sort_keys=True))


if __name__ == '__main__':
    main()
