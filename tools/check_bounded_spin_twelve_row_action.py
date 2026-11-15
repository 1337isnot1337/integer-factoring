#!/usr/bin/env python3
"""Fast row-cut monomial action for the bounded-spin coefficient algebra."""
from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import product
from math import comb
import json
from pathlib import Path

import numpy as np

from bounded_spin_crossing import paths, path_entry, path_norm
from check_bounded_spin_six_vertex import operators
from bounded_spin_batched_blocks import coefficient_tensor
from check_bounded_spin_growth import original_coefficients
from native_gf5_flint import rank

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'evidence/bounded-spin-twelve-coefficient-idempotents.json'
TEN_EVIDENCE = ROOT / 'evidence/bounded-spin-ten-algebra-model.json'


def fwht(a):
    n = len(a)
    h = 1
    while h < n:
        for i in range(0, n, 2*h):
            x = a[i:i+h].copy()
            y = a[i+h:i+2*h].copy()
            a[i:i+h] = (x+y) % 5
            a[i+h:i+2*h] = (x-y) % 5
        h *= 2
    return a


def inverse_fwht(a):
    n = len(a)
    fwht(a)
    a[:] = a*pow(n % 5, -1, 5) % 5
    return a


def prepare(n):
    ps = paths(n)
    words = [w for w in range(1 << n) if w.bit_count() == n//2]
    indices = {w:i for i,w in enumerate(words)}
    frame = np.array([[path_entry(p, w) % 5 for w in words] for p in ps], dtype=np.int64)
    groups = defaultdict(list)
    for x, y in product(words, repeat=2):
        groups[x & y, x ^ y].append((x, y))
    info = []
    for (a, d), pairs in groups.items():
        s = d.bit_count()
        t = s//2
        dpos = [i for i in range(n) if d >> i & 1]
        expand = np.zeros(1 << s, dtype=np.int64)
        for e in range(1 << s):
            value = 0
            for k, pos in enumerate(dpos):
                if e >> k & 1:
                    value |= 1 << pos
            expand[e] = value
        pc = np.array([e.bit_count() for e in range(1 << s)], dtype=np.int64)
        support = np.flatnonzero(pc == t)
        xi = np.array([indices[int(a | expand[e])] for e in support])
        yi = np.array([indices[int(a | (d ^ expand[e]))] for e in support])
        info.append((s, t, xi, yi, support, pc))
    return ps, frame, info


def apply_index(value, info, frame, q, index):
    w_total = (frame.T @ value.T @ frame) % 5
    out = np.zeros((q*q,), dtype=np.int64)
    for s, t, xi, yi, support, pc in info:
        if index > s:
            continue
        w = np.zeros(1 << s, dtype=np.int64)
        w[support] = w_total[xi, yi]
        kw = np.zeros(s+1, dtype=np.int64)
        for r in range(t+1):
            j = index-s+2*r
            if j >= 0:
                kw[2*r] = comb(2*r, j)*((-1)**j % 5) % 5
        kernel = kw[pc]
        fwht(kernel)
        fwht(w)
        kernel = (kernel*w) % 5
        inverse_fwht(kernel)
        local = kernel[support]
        out += ((frame[:, xi]*local[None, :]) @ frame[:, yi].T).reshape(q*q)
    return out % 5


def row_action_matrix(n, index):
    ps, frame, info = prepare(n)
    q = len(ps)
    matrix = np.zeros((q*q, q*q), dtype=np.int64)
    for j in range(q*q):
        value = np.zeros(q*q, dtype=np.int64)
        value[j] = 1
        matrix[:, j] = apply_index(value.reshape(q, q), info, frame, q, index)
    return matrix


def scaled_operator(coeffs, n, value, info, frame, q):
    norms = [path_norm(p) % 5 for p in paths(n)]
    inverse = np.array([pow(x*y, -1, 5) for x, y in product(norms, repeat=2)],
                       dtype=np.int64).reshape(q, q)
    scaled = (value*inverse) % 5
    out = np.zeros((q, q), dtype=np.int64)
    for index, coefficient in coeffs.items():
        if coefficient:
            out = (out+coefficient*apply_index(scaled, info, frame, q, index).reshape(q, q)) % 5
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=8)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    controls = []
    for n in (6, 8):
        monomials, tensor = coefficient_tensor(n)
        old = original_coefficients(monomials, tensor, n)
        dense, _ = operators(old, n)
        ps = paths(n)
        q = len(ps)
        for index in range(len(old)):
            fast = row_action_matrix(n, index)
            assert rank(fast) == rank(dense[index]), (n, index)
            controls.append({'length': n, 'index': index,
                             'rank': rank(fast), 'row_rank': rank(dense[index])})

    ten = json.loads(TEN_EVIDENCE.read_text())
    basis = ten['basis_indices']
    idempotent = ten['idempotents'][0]
    coefficients = {index: value for index, value in zip(basis, idempotent) if value}

    n = 10
    ps, frame, info = prepare(n)
    q = len(ps)
    rng = np.random.default_rng(20261114)
    value = rng.integers(0, 5, size=(q, q), dtype=np.int64)
    first = scaled_operator(coefficients, n, value, info, frame, q)
    second = scaled_operator(coefficients, n, first, info, frame, q)
    assert np.array_equal(first, second)
    assert rank(first) == 34
    controls.append({'length': n, 'idempotent': 0, 'rank': rank(first), 'idempotent_check': True})

    n = 12
    ps, frame, info = prepare(n)
    q = len(ps)
    coefficients = {9: 4}
    rng = np.random.default_rng(1701)
    columns = []
    for _ in range(args.samples):
        value = rng.integers(0, 5, size=(q, q), dtype=np.int64)
        image = scaled_operator(coefficients, n, value, info, frame, q)
        columns.append(image.reshape(q*q))
    matrix = np.array(columns, dtype=np.int64).T
    sample_rank = rank(matrix)
    first = columns[0].reshape(q, q)
    assert int(np.trace(first)) % 5 == 0
    assert rank(first) == 88
    result = {
        'status': 'pass',
        'prime': 5,
        'length': n,
        'coefficient': {'9': 4},
        'row_action_controls': controls,
        'sampled_columns': args.samples,
        'sampled_row_rank_lower_bound': sample_rank,
        'first_image_matrix_rank': rank(first),
        'first_image_trace_mod_5': int(np.trace(first)) % 5,
        'scope': ('Fast row-cut monomial action; controls reproduce the dense row ranks '
                  'at lengths six and eight and the published ten-site idempotent-0 '
                  'dimension. The twelve-site number is a randomized lower bound for '
                  'the row action of 4 C_9, not an exact rank or a column result.'),
    }
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text+'\n')
    print(text)


if __name__ == '__main__':
    main()
