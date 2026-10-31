#!/usr/bin/env python3
"""Independent full M4 Fock controls over the sixteen-element field.

The optional active-matrix comparison also audits the compiled CAR program,
after restoring its recorded common unit denominator.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'tools'))
from check_rectangular_f4_repair import build
from check_rectangular_defect_support import (
    infinity_matrices, concentrate_infinity, concentrate_zero,
)


def mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 19
    return result


TABLE = np.array([[mul(a, b) for b in range(16)] for a in range(16)], dtype=np.uint8)


def power(a, e):
    if not a and e < 0:
        raise ZeroDivisionError
    if a:
        e %= 15
    result = 1
    while e:
        if e & 1:
            result = mul(result, a)
        a = mul(a, a)
        e >>= 1
    return result


def eval_binary(p, z):
    result = 0
    for i in range(p.bit_length()-1, -1, -1):
        result = mul(result, z)^((p >> i)&1)
    return result


def eval_rat(r, z):
    return mul(eval_binary(r.n, z), power(eval_binary(r.d, z), -1))


def matrix_product(A, B):
    result = np.zeros((A.shape[0], B.shape[1]), dtype=np.uint8)
    for k in range(A.shape[1]):
        result ^= TABLE[A[:, k, None], B[None, k, :]]
    return result


def compound(A):
    n = len(A)
    @lru_cache(None)
    def actual(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        rb = rows & -rows
        i = rb.bit_length()-1
        result, remaining = 0, columns
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            j = cb.bit_length()-1
            if A[i, j]:
                result ^= mul(int(A[i, j]), actual(rows ^ rb, columns ^ cb))
        return result
    return np.array([[actual(i, j) for j in range(1 << n)] for i in range(1 << n)], dtype=np.uint8)


def rank(A):
    A = A.copy()
    pivot = 0
    for column in range(A.shape[1]):
        possible = np.flatnonzero(A[pivot:, column])
        if not len(possible):
            continue
        row = pivot+int(possible[0])
        A[[pivot, row]] = A[[row, pivot]]
        A[pivot] = TABLE[A[pivot], power(int(A[pivot, column]), -1)]
        for i in range(pivot+1, A.shape[0]):
            if A[i, column]:
                A[i] ^= TABLE[A[pivot], A[i, column]]
        pivot += 1
        if pivot == A.shape[0]:
            break
    return pivot


@lru_cache(None)
def support_data(place, M):
    procedure = concentrate_infinity if place == 'infinity' else concentrate_zero
    images, active, orders, trace = procedure(M)
    return {'place': place, 'M': M, 'active_indices': active,
            'active_orders': [orders[i] for i in active],
            'spectator_orders': [orders[i] for i in range(2*M) if i not in active],
            'active_images': {name: {'creation': [a[i].pair() for i in active],
                                    'annihilation': [b[i].pair() for i in active]}
                              for name, (a, b) in images.items()},
            'regular_shears': trace}


@lru_cache(None)
def compiled_active(place, M):
    from rectangular_defect_profiles import active_matrix
    return active_matrix(support_data(place, M), 1, 2)


def original_matrix(z, place):
    M, beta = 4, 6
    assert mul(beta, beta)^beta^1 == 0
    seed = build(tuple(range(M)), M, 1, 2, check_types=False)
    w = power(z, -1) if place == 'infinity' else z
    denominator = power(1 ^ w, seed['denominator_degree'])
    bulk, top = 1 << (M-2), 1 << (M-1)
    mx = [((a % bulk) << 1)|(1, 0, top, top|1)[a // bulk] for a in range(seed['q_x'])]
    mu = [((a % bulk) << 1)|(1, 0, top)[a // bulk] for a in range(seed['q_u'])]
    K = np.zeros((1 << (2*M), 1 << (2*M)), dtype=np.uint8)
    if place == 'infinity':
        hx, hy = [M]*M, [-M]*M
        hu = hv = [2*a-M+1 for a in range(M)]
    else:
        hx, hy = [-2*M+1+2*a for a in range(M)], [1+2*a for a in range(M)]
        hu = hv = [0]*M
    qx, qu = seed['q_x'], seed['q_u']
    for x in range(qx):
        for y in range(qu):
            for u in range(qu):
                for v in range(qx):
                    lo, hi = seed['values'][x+qx*(y+qu*(u+qu*v))]
                    coefficient = eval_binary(lo, w)^mul(beta, eval_binary(hi, w))
                    coefficient = mul(coefficient, power(denominator, -1))
                    exponent = sum(weights[a] for mask, weights in ((mx[x], hx), (mu[y], hy), (mu[u], hu), (mx[v], hv)) for a in range(M) if mask >> a & 1)
                    if place == 'infinity':
                        row, column = mx[v]|(mu[y] << M), mu[u]|(mx[x] << M)
                    else:
                        row, column = mx[x]|(mu[y] << M), mx[v]|(mu[u] << M)
                    K[row, column] = mul(coefficient, power(z, exponent))
    return K


def infinity_transform(K, z, full_core=False):
    M = 4
    _, U, Ui, V, Vi, sigma = infinity_matrices(M)
    _, active, orders, trace = concentrate_infinity(M)
    if full_core:
        trace = []  # The full-core M4 construction keeps all eight modes.
    U = np.array([[eval_rat(r, z) for r in row] for row in U], dtype=np.uint8)
    V = np.array([[eval_rat(r, z) for r in row] for row in V], dtype=np.uint8)
    for step in trace:
        row, column = step['row'], step['column']
        numerator, denominator = step['coefficient']
        coefficient = mul(eval_binary(numerator, z), power(eval_binary(denominator, z), -1))
        U[row] ^= TABLE[U[column], coefficient]
        companion = mul(coefficient, power(z, orders[column]-orders[row]))
        V[:, column] ^= TABLE[V[:, row], companion]
    return matrix_product(matrix_product(compound(U), K), compound(V))


def zero_transform(K, z):
    """Remove the two full internal Gaussians independently of CAR images."""
    M = 4
    transformed = K.copy()
    t, u = mul(z, z), 1 ^ z
    q = mul(u, u)
    for a in range(M):
        for b in range(M):
            bits = (1 << a)|(1 << (M+b))
            low = [mask for mask in range(1 << (2*M)) if not mask & bits]
            high = [mask | bits for mask in low]
            coefficient = mul(u, power(t, a & b))
            transformed[high] ^= TABLE[transformed[low], coefficient]
            if a & b == a:
                # The uv coefficient has v-index a and u-index b.
                bits = (1 << a)|(1 << (M+b))
                low = [mask for mask in range(1 << (2*M)) if not mask & bits]
                high = [mask | bits for mask in low]
                coefficient = mul(u, power(q, a))
                transformed[:, high] ^= TABLE[transformed[:, low], coefficient]
    V = np.eye(2*M, dtype=np.uint8)
    for a in range(M):
        for b in range(M):
            V[a, b] = int(a & b == a)
    return matrix_product(transformed, compound(V))


def check(z, place='infinity', compare_active=False):
    M, beta = 4, 6
    K = original_matrix(z, place)
    transformed = infinity_transform(K, z) if place == 'infinity' else zero_transform(K, z)
    data = support_data(place, M)
    active = data['active_indices']
    orders = (list(range(1-2*M, 2*M, 2)) if place == 'infinity'
              else list(range(1, 2*M, 2))*2)
    masks = [sum(((i >> j)&1) << a for j, a in enumerate(active)) for i in range(1 << len(active))]
    block = transformed[np.ix_(masks, masks)]
    if place == 'infinity':
        spectator = next(i for i in range(2*M) if i not in active)
        high = [mask|(1 << spectator) for mask in masks]
        assert not np.any(transformed[np.ix_(masks, high)])
        assert not np.any(transformed[np.ix_(high, masks)])
        assert np.array_equal(transformed[np.ix_(high, high)], TABLE[block, power(z, orders[spectator])])
    record = {'field_modulus': 'X^4+X+1', 'omega_image': beta,
              'place': place, 'z': z, 'M': M,
              'full_fock_dimension': 256, 'defect_dimension': len(block),
              'spectator_orders': data['spectator_orders'],
              'all_matrix_entries_checked': 65536,
              'specialized_full_rank': rank(K), 'specialized_defect_rank': rank(block)}
    if compare_active:
        matrix, shifts, metadata = compiled_active(place, M)
        denominator = eval_binary(int(metadata['unit_denominator_hex'], 16), z)
        inverse = power(denominator, -1)
        actual = np.zeros_like(block)
        for i, row in enumerate(matrix):
            for j, (lo, hi) in enumerate(row):
                coefficient = eval_binary(lo, z)^mul(beta, eval_binary(hi, z))
                actual[i, j] = mul(mul(coefficient, power(z, shifts[i][j])), inverse)
        mismatch = np.argwhere(actual != block)
        assert not len(mismatch), (place, z, mismatch[:8].tolist())
        record['compiled_active_entries_checked'] = len(block)**2
        record['restored_unit_denominator_degree'] = metadata['unit_denominator_degree']
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--compare-active', action='store_true',
                        help='Also compare both cuts with the compiled CAR program (requires C++17).')
    args = parser.parse_args()
    result = {'status': 'independent full Fock controls at finite field points',
              'cases': [check(z, place, args.compare_active)
                        for place in (('infinity', 'zero') if args.compare_active else ('infinity',))
                        for z in (2, 3, 4, 5)]}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)
