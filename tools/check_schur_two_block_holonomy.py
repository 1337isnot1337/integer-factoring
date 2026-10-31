#!/usr/bin/env python3
"""Exact two-block Schur holonomy decomposition and direct Gaussian seeds."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_schur_downset_crossing import binomial_polynomial, poly_matrix_multiply
from check_schur_gaussian_reduction import (
    evaluate, flip, matrix_multiply, polynomial_determinant,
    polynomial_multiply, scale, transpose,
)
from check_schur_offset_interval import block_entry
from schur_fermion import BinaryField


ROOT = Path(__file__).resolve().parent.parent


def direct_kernel(m, e, a=None, b=None):
    assert m > 0 and e > 0
    if a is None:
        assert m % 2 == 0
        a = m // 2
    if b is None:
        assert e % 2 == 0
        b = e // 2
    assert a > 0 and b > 0
    qa, qb = binomial_polynomial(a), binomial_polynomial(b)
    return [[0, qb, qa, 1], [qa, 0, 1 << m, qa],
            [qb, 1 << e, 0, qb << m], [1, qb, qa << m, 0]]


def seed_data(m, e, a=None, b=None):
    kernel = direct_kernel(m, e, a, b)
    local = []
    for x, y, u, v in product(range(4), repeat=4):
        value = block_entry(kernel, (flip(u), flip(x), v, y))
        assert value.bit_count() % 2 == int(x == v and y == u)
        local.append({'indices': [x, y, u, v], 'matrix': [[value]]})
    degree = max(item['matrix'][0][0].bit_length() - 1 for item in local)
    if a is None and b is None:
        assert degree == 3 * m + e
    return {'physical_alphabet': 4, 'auxiliary_dimension': 1, 'field': 2,
            'M': m, 'N': e, 'a': m // 2 if a is None else a,
            'b': e // 2 if b is None else b, 'degree': degree,
            'kernel_polynomial_bits': kernel, 'local': local}


def eigenvectors(n, l):
    core = [[1, 1, 1], [1, 0, 1 << n], [0, 1, 1 << l]]
    return [[core[i // n][j // n] if i % n == j % n else 0
             for j in range(3 * n)] for i in range(3 * n)]


def check_polynomial(n, l, m):
    indices = tuple(range(2 * n)) + tuple(range(l, l + n))
    size = len(indices)
    g = [[1 << (a & b) for b in indices] for a in indices]
    h = [[1 << (m - 1 - (a | b)) for b in indices] for a in indices]
    c = eigenvectors(n, l)
    common = binomial_polynomial(l + 2 * n - 1) << (m - l - 2 * n)
    exceptional = binomial_polynomial(n - 1) << (m - n)
    actual = poly_matrix_multiply(poly_matrix_multiply(g, h), c)
    expected = [[polynomial_multiply(c[i][j], common if j < 2 * n else exceptional)
                 for j in range(size)] for i in range(size)]
    assert actual == expected
    core = [[1, 1, 1], [1, 0, 1 << n], [0, 1, 1 << l]]
    assert polynomial_determinant(core) == 1 ^ (1 << n) ^ (1 << l)
    return size * size + 1


def inverse(matrix, field):
    n = len(matrix)
    a = [list(row) + [int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    mul = field.multiply
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        assert pivot is not None
        a[j], a[pivot] = a[pivot], a[j]
        iv = field.inverse[a[j][j]]
        a[j] = [mul(value, iv) for value in a[j]]
        for i in range(n):
            if i == j or not a[i][j]:
                continue
            scalar = a[i][j]
            a[i] = [x ^ mul(scalar, y) for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def check_gauge(n, l, m, z, field):
    indices = tuple(range(2 * n)) + tuple(range(l, l + n))
    size = len(indices)
    q = z ^ 1
    ident = [[int(i == j) for j in range(size)] for i in range(size)]
    p = [[int(b & ~a == 0) for b in indices] for a in indices]
    pt = transpose(p)
    g = [[field.power(z, a & b) for b in indices] for a in indices]
    h = [[field.power(z, m - 1 - (a | b)) for b in indices] for a in indices]
    xv = [[field.power(q, -m + b) if b & ~a == 0 else 0 for b in indices]
          for a in indices]
    vx = [[field.power(q, -1 - a) if b & ~a == 0 else 0 for b in indices]
          for a in indices]
    raw = {(0, 1): h, (1, 0): g, (0, 2): pt, (2, 0): pt,
           (0, 3): ident, (3, 0): ident, (1, 2): xv, (2, 1): vx,
           (1, 3): scale(g, field.power(z, -m), field),
           (3, 1): scale(h, field.power(z, -m), field),
           (2, 3): pt, (3, 2): pt}
    c = [[evaluate(value, z, field) for value in row] for row in eigenvectors(n, l)]
    assert field.determinant(c)
    ct, ti = transpose(c), transpose(inverse(c, field))
    gi, hi = inverse(g, field), inverse(h, field)
    qa = field.power(q, m // 2)
    zm = field.power(z, m)
    b = [[field.power(q, (l + 2 * n if i < 2 * n else n) // 2) if i == j else 0
          for j in range(size)] for i in range(size)]
    u0 = [ident, scale(gi, field.multiply(qa, zm), field), scale(pt, zm, field),
          scale(ident, zm, field)]
    v0 = [scale(ident, field.power(z, -m), field), hi, scale(pt, qa, field), ident]
    u = [matrix_multiply(ct, value, field) for value in u0]
    v = [matrix_multiply(value, ti, field) for value in v0]
    u[2] = matrix_multiply(b, u[2], field)
    v[1] = matrix_multiply(v[1], b, field)
    kernels = [direct_kernel(m, l + 2 * n)] * (2 * n) + [direct_kernel(m, n)] * n
    for (a, bindex), block in raw.items():
        actual = matrix_multiply(matrix_multiply(u[a], block, field), v[bindex], field)
        expected = [[evaluate(kernels[i][a][bindex], z, field) if i == j else 0
                     for j in range(size)] for i in range(size)]
        assert actual == expected, (n, l, m, z, a, bindex)
    return 12 * size * size


def check_determinants(m, e, a, b):
    k = direct_kernel(m, e, a, b)
    cuts = (((1, 3), (2, 0)), ((2, 0), (1, 3)),
            ((1, 0), (2, 3)), ((2, 3), (1, 0)))
    actual = [polynomial_determinant([[k[i][j] for j in cols] for i in rows])
              for rows, cols in cuts]
    qa, qb = binomial_polynomial(2 * a), binomial_polynomial(2 * b)
    expected = [(1 ^ qa) << m, (1 << e) ^ (qb << m), (1 << m) ^ qa, (1 << e) ^ qb]
    assert actual == expected
    return 4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--seed-output', type=Path)
    args = parser.parse_args()
    families = ((2, 4, 8), (2, 8, 16), (4, 8, 16), (4, 16, 32), (8, 32, 64))
    field = BinaryField(7)
    seeds = [seed_data(m, e) for m, e in ((8, 6), (16, 12), (32, 20), (32, 24), (64, 48))]
    for params in ((1, 1, 1, 1), (2, 3, 4, 5), (7, 9, 3, 4), (16, 6, 7, 3)):
        seeds.append(seed_data(*params))
    report = {
        'status': 'exact two-block holonomy and scalar Gaussian-sector construction; no improved factoring exponent',
        'families': [{'N': n, 'L': l, 'M': m, 'dimension': 3 * n,
                      'scalar_sector_lengths': [l + 2 * n, n], 'multiplicities': [2 * n, n]}
                     for n, l, m in families],
        'polynomial_identity_entries': sum(check_polynomial(n, l, m) for n, l, m in families),
        'finite_field_block_identity_entries': sum(check_gauge(n, l, m, z, field)
                                                  for n, l, m in families for z in (2, 3)),
        'direct_seed_endpoint_checks': 256 * len(seeds),
        'direct_seeds': [{'M': s['M'], 'N': s['N'], 'a': s['a'], 'b': s['b'], 'degree': s['degree']}
                         for s in seeds],
        'cross_determinant_identity_checks': sum(check_determinants(m, e, a, b)
                                                for m, e, a, b in product(range(1, 5), repeat=4)),
        'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ('tools/check_schur_two_block_holonomy.py',
                                       'tools/check_schur_downset_crossing.py',
                                       'tools/check_schur_gaussian_reduction.py',
                                       'tools/check_schur_offset_interval.py')},
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')
    if args.seed_output:
        args.seed_output.parent.mkdir(parents=True, exist_ok=True)
        args.seed_output.write_text(json.dumps(seed_data(16, 12), indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
