#!/usr/bin/env python3
"""Exact offset-interval Schur reduction and its polynomial Gaussian seed.

All checks concern identities and effective constructions. They do not
certify an improved rank-width or integer-factoring exponent.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path

from check_schur_downset_crossing import holonomy_matrices, poly_matrix_multiply
from check_schur_gaussian_reduction import (
    apply_axis, evaluate, flip, lift, matrix_multiply, partition_from_frobenius,
    polynomial_determinant, scale, transpose,
)
from schur_fermion import BinaryField, crossing_value_binary


ROOT = Path(__file__).resolve().parent.parent


def offset_kernel(m, n):
    assert m >= n >= 2 and m & (m - 1) == n & (n - 1) == 0
    qm = 1 | (1 << (m // 2))
    qn = 1 | (1 << (n // 2))
    zm, zn = 1 << m, 1 << n
    return [[0, qn, qm, 1], [qm, 0, zm, qm],
            [qn, zn, 0, zm * qn], [1, qn, zm * qm, 0]]


def block_entry(kernel, states):
    rows = [i for i, state in enumerate(states) if state & 1]
    columns = [i for i, state in enumerate(states) if state & 2]
    if len(rows) != len(columns):
        return 0
    return polynomial_determinant([[kernel[i][j] for j in columns] for i in rows])


def seed_data(m, n):
    kernel = offset_kernel(m, n)
    local = []
    for x, y, u, v in product(range(4), repeat=4):
        value = block_entry(kernel, (flip(u), flip(x), v, y))
        assert value.bit_count() % 2 == int(x == v and y == u)
        local.append({"indices": [x, y, u, v], "matrix": [[value]]})
    degree = max(item["matrix"][0][0].bit_length() - 1 for item in local)
    assert degree == 3 * m + n
    determinant = polynomial_determinant(kernel)
    expected = 1 ^ (1 << (m + n)) ^ (1 << (2 * m + n)) ^ (1 << (3 * m)) ^ (1 << (3 * m + n))
    assert determinant == expected
    return {"physical_alphabet": 4, "auxiliary_dimension": 1,
            "field": 2, "M": m, "N": n, "degree": degree,
            "label_encoding": "0=empty, 1=plus, 2=minus, 3=plus_minus",
            "kernel_polynomial_bits": kernel, "local": local}


def raw_kernel(m, n, z, field):
    # Direct contractions after relabelling opposite-species indices M-1-a.
    q = z ^ 1
    ident = [[int(a == b) for b in range(n)] for a in range(n)]
    p = [[int(b & ~a == 0) for b in range(n)] for a in range(n)]
    g = [[field.power(z, a & b) for b in range(n)] for a in range(n)]
    h = [[field.power(z, m - 1 - (a | b)) for b in range(n)] for a in range(n)]
    xv = [[field.power(q, -m + b) if b & ~a == 0 else 0
           for b in range(n)] for a in range(n)]
    vx = [[field.power(q, -1 - a) if b & ~a == 0 else 0
           for b in range(n)] for a in range(n)]
    pt = transpose(p)
    kernel = {(0, 1): h, (1, 0): g, (0, 2): pt, (2, 0): pt,
              (0, 3): ident, (3, 0): ident, (1, 2): xv, (2, 1): vx,
              (1, 3): scale(g, field.power(z, -m), field),
              (3, 1): scale(h, field.power(z, -m), field),
              (2, 3): pt, (3, 2): pt}
    return ident, p, g, h, kernel


def check_gauges(m, n, z, field):
    q = z ^ 1
    ident, p, g, h, raw = raw_kernel(m, n, z, field)
    # G H = z^(M-N) q^(N-1) I gives both inverses directly.
    scalar = field.multiply(field.power(z, m - n), field.power(q, n - 1))
    assert matrix_multiply(g, h, field) == scale(ident, scalar, field)
    assert matrix_multiply(h, g, field) == scale(ident, scalar, field)
    gi = scale(h, field.inverse[scalar], field)
    hi = scale(g, field.inverse[scalar], field)
    pt = transpose(p)
    u = [ident, scale(gi, field.multiply(field.power(q, m // 2), field.power(z, m)), field),
         scale(pt, field.multiply(field.power(q, n // 2), field.power(z, m)), field),
         scale(ident, field.power(z, m), field)]
    v = [scale(ident, field.power(z, -m), field), scale(hi, field.power(q, n // 2), field),
         scale(pt, field.power(q, m // 2), field), ident]
    small = offset_kernel(m, n)
    checks = 2 * n * n
    for (a, b), block in raw.items():
        actual = matrix_multiply(matrix_multiply(u[a], block, field), v[b], field)
        expected = scale(ident, evaluate(small[a][b], z, field), field)
        assert actual == expected, (m, n, z, a, b)
        checks += n * n
    return checks


def check_polynomial_holonomy(m, n):
    _, g, h, _, f, linfinity = holonomy_matrices(tuple(range(n)), m)
    ident = [[int(i == j) for j in range(n)] for i in range(n)]
    scalar = ((1 << n) - 1) << (m - n)  # z^(M-N) (1+z)^(N-1).
    expected = [[value * scalar for value in row] for row in ident]
    assert poly_matrix_multiply(g, h) == expected
    assert poly_matrix_multiply(h, g) == expected
    assert f == [[value << n for value in row] for row in ident]
    assert linfinity == [[value << (m - n) for value in row] for row in ident]
    return 4 * n * n


def check_projected_identity(m, n, r, z, field):
    # Apply all exterior maps to the full coherent tensor before projecting.
    q = z ^ 1
    ident, p, g, _, _ = raw_kernel(m, n, z, field)
    hn = [[field.power(z, n - 1 - (a | b)) for b in range(n)] for a in range(n)]
    subsets = tuple(combinations(range(n), r))
    labels = tuple(product(subsets, repeat=2))
    dimension = len(labels)
    x_map = lift(g, hn, labels, field)
    v_map = lift(transpose(p), transpose(p), labels, field)
    small = offset_kernel(m, n)
    site = {states: evaluate(block_entry(small, states), z, field)
            for states in product(range(4), repeat=4)}
    tensor = []
    for positions in product(range(dimension), repeat=4):
        value = 1
        for index in range(n):
            states = tuple(int(index in labels[k][0]) + 2 * int(index in labels[k][1])
                           for k in positions)
            value = field.multiply(value, site[states])
        tensor.append(value)
    tensor = apply_axis(tensor, x_map, 1, dimension, field)
    tensor = apply_axis(tensor, v_map, 2, dimension, field)
    lookup = {label: i for i, label in enumerate(labels)}
    physical = [(a, b, partition_from_frobenius(a, [m - 1 - c for c in b]))
                for a, b in labels if sum(a) == sum(b)]
    assert len(physical) * (r * (n - r) + 1) >= comb(n, r) ** 2
    scalar = field.multiply(field.power(z, (m - n) * r), field.power(q, -r * (n - 1)))
    checks = 0
    for x, y, u, v in product(physical, repeat=4):
        positions = [lookup[(u[1], u[0])], lookup[(x[1], x[0])],
                     lookup[(v[0], v[1])], lookup[(y[0], y[1])]]
        index = 0
        for position in positions:
            index = index * dimension + position
        actual = field.multiply(tensor[index], scalar)
        expected = crossing_value_binary(x[2], y[2], u[2], v[2], z, field)
        assert actual == expected, (m, n, r, z, x, y, u, v)
        checks += 1
    return {"M": m, "N": n, "r": r, "h": m * r,
            "full_occupation_dimension": dimension, "energy_dimension": len(physical),
            "complete_projected_schur_checks": checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seed-output", type=Path)
    args = parser.parse_args()
    field = BinaryField(7)
    pairs = ((2, 2), (4, 2), (8, 2), (8, 4), (16, 4), (32, 8), (64, 8), (256, 16))
    seeds = [seed_data(m, n) for m, n in pairs]
    result = {
        "status": "exact offset-interval construction and projection identity; no improved factoring exponent",
        "seed_checks": [{"M": seed["M"], "N": seed["N"], "degree": seed["degree"],
                         "endpoint_entries": 256} for seed in seeds],
        "polynomial_holonomy_identity_entries": sum(check_polynomial_holonomy(m, n) for m, n in pairs),
        "finite_field_one_particle_identity_entries": sum(check_gauges(m, n, z, field)
                                                           for m, n in pairs for z in (2, 3, 5)),
        "projected_schur_checks": [check_projected_identity(m, n, r, z, field)
                                    for m, n, r, z in ((4, 2, 1, 2), (8, 4, 1, 2),
                                                       (8, 4, 3, 3), (16, 4, 1, 5))],
        "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ("tools/check_schur_offset_interval.py",
                                       "tools/check_schur_downset_crossing.py",
                                       "tools/check_schur_gaussian_reduction.py",
                                       "tools/schur_fermion.py")},
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    else:
        print(text, end="")
    if args.seed_output:
        args.seed_output.parent.mkdir(parents=True, exist_ok=True)
        args.seed_output.write_text(json.dumps(seed_data(4, 2), indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
