#!/usr/bin/env python3
"""Exact residual quartic codes, direct polynomial seed, and typed endpoints.

The resulting seed is a candidate for pressure optimization.  These checks
do not certify an improved asymptotic factoring exponent.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

from check_schur_downset_crossing import sparse_kernel
from weighted_smith_oracle import polynomial_exact_divide, polynomial_multiply


ROOT = Path(__file__).resolve().parent.parent


def residual_operator():
    # Variable order x0,xT,y0,yT,u0,uT,v0,vT.
    factors = [(1 << 0) | (1 << 6), (1 << 1) | (1 << 7),
               (1 << 4) | (1 << 2), (1 << 5) | (1 << 3),
               (1 << 5) | (1 << 6), (1 << 1) | (1 << 2),
               (1 << 0) | (1 << 3) | (1 << 4) | (1 << 7)]
    support = {0}
    for monomial in factors:
        added = {mask | monomial for mask in support if not mask & monomial}
        support ^= added
    rows = [0] * 16
    for x, y, u, v in product(range(4), repeat=4):
        if (x | y << 2 | u << 4 | v << 6) in support:
            rows[4 * y + v] |= 1 << (4 * u + x)
    return rows


def apply(rows, vector):
    return sum(((row & vector).bit_count() & 1) << i for i, row in enumerate(rows))


def kron(a, b):
    return sum(1 << (4 * i + j) for i, j in product(range(4), repeat=2)
               if a >> i & 1 and b >> j & 1)


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def subspaces(dimension):
    result = {}
    for basis in combinations(range(1, 16), dimension):
        span = {0}
        for vector in basis:
            span |= {value ^ vector for value in span}
        if len(span) == 1 << dimension:
            result.setdefault(tuple(sorted(span)), basis)
    return list(result.values())


def duals(basis):
    choices = [[row for row in range(16)
                if all(((row & vector).bit_count() & 1) == int(i == j)
                       for j, vector in enumerate(basis))]
               for i in range(len(basis))]
    return product(*choices)


def restrict(rows, input_first, input_second, output_first, output_second):
    images = [apply(rows, kron(a, b)) for a in input_first for b in input_second]
    return [sum(((kron(a, b) & image).bit_count() & 1) << i
                for i, image in enumerate(images))
            for a in output_first for b in output_second]


def is_invertible_product(rows, dimension):
    """Exact realignment-rank-one test, including invertibility of both factors."""
    realigned = [0] * (dimension * dimension)
    for i, j, k, ell in product(range(dimension), repeat=4):
        if rows[dimension * i + j] >> (dimension * k + ell) & 1:
            realigned[dimension * i + k] |= 1 << (dimension * j + ell)
    second = next((row for row in realigned if row), 0)
    if not second or any(row not in (0, second) for row in realigned):
        return False
    first = sum(bool(row) << i for i, row in enumerate(realigned))
    mask = (1 << dimension) - 1
    return (rank([(first >> (dimension * i)) & mask for i in range(dimension)]) == dimension
            and rank([(second >> (dimension * i)) & mask for i in range(dimension)]) == dimension)


def classify_residual():
    rows = residual_operator()
    assert rank(rows) == 12
    assert restrict(rows, (1, 10), (1, 10), (1, 12), (1, 12)) == [1, 2, 4, 8]
    two_letter = []
    dual_checks = 0
    for basis in subspaces(2):
        for dual in duals(basis):
            dual_checks += 1
            if restrict(rows, basis, basis, dual, dual) == [1, 2, 4, 8]:
                two_letter.append({'input_columns': basis, 'output_rows': dual})
    assert dual_checks == 560 and len(two_letter) == 13
    spaces = subspaces(3)
    assert len(spaces) == 15
    output_pairs = [(first, second, [kron(a, b) for a in first for b in second])
                    for first, second in product(spaces, repeat=2)]
    checks = 0
    for first, second in product(spaces, repeat=2):
        images = [apply(rows, kron(a, b)) for a in first for b in second]
        for _, _, outputs in output_pairs:
            restricted = [sum(((row & image).bit_count() & 1) << i
                              for i, image in enumerate(images)) for row in outputs]
            assert not is_invertible_product(restricted, 3)
            checks += 1
    assert checks == 50625
    return {'operator_rank': 12, 'two_letter_same_dual_checks': dual_checks,
            'two_letter_same_dual_codes': two_letter,
            'three_letter_four_subspace_checks': checks,
            'three_letter_binary_restrictions': 0}


def one_particle(m):
    kernel = sparse_kernel(tuple(range(m)), m)
    return [row[:2 * m] for row in kernel[2 * m:]]


def code_maps(label, m):
    bulk = 1 << (m - 2)
    body = (label & (bulk - 1)) << 1
    if label < bulk:
        return [(body, 1)], [body]
    return [(body | 1, 3), (body | 1 | (1 << (m - 1)), 1)], \
           [body | (1 << (m - 1)), body | 1 | (1 << (m - 1))]


def build(m=4):
    assert m >= 2 and m & (m - 1) == 0
    matrix = one_particle(m)

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i = bit.bit_length() - 1
        result = 0
        remaining = columns
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            j = cb.bit_length() - 1
            if matrix[i][j]:
                result ^= polynomial_multiply(matrix[i][j], minor(rows ^ bit, columns ^ cb))
        return result

    q = 1 << (m - 1)
    values = [0] * q ** 4
    perturbation = [0] * q ** 4
    top = m - 1
    divisions = 0
    for x, y, u, v in product(range(q), repeat=4):
        g = h = 0
        for (xx, bx), yy, (uu, bu), vv in product(code_maps(x, m)[0], code_maps(y, m)[1],
                                                  code_maps(u, m)[0], code_maps(v, m)[1]):
            factor = polynomial_multiply(bx, bu)
            rows, columns = vv | (yy << m), uu | (xx << m)
            g ^= polynomial_multiply(factor, minor(rows, columns))
            if xx & 1 and uu & 1 and yy >> top & 1 and vv >> top & 1:
                numerator = polynomial_multiply(
                    factor, minor(rows ^ (1 << top) ^ (1 << (m + top)),
                                  columns ^ 1 ^ (1 << m)))
                h ^= polynomial_exact_divide(numerator, 5)
                divisions += bool(numerator)
        index = x + q * y + q * q * u + q ** 3 * v
        values[index], perturbation[index] = g ^ h, h
    height = [sum(a + 1 for a in range(m - 2) if label >> a & 1)
              + (m - 1) * (label >> (m - 2)) for label in range(q)]
    return {'q': q, 'original_modes': m, 'degree': max(p.bit_length() - 1 for p in values),
            'endpoint': 1, 'index_order': 'x+q*y+q^2*u+q^3*v', 'values': values,
            'logical_energy': height, 'input_code': ['vacuum', '(1+w)low+both'],
            'output_code': ['vacuum*', 'top*+both*']}, perturbation, divisions


def independent_exterior_build(m):
    """Expand quadratic factors directly, independently of determinant recursion."""
    matrix = one_particle(m)
    polynomial = {0: 1}
    for i, j in product(range(2 * m), repeat=2):
        value = matrix[i][j]
        if not value:
            continue
        row_bit = 3 * m + i if i < m else i
        column_bit = 2 * m + j if j < m else j - m
        edge = (1 << row_bit) | (1 << column_bit)
        added = [(mask | edge, polynomial_multiply(p, value))
                 for mask, p in polynomial.items() if not mask & edge]
        for mask, p in added:
            polynomial[mask] = polynomial.get(mask, 0) ^ p
            if not polynomial[mask]:
                del polynomial[mask]
    quartic = 1 | (1 << (2 * m - 1)) | (1 << (2 * m)) | (1 << (4 * m - 1))
    numerator = {mask: polynomial_multiply(5, p) for mask, p in polynomial.items()}
    for mask, p in polynomial.items():
        if not mask & quartic:
            numerator[mask | quartic] = numerator.get(mask | quartic, 0) ^ p
    q = 1 << (m - 1)
    values = [0] * q ** 4
    for mask, p in numerator.items():
        logical = 0
        low_factors = 0
        for leg in range(4):
            old = (mask >> (m * leg)) & ((1 << m) - 1)
            low, top = old & 1, old >> (m - 1)
            if leg in (0, 2):
                if top and not low:
                    break
                core = low
                low_factors += low and not top
            else:
                if low and not top:
                    break
                core = top
            new = ((old >> 1) & ((1 << (m - 2)) - 1)) | (core << (m - 2))
            logical |= new << ((m - 1) * leg)
        else:
            for _ in range(low_factors):
                p ^= p << 1
            values[logical] ^= p
    return [polynomial_exact_divide(p, 5) for p in values]


def endpoint_checks(data, perturbation):
    q, values, height = data['q'], data['values'], data['logical_energy']
    endpoint = [p.bit_count() & 1 for p in values]
    positive = 0
    diagonal_perturbations = []

    def index(x, y, u, v):
        return x + q * y + q * q * u + q ** 3 * v

    for x, y, u, v in product(range(q), repeat=4):
        i = index(x, y, u, v)
        delta = height[y] + height[v] - height[x] - height[u]
        if endpoint[i]:
            assert delta >= 0
        if delta == 0:
            assert endpoint[i] == int(x == v and y == u)
        positive += endpoint[i] and delta > 0
        if x == y == u == v and perturbation[i]:
            diagonal_perturbations.append(x)
    words = defaultdict(list)
    for a, b in product(range(q), repeat=2):
        counts = tuple(((a >> i) & 1) + ((b >> i) & 1)
                       for i in range(q.bit_length() - 1))
        words[counts].append((a, b))
    checks = 0
    for group in words.values():
        for x, y, u, v in product(group, repeat=4):
            actual = endpoint[index(x[0], y[0], u[0], v[0])] & endpoint[index(x[1], y[1], u[1], v[1])]
            assert actual == int(x == v and y == u)
            checks += 1
    return {'alphabet': q, 'degree': data['degree'], 'local_entries': q ** 4,
            'strictly_positive_endpoint_entries': positive,
            'two_copy_fixed_mode_types': len(words), 'two_copy_endpoint_checks': checks,
            'nonzero_quartic_entries': sum(bool(p) for p in perturbation),
            'quartic_nonzero_diagonal_labels': diagonal_perturbations}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--seed-output', type=Path)
    args = parser.parse_args()
    records = []
    seed = None
    for m in (2, 4):
        data, perturbation, divisions = build(m)
        assert independent_exterior_build(m) == data['values']
        record = endpoint_checks(data, perturbation)
        record['nonzero_termwise_polynomial_divisions'] = divisions
        record['independent_exterior_coefficients_checked'] = len(data['values'])
        records.append(record)
        if m == 4:
            seed = data
    report = {'status': 'exact polynomial and typed-endpoint certificates; candidate only',
              'field': 2, 'direct_code_controls': records, 'residual_classification': classify_residual(),
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest()
                                for name in ('tools/check_residual_quartic_codes.py',
                                             'tools/check_schur_downset_crossing.py',
                                             'tools/weighted_smith_oracle.py')}}
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(rendered)
    if args.seed_output:
        args.seed_output.write_text(json.dumps(seed, separators=(',', ':')) + '\n')
    print(rendered)


if __name__ == '__main__':
    main()
