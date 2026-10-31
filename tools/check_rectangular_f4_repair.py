#!/usr/bin/env python3
"""Exact mixed-field partial repairs of the rectangular residual crossing.

The constructor and its endpoint checks use exact F4 polynomial arithmetic.
No pressure or improved factoring bound is asserted.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement, permutations, product
import json
from pathlib import Path


ZERO = (0, 0)
ONE = (1, 0)
F0_BITS = (1153, 66, 4, 8, 528, 32, 64, 1152, 256, 512, 1025, 2048)
HX = (0, 0, 1, 1)
HU = (0, 0, 1)
FA = ((0, 9), (3, 10), (11, 8))
FB = ((8, 11), (9, 10), (10, 3))
FAB = ((8, 8), (9, 9), (10, 10))


def fmul(a, b):
    """F4 = F2[omega]/(omega^2+omega+1), encoded by 0,1,2,3."""
    a0, a1, b0, b1 = a & 1, a >> 1, b & 1, b >> 1
    return ((a0*b0) ^ (a1*b1)) | (((a0*b1) ^ (a1*b0) ^ (a1*b1)) << 1)


def finv(a):
    if not a:
        raise ZeroDivisionError('Zero has no inverse in F4.')
    return fmul(a, a)


def binary_mul(a, b):
    answer = 0
    while b:
        bit = b & -b
        answer ^= a << (bit.bit_length() - 1)
        b ^= bit
    return answer


def binary_power(base, exponent):
    answer = 1
    while exponent:
        if exponent & 1:
            answer = binary_mul(answer, base)
        base = binary_mul(base, base)
        exponent >>= 1
    return answer


def binary_divide(a, b):
    if not b:
        raise ZeroDivisionError('Polynomial denominator is zero.')
    answer = 0
    while a and a.bit_length() >= b.bit_length():
        shift = a.bit_length() - b.bit_length()
        a ^= b << shift
        answer ^= 1 << shift
    if a:
        raise ArithmeticError('Polynomial division is not exact.')
    return answer


def add(a, b):
    return a[0] ^ b[0], a[1] ^ b[1]


def scale(a, scalar):
    if scalar == 0:
        return ZERO
    if scalar == 1:
        return a
    if scalar == 2:
        return a[1], a[0] ^ a[1]
    if scalar == 3:
        return a[0] ^ a[1], a[0]
    raise ValueError('Not an F4 coefficient.')


def multiply(a, b):
    low = binary_mul(a[0], b[0])
    high = binary_mul(a[1], b[1])
    cross = binary_mul(a[0] ^ a[1], b[0] ^ b[1])
    return low ^ high, cross ^ low


def times_binary(a, b):
    return binary_mul(a[0], b), binary_mul(a[1], b)


def divide_binary(a, b):
    return binary_divide(a[0], b), binary_divide(a[1], b)


def evaluate_one(a):
    return (a[0].bit_count() & 1) | ((a[1].bit_count() & 1) << 1)


def c_order(a):
    if a == ZERO:
        return None
    count = 0
    while not evaluate_one(a):
        a = divide_binary(a, 3)
        count += 1
    return count


def ramify_binary(a, top):
    answer = 0
    while a:
        bit = a & -a
        exponent = bit.bit_length() - 1
        a ^= bit
        sub = exponent
        while True:
            answer ^= 1 << (top * sub)
            if not sub:
                break
            sub = (sub - 1) & exponent
    return answer


def matrix_multiply(left, right):
    answer = [[0]*len(right[0]) for _ in left]
    for i, row in enumerate(left):
        for k, a in enumerate(row):
            if a:
                for j, b in enumerate(right[k]):
                    answer[i][j] ^= fmul(a, b)
    return answer


def residual(alpha, beta):
    result = [[row >> j & 1 for j in range(12)] for row in F0_BITS]
    for coefficient, entries in ((alpha, FA), (beta, FB), (fmul(alpha, beta), FAB)):
        for i, j in entries:
            result[i][j] ^= coefficient
    return result


def partial_repair(alpha, beta):
    """The explicit block elimination, valid exactly when 1+alpha*beta != 0."""
    if alpha not in range(4) or beta not in range(4):
        raise ValueError('alpha and beta must be F4 coefficients.')
    delta = 1 ^ fmul(alpha, beta)
    if not delta:
        raise ValueError('The middle height block is singular: alpha*beta=1.')
    result = [[int(i == j) for j in range(12)] for i in range(12)]
    result[8][8] = result[9][9] = finv(delta)
    result[10] = [0]*12
    for column in (0, 7, 10):
        result[10][column] = delta
    result[10][3] = fmul(beta, delta)
    result[10][9] = alpha
    result[11] = [0]*12
    result[11][8], result[11][11] = alpha, delta
    return result


def check_literal_residual(alpha, beta):
    """Independent seven-factor exterior endpoint, before block elimination."""
    factors = (65, 130, 20, 40, 96, 6, 153)
    support = {0}
    for factor in factors:
        support ^= {mask | factor for mask in support if not mask & factor}
    xmap = (1, 0, 2, 3)
    umap = (((1, 1),), ((0, 1),), ((2, 1), (3, alpha)))
    ymap = (((1, 1),), ((0, 1),), ((2, 1), (3, beta)))
    literal = [[0]*12 for _ in range(12)]
    for x, y, u, v in product(range(4), range(3), range(3), range(4)):
        value = 0
        for (yy, cy), (uu, cu) in product(ymap[y], umap[u]):
            mask = xmap[x] | (yy << 2) | (uu << 4) | (xmap[v] << 6)
            if mask in support:
                value ^= fmul(cy, cu)
        literal[v+4*y][x+4*u] = value
    assert literal == residual(alpha, beta)
    return literal


def check_repair(alpha, beta):
    original = check_literal_residual(alpha, beta)
    correction = partial_repair(alpha, beta)
    endpoint = matrix_multiply(correction, original)
    heights = [HX[i % 4] + HU[i // 4] for i in range(12)]
    for row in range(12):
        assert endpoint[row][row] == 1
        for column in range(12):
            if row != column and endpoint[row][column]:
                assert heights[row] < heights[column]
    assert correction[0] == [int(j == 0) for j in range(12)]
    return original, correction, endpoint


def check_symbolic_repair():
    """Verify delta*C*F = delta*U in F2[alpha,beta], not only at F4 points."""
    zero, one = frozenset(), frozenset({(0, 0)})
    a, b = frozenset({(1, 0)}), frozenset({(0, 1)})
    def pmul(left, right):
        answer = set()
        for i, j in left:
            for k, l in right:
                monomial = i+k, j+l
                if monomial in answer:
                    answer.remove(monomial)
                else:
                    answer.add(monomial)
        return frozenset(answer)
    delta = one ^ pmul(a, b)
    delta2 = pmul(delta, delta)
    f = [[one if row >> j & 1 else zero for j in range(12)] for row in F0_BITS]
    for coefficient, entries in ((a, FA), (b, FB), (pmul(a, b), FAB)):
        for i, j in entries:
            f[i][j] ^= coefficient
    scaled_c = [[delta if i == j else zero for j in range(12)] for i in range(12)]
    scaled_c[8][8] = scaled_c[9][9] = one
    scaled_c[10] = [zero]*12
    for j in (0, 7, 10):
        scaled_c[10][j] = delta2
    scaled_c[10][3], scaled_c[10][9] = pmul(b, delta2), pmul(a, delta)
    scaled_c[11] = [zero]*12
    scaled_c[11][8], scaled_c[11][11] = pmul(a, delta), delta2
    scaled_u = [[delta if i == j else zero for j in range(12)] for i in range(12)]
    for i, j in ((0, 7), (0, 10), (1, 6), (4, 9), (7, 10)):
        scaled_u[i][j] = delta
    for i, j in ((0, 9), (3, 10)):
        scaled_u[i][j] = pmul(a, delta)
    for i, j in ((8, 11), (9, 10)):
        scaled_u[i][j] = b
    for i, j in product(range(12), repeat=2):
        answer = zero
        for k in range(12):
            answer ^= pmul(scaled_c[i][k], f[k][j])
        assert answer == scaled_u[i][j]
    return {'identity': 'delta*C(alpha,beta)*F(alpha,beta)=delta*U(alpha,beta)',
            'coefficient_ring': 'F2[alpha,beta]', 'polynomial_entries_checked': 144,
            'hypothesis_for_division': 'delta=1+alpha*beta != 0'}


def check_arithmetic():
    """Compare bit-plane polynomial multiplication with scalar convolution."""
    def direct_field_product(a, b):
        raw = 0
        for i, j in product(range(2), repeat=2):
            if (a >> i & 1) and (b >> j & 1):
                raw ^= 1 << (i+j)
        if raw & 4:
            raw ^= 7
        return raw
    def planes(coefficients):
        return (sum((a & 1) << i for i, a in enumerate(coefficients)),
                sum((a >> 1) << i for i, a in enumerate(coefficients)))
    polynomials = list(product(range(4), repeat=3))
    for left, right in product(polynomials, repeat=2):
        expected = [0]*5
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                expected[i+j] ^= direct_field_product(a, b)
        assert multiply(planes(left), planes(right)) == planes(expected)
    return 4096


def one_particle(indices, ambient):
    indices = tuple(indices)
    n = len(indices)
    if (n < 2 or indices != tuple(sorted(set(indices))) or indices[0] != 0
            or ambient <= indices[-1] or ambient & (ambient-1)):
        raise ValueError('Use sorted distinct modes containing zero and a larger power-of-two ambient.')
    answer = [[0]*(2*n) for _ in range(2*n)]
    for i, a in enumerate(indices):
        for j, b in enumerate(indices):
            answer[i][j] = binary_mul(3, binary_power(5, a)) if a & b == a else 0
            answer[i][n+j] = (1 << (2*ambient)) if a & b == b else 0
            answer[n+i][j] = int(i == j)
            answer[n+i][n+j] = 3 << (2*(ambient-1-(a | b)))
    return answer


def build(indices=(0, 1, 2), ambient=4, alpha=1, beta=2, check_types=True,
          independent_exterior=False):
    indices = tuple(indices)
    original, correction, endpoint = check_repair(alpha, beta)
    n, top_index = len(indices), indices[-1]
    bulk, top = 1 << (n-2), 1 << (n-1)
    qx, qu = 4*bulk, 3*bulk
    matrix = one_particle(indices, ambient)

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i, answer = bit.bit_length()-1, 0
        remaining = columns
        while remaining:
            column = remaining & -remaining
            remaining ^= column
            j = column.bit_length()-1
            if matrix[i][j]:
                answer ^= binary_mul(matrix[i][j], minor(rows ^ bit, columns ^ column))
        return answer

    xmap = [((a % bulk) << 1) | (1, 0, top, top | 1)[a // bulk] for a in range(qx)]
    def mixed_map(coefficient):
        result = []
        for a in range(qu):
            low = (a % bulk) << 1
            entries = (((1, 1),), ((0, 1),), ((top, 1), (top | 1, coefficient)))[a // bulk]
            result.append(tuple((low | mask, value) for mask, value in entries if value))
        return result
    umap, ymap = mixed_map(alpha), mixed_map(beta)
    def position(x, y, u, v):
        return x + qx*(y + qu*(u + qu*v))
    values = [ZERO]*(qx*qx*qu*qu)
    numerator = [ZERO]*len(values)
    raw_support = None
    if independent_exterior:
        raw_support = {0: 1}
        for i, row in enumerate(matrix):
            output = 3*n+i if i < n else i
            for j, coefficient in enumerate(row):
                if not coefficient:
                    continue
                input_index = 2*n+j if j < n else j-n
                monomial = (1 << output) | (1 << input_index)
                for mask, value in list(raw_support.items()):
                    if mask & monomial:
                        continue
                    target = mask | monomial
                    raw_support[target] = raw_support.get(target, 0) ^ binary_mul(coefficient, value)
                    if not raw_support[target]:
                        del raw_support[target]
        quartic = 1 | (1 << (2*n-1)) | (1 << (2*n)) | (1 << (4*n-1))
        exterior = {mask: binary_mul(5, value) for mask, value in raw_support.items()}
        for mask, value in raw_support.items():
            if not mask & quartic:
                target = mask | quartic
                exterior[target] = exterior.get(target, 0) ^ value

    for x, y, u, v in product(range(qx), range(qu), range(qu), range(qx)):
        xx, vv, answer, control = xmap[x], xmap[v], ZERO, ZERO
        for (yy, cy), (uu, cu) in product(ymap[y], umap[u]):
            p = binary_mul(5, minor(vv | (yy << n), uu | (xx << n)))
            if xx & 1 and uu & 1 and yy & top and vv & top:
                p ^= minor((vv ^ top) | ((yy ^ top) << n), (uu ^ 1) | ((xx ^ 1) << n))
            coefficient = fmul(cy, cu)
            answer = add(answer, scale((p, 0), coefficient))
            if independent_exterior:
                mask = xx | (yy << n) | (uu << (2*n)) | (vv << (3*n))
                control = add(control, scale((exterior.get(mask, 0), 0), coefficient))
        if independent_exterior:
            assert answer == control
        numerator[position(x, y, u, v)] = answer

    gate = []
    for row in range(12):
        ht = HX[row % 4] + HU[row // 4]
        entries = []
        for column, coefficient in enumerate(correction[row]):
            if coefficient:
                vs, ys = column % 4, column // 4
                exponent = 2 + HX[vs] + HU[ys] - ht
                assert exponent >= 0
                entries.append((vs, ys, coefficient, binary_power(3, exponent)))
        gate.append(entries)
    for x, y, u, v in product(range(qx), range(qu), range(qu), range(qx)):
        answer = ZERO
        for vs, ys, coefficient, p in gate[v // bulk + 4*(y // bulk)]:
            source = numerator[position(x, ys*bulk + y % bulk, u, vs*bulk + v % bulk)]
            answer = add(answer, scale(times_binary(source, p), coefficient))
        values[position(x, y, u, v)] = answer

    cancellation = min(4, min(c_order(value) for value in values if value != ZERO))
    denominator_degree = 4-cancellation
    assert denominator_degree <= 1  # Uniform residue cancellation proved in the note.
    if cancellation:
        denominator = binary_power(3, cancellation)
        values = [divide_binary(value, denominator) for value in values]
    hx = [sum(indices[j+1] for j in range(n-2) if a >> j & 1) + top_index*HX[a // bulk]
          for a in range(qx)]
    hu = [sum(indices[j+1] for j in range(n-2) if a >> j & 1) + top_index*HU[a // bulk]
          for a in range(qu)]
    for x, y, u, v in product(range(qx), range(qu), range(qu), range(qx)):
        p = tuple(ramify_binary(plane, top_index) for plane in values[position(x, y, u, v)])
        shift = hu[y]+hx[v]-hx[x]-hu[u]-denominator_degree*top_index
        if shift < 0:
            assert not ((p[0] | p[1]) & ((1 << (-shift))-1))
            p = p[0] >> (-shift), p[1] >> (-shift)
        else:
            p = p[0] << shift, p[1] << shift
        expected = endpoint[v // bulk+4*(y // bulk)][x // bulk+4*(u // bulk)]
        expected *= int(x % bulk == v % bulk and y % bulk == u % bulk)
        assert (p[0] & 1) | ((p[1] & 1) << 1) == expected

    checks = 0
    if check_types:
        denominator = binary_power(3, 2*denominator_degree)
        xtypes = [sorted(set(permutations(pair))) for pair in combinations_with_replacement(range(qx), 2)]
        utypes = [sorted(set(permutations(pair))) for pair in combinations_with_replacement(range(qu), 2)]
        for X, U in product(xtypes, utypes):
            for x, v, y, u in product(X, X, U, U):
                a = values[position(x[0], y[0], u[0], v[0])]
                b = values[position(x[1], y[1], u[1], v[1])]
                polynomial = divide_binary(multiply(a, b), denominator)
                assert evaluate_one(polynomial) == int(x == v and y == u)
                checks += 1
    histogram = Counter(c_order(value) for value in values if value != ZERO)
    return {'field': 4, 'field_polynomial': 'omega^2+omega+1', 'indices': list(indices),
            'ambient': ambient, 'alpha': alpha, 'beta': beta, 'delta': 1 ^ fmul(alpha, beta),
            'dimensions': [qx, qu, qu, qx], 'q_x': qx, 'q_u': qu, 'bulk_size': bulk,
            'repair_kind': 'mixed_height_partial', 'degree': max(max(value).bit_length()-1 for value in values),
            'denominator_degree': denominator_degree, 'cancelled_c_order': cancellation,
            'strict_endpoint': False, 'typed_endpoint': 'independent common X/V and U/Y full-letter types',
            'logical_heights_x_v': hx, 'logical_heights_u_y': hu,
            'residual_matrix': original, 'repair_matrix': correction, 'ramified_endpoint_matrix': endpoint,
            'ramified_endpoint_checks': len(values), 'two_copy_type_checks': checks,
            'independent_exterior_checks': len(values) if independent_exterior else 0,
            'nonzero_entries': sum(value != ZERO for value in values),
            'c_order_histogram': dict(sorted(histogram.items())),
            'coefficient_sha256': sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest(),
            'values': values}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--seed-directory', type=Path)
    parser.add_argument('--full', action='store_true', help='Include all four-mode two-copy type checks.')
    args = parser.parse_args()
    symbolic = check_symbolic_repair()
    arithmetic_checks = check_arithmetic()
    valid, rejected = [], []
    for alpha, beta in product(range(4), repeat=2):
        check_literal_residual(alpha, beta)
        if fmul(alpha, beta) == 1:
            try:
                partial_repair(alpha, beta)
            except ValueError:
                rejected.append([alpha, beta])
            else:
                raise AssertionError('Singular height block was accepted.')
        else:
            check_repair(alpha, beta)
            valid.append([alpha, beta])
    cases = []
    for indices, ambient in (((0, 1), 2), ((0, 1, 2), 4), ((0, 1, 2, 3), 4), ((0, 2, 5), 8)):
        data = build(indices, ambient, check_types=args.full or len(indices) <= 3,
                     independent_exterior=len(indices) <= 3)
        cases.append({key: value for key, value in data.items() if key != 'values'})
        if args.seed_directory:
            args.seed_directory.mkdir(parents=True, exist_ok=True)
            name = 'rectangular-f4-12-'+'-'.join(map(str, indices))+f'-M{ambient}.json'
            (args.seed_directory/name).write_text(json.dumps(data, separators=(',', ':'))+'\n')
    result = {'status': 'exact mixed-F4 partial repair and typed crossing checks; no pressure improvement claimed',
              'symbolic_repair': symbolic, 'independent_polynomial_products_checked': arithmetic_checks,
              'valid_parameter_pairs': valid, 'singular_parameter_pairs_rejected': rejected, 'cases': cases}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
