#!/usr/bin/env python3
"""Exact integral spin-11 seeds and polynomial inverse witness generation."""
from itertools import product
from math import gcd
from functools import reduce

from check_modular_spin_crossings import rational_casimir_crossing
from check_all_field_pascal_crossings import trim, exact_divide

Q, N = 12, 11
PRIMES = (3, 5, 7, 11)


def polynomial_gcd(a, b, p):
    while b:
        r = list(a)
        while len(r) >= len(b):
            shift = len(r)-len(b)
            coefficient = r[-1]*pow(b[-1], -1, p) % p
            for j, value in enumerate(b):
                r[j+shift] = (r[j+shift]-coefficient*value) % p
            while r and not r[-1]:
                r.pop()
        a, b = b, tuple(r)
    return trim([x*pow(a[-1], -1, p) for x in a], p) if a else ()


def integer_source():
    rational, denominator = rational_casimir_crossing(N)
    assert denominator == 39916800
    assert all((x*denominator).denominator == 1
               for value in rational.values() for x in value)
    integer = {key: tuple(int(x*denominator) for x in value)
               for key, value in rational.items()}
    assert reduce(gcd, (abs(x) for v in integer.values() for x in v)) == 1
    assert all(v[0] == denominator*(k[0] == k[3] and k[1] == k[2])
               and v[-1] == (k[0] == k[2] and k[1] == k[3])
               for k, v in integer.items())
    return integer


def primitive_source(integer, p):
    reduced = {key: trim(value, p) for key, value in integer.items()}
    common = ()
    for value in reduced.values():
        common = polynomial_gcd(common, value, p)
    primitive = {key: exact_divide(value, common, p)
                 for key, value in reduced.items()}
    return common, primitive


def cut_blocks(primitive, cut):
    # T_xyuv=F_xyvu: the leading identity of F becomes a strict crossing.
    def entry(x, y, u, v):
        return primitive.get((x, y, v, u), ())
    if cut == 'row':
        for occupation in range(2*N+1):
            labels = range(max(0, occupation-N), min(N, occupation)+1)
            yield occupation, [[entry(x, occupation-x, u, occupation-u)
                                for u in labels] for x in labels]
    else:
        assert cut == 'opposite'
        for d in range(-N, N+1):
            left = range(max(0, -d), min(N, N-d)+1)
            right = range(max(0, d), min(N, N+d)+1)
            yield d, [[entry(u+d, v-d, u, v) for v in right] for u in left]


def matrix_product(a, b, p):
    size = len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(size)) % p
             for j in range(size)] for i in range(size)]


def inverse_witness(matrix, degree, p):
    """Generate A with M(u)A(u)=(u^p-u)^2 I by a reversed series."""
    size = len(matrix)
    coefficients = [[[value[k] if k < len(value) else 0 for value in row]
                     for row in matrix] for k in range(degree+1)]
    leading = coefficients[-1]
    assert all(sum(row) == 1 for row in leading)
    assert all(sum(leading[i][j] for i in range(size)) == 1 for j in range(size))
    assert all(x in (0, 1) for row in leading for x in row)
    inverse_leading = [list(row) for row in zip(*leading)]
    maximum = 2*p-degree
    series = [inverse_leading]
    for k in range(1, maximum+1):
        total = [[0]*size for _ in range(size)]
        for j in range(1, min(k, degree)+1):
            term = matrix_product(coefficients[degree-j], series[k-j], p)
            total = [[(x+y) % p for x, y in zip(a, b)] for a, b in zip(total, term)]
        series.append([[-x % p for x in row]
                       for row in matrix_product(inverse_leading, total, p)])
    reversed_q = {0: 1, p-1: -2, 2*p-2: 1}
    answer = [[()] * size for _ in range(size)]
    for i, j in product(range(size), repeat=2):
        reversed_a = [sum(c*series[k-h][i][j]
                          for h, c in reversed_q.items() if h <= k) % p
                      for k in range(maximum+1)]
        answer[i][j] = trim(reversed_a[::-1], p)
    return answer


def generate():
    integer = integer_source()
    cases = []
    for p in PRIMES:
        common, primitive = primitive_source(integer, p)
        degree = max(len(value)-1 for value in primitive.values())
        case = dict(prime=p, common_polynomial=common, degree=degree)
        if p != 3:
            case['witnesses'] = [dict(cut=cut, block=key,
                                     inverse=inverse_witness(matrix, degree, p))
                                 for cut in ('row', 'opposite')
                                 for key, matrix in cut_blocks(primitive, cut)]
        cases.append(case)
    return dict(spin=N, alphabet=Q, denominator=39916800, cases=cases)
