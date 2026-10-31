#!/usr/bin/env python3
"""Odd cyclic crossings from arbitrary endpoint matrix polynomials.

The two cuts are a diagonal symbol vector and q repeated matrix blocks.
All algebra and profile computations are exact; pressure minimization is a
numerical search aid. See research/odd-cyclic-matrix-crossings.md.
"""
from __future__ import annotations

from fractions import Fraction
from math import exp, log


def order_root(q, prime):
    if q < 3 or q % 2 == 0 or (prime - 1) % q:
        raise ValueError("need odd q>=3 dividing prime-1")
    for value in range(2, prime):
        if pow(value, q, prime) == 1 and all(
                pow(value, j, prime) != 1 for j in range(1, q)):
            return value
    raise ValueError("no root of the required order")


def matrix_rank(matrix, prime):
    work = [[int(x) % prime for x in row] for row in matrix]
    pivots = 0
    for col in range(len(work[0]) if work else 0):
        pivot = next((r for r in range(pivots, len(work)) if work[r][col]), None)
        if pivot is None:
            continue
        work[pivots], work[pivot] = work[pivot], work[pivots]
        inv = pow(work[pivots][col], -1, prime)
        work[pivots] = [x * inv % prime for x in work[pivots]]
        for r in range(pivots + 1, len(work)):
            scale = work[r][col]
            if scale:
                work[r] = [(x - scale * y) % prime
                           for x, y in zip(work[r], work[pivots])]
        pivots += 1
        if pivots == len(work):
            break
    return pivots


def symbols(coefficients, prime, root=None):
    q = len(coefficients[0])
    root = order_root(q, prime) if root is None else root
    phase = [[pow(root, (-b * r) % q, prime) for r in range(q)]
             for b in range(q)]
    return [[sum(phase[b][(2*x-a) % q] * matrix[x][(a-x) % q]
                  for x in range(q)) % prime
             for matrix in coefficients]
            for a in range(q) for b in range(q)]


def phase_point(q, prime, a=0, b=0, scale=1):
    root = order_root(q, prime)
    inv = pow(q, -1, prime)
    return [[scale * inv * pow(root, b*(x-v) % q, prime) % prime
             if (x+v) % q == a else 0 for v in range(q)] for x in range(q)]


def tensor_entry(matrix, x, y, u, v):
    q = len(matrix)
    if (x+y-u-v) % q:
        return 0
    half = pow(2, -1, q)
    return matrix[((x+u)*half) % q][(v-(x-u)*half) % q]


def flatten(matrix, cut):
    q = len(matrix)
    pairs = [(x, y) for x in range(q) for y in range(q)]
    if cut == 0:
        return [[tensor_entry(matrix, x, y, u, v) for u, v in pairs]
                for x, y in pairs]
    if cut == 1:
        return [[tensor_entry(matrix, x, y, u, v) for y, v in pairs]
                for x, u in pairs]
    raise ValueError("cut must be zero or one")


def _mul(a, b, length, prime):
    out = [0] * length
    for i, x in enumerate(a[:length]):
        if x:
            for j, y in enumerate(b[:length-i]):
                out[i+j] = (out[i+j] + x*y) % prime
    return out


def _inverse(a, length, prime):
    out = [0] * length
    out[0] = inv = pow(a[0], -1, prime)
    for n in range(1, length):
        out[n] = -inv * sum(a[k] * out[n-k]
                            for k in range(1, min(n+1, len(a)))) % prime
    return out


def smith_orders(coefficients, prime):
    """Exact orders of a square polynomial matrix with nonzero determinant.

    Initial precision q*degree+1 exceeds every possible determinant order.
    Every common valuation division decreases the retained precision. The
    residual pivot operations are regular, so no unproved precision guess
    is made. Endpoint-invertible input matrices meet the determinant premise.
    """
    q = len(coefficients[0])
    degree = len(coefficients) - 1
    length = q*degree + 1
    matrix = [[[int(coefficients[d][i][j]) % prime if d <= degree else 0
                for d in range(length)] for j in range(q)] for i in range(q)]
    offset = 0
    orders = []
    for remaining in range(q, 0, -1):
        val, row, col = min((next((k for k, x in enumerate(matrix[i][j]) if x), length), i, j)
                            for i in range(remaining) for j in range(remaining))
        if val == length:
            raise ValueError("zero determinant or insufficient polynomial degree bound")
        matrix[0], matrix[row] = matrix[row], matrix[0]
        for line in matrix:
            line[0], line[col] = line[col], line[0]
        if val:
            matrix = [[series[val:] for series in line] for line in matrix]
            length -= val
            offset += val
        orders.append(offset)
        inverse = _inverse(matrix[0][0], length, prime)
        residual = []
        for i in range(1, remaining):
            ratio = _mul(matrix[i][0], inverse, length, prime)
            line = []
            for j in range(1, remaining):
                subtract = _mul(ratio, matrix[0][j], length, prime)
                line.append([(x-y) % prime for x, y in zip(matrix[i][j], subtract)])
            residual.append(line)
        matrix = residual
    return orders


def profiles(coefficients, prime):
    q = len(coefficients[0])
    endpoint = [[sum(matrix[i][j] for matrix in coefficients) % prime
                 for j in range(q)] for i in range(q)]
    if endpoint != [[int(i == j) for j in range(q)] for i in range(q)]:
        raise ValueError("matrix endpoint must be exactly the identity")
    sym = symbols(coefficients, prime)
    zero = [next(i for i, x in enumerate(p) if x) for p in sym]
    degrees = [max(i for i, x in enumerate(p) if x) for p in sym]
    matrix_zero = smith_orders(coefficients, prime)
    degree = len(coefficients)-1
    matrix_degrees = [degree-x for x in smith_orders(coefficients[::-1], prime)]
    return {"symbol_orders": zero, "symbol_degrees": degrees,
            "matrix_orders": matrix_zero, "matrix_degrees": matrix_degrees,
            "block_multiplicity": q}


def pressure(profile, q, orientation=0):
    if orientation == 0:
        degrees, orders = profile['matrix_degrees'], profile['symbol_orders']
    else:
        degrees, orders = profile['symbol_degrees'], profile['matrix_orders']
    slopes = [a-b for a in degrees for b in orders]

    def value(t):
        shift = max(slopes)*t
        weights = [exp(a*t-shift) for a in slopes]
        total = sum(weights)
        return log(q)+shift+log(total), sum(a*w for a, w in zip(slopes, weights))/total

    if sum(slopes) >= 0:
        tau = 0.0
    elif max(slopes) <= 0:
        tau = 100.0
    else:
        lo, hi = 0.0, 1.0
        while value(hi)[1] < 0:
            hi *= 2
        for _ in range(70):
            mid = (lo+hi)/2
            if value(mid)[1] < 0:
                lo = mid
            else:
                hi = mid
        tau = (lo+hi)/2
    logb, _ = value(tau)
    return {"rho_numerical": logb/(4*log(q)), "log_scale": tau,
            "orientation": orientation}


def rational_pressure(profile, q, scale, orientation=0):
    scale = Fraction(scale)
    if scale < 1:
        raise ValueError("scale must be at least one")
    if orientation == 0:
        degrees, orders = profile['matrix_degrees'], profile['symbol_orders']
    else:
        degrees, orders = profile['symbol_degrees'], profile['matrix_orders']
    return q * sum(scale**a for a in degrees) * sum(scale**(-b) for b in orders)
