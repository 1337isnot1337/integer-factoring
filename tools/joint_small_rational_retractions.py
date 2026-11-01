#!/usr/bin/env python3
"""Rational joint retraction module for the five-to-four-letter crossing.

Uses the saturation, predictable-degree and right-inverse method of
joint_hot_retraction_module.py, applied to a different 9-by-16 row space.
"""
from functools import lru_cache
from itertools import product

from joint_five_letter_splits import HEIGHTS, coordinate_start, source
from joint_hot_retraction_module import degree, dependency, divide_linear, polynomial_multiply

H = tuple(a+b for a in HEIGHTS for b in HEIGHTS)


@lru_cache(None)
def source_rows():
    matrix, _, _ = source()
    selected = tuple(row.index(1) for row in coordinate_start()[0])
    excluded = tuple(i for i in range(25) if i not in selected)
    assert all(high == 0 for row in matrix for _, high in row)
    return tuple(tuple(matrix[i][j][0] for j in selected) for i in excluded)


@lru_cache(None)
def saturated_basis():
    rows = [list(row) for row in source_rows()]
    steps = []
    for root in (0, 1):
        while True:
            evaluations = [sum(((v & 1) if root == 0 else v.bit_count() % 2) << j
                               for j, v in enumerate(row)) for row in rows]
            relation = dependency(evaluations)
            if not relation:
                break
            selected = [i for i in range(len(rows)) if relation >> i & 1]
            pivot = max(selected, key=lambda i: degree(rows[i]))
            new = [0]*16
            for i in selected:
                new = [a ^ b for a, b in zip(new, rows[i])]
            if any(new):
                rows[pivot] = [divide_linear(v, root) for v in new]
            else:
                rows.pop(pivot)
            steps.append(dict(root=root, rows=selected, pivot=pivot))
    reductions = 0
    while True:
        ds = [degree(row) for row in rows]
        relation = dependency([sum(((v >> d) & 1) << j for j, v in enumerate(row))
                               for row, d in zip(rows, ds)])
        if not relation:
            break
        selected = [i for i in range(len(rows)) if relation >> i & 1]
        pivot = max(selected, key=lambda i: ds[i])
        new = [0]*16
        for i in selected:
            new = [a ^ (b << (ds[pivot]-ds[i])) for a, b in zip(new, rows[i])]
        if any(new):
            assert degree(new) < ds[pivot]
            rows[pivot] = new
        else:
            rows.pop(pivot)
        reductions += 1
    return rows, dict(saturation_steps=steps, leading_reductions=reductions)


@lru_cache(None)
def row_space(row_height, protected):
    """All strict bounded corrections fixing the first protected raw jets."""
    if row_height not in range(9) or protected not in range(7):
        raise ValueError('Use a row height from zero to eight and zero to six protected jets.')
    basis, _ = saturated_basis()
    strict = [[v ^ (v << 1) for v in row] for row in basis]
    lower = 8+protected-row_height
    upper = [8+h-protected for h in H]
    width = max(upper)+1
    allowed = [((1 << (u+1))-1) ^ ((1 << max(0, lower))-1) if u >= lower else 0 for u in upper]
    generators = [[v << shift for v in row] for row in strict
                  for shift in range(max(upper)-degree(row)+1)]
    pivots, output = {}, []
    for i, generator in enumerate(generators):
        image = sum((v & ~mask) << (j*width) for j, (v, mask) in enumerate(zip(generator, allowed)))
        combo = 1 << i
        while image:
            k = image.bit_length()-1
            if k in pivots:
                image ^= pivots[k][0]
                combo ^= pivots[k][1]
            else:
                pivots[k] = image, combo
                break
        if not image:
            row = [0]*16
            for j, g in enumerate(generators):
                if combo >> j & 1:
                    row = [a ^ b for a, b in zip(row, g)]
            assert any(row) and all(not (v & ~m) for v, m in zip(row, allowed))
            output.append(row)
    return output


def regular_zero(matrix):
    result = []
    for h, row in zip(H, matrix):
        if h >= 8:
            result.append([p << (h-8) for p in row])
        else:
            assert all(not (p & ((1 << (8-h))-1)) for p in row)
            result.append([p >> (8-h) for p in row])
    return result


def regular_infinity(matrix):
    result = [[0]*16 for _ in range(16)]
    for x, y, u, v in product(range(4), repeat=4):
        i, j = 4*x+y, 4*v+u
        bound, p, value = 8+H[j], matrix[i][j], 0
        for k in range(p.bit_length()):
            if p >> k & 1:
                assert k <= bound
                value ^= 1 << (bound-k)
        result[4*x+u][4*v+y] = value
    return result


def two_pivot_flag(matrix):
    """Exact truncated binary series elimination, recording constant gauges."""
    precision, mask = 8, 255
    a = [[p & mask for p in row] for row in matrix]
    left = [[int(i == j) for j in range(16)] for i in range(16)]
    right = [row[:] for row in left]
    for valuation, r, c in ((0, 0, 0), (1, 1, 1)):
        assert (a[r][c] & -a[r][c]).bit_length()-1 == valuation
        unit, inverse = a[r][c] >> valuation, 1
        for k in range(1, precision):
            if polynomial_multiply(unit, inverse) >> k & 1:
                inverse ^= 1 << k
        a[r] = [polynomial_multiply(p, inverse) & mask for p in a[r]]
        assert a[r][c] == 1 << valuation
        for i in range(r+1, 16):
            assert not (a[i][c] & ((1 << valuation)-1))
            factor = a[i][c] >> valuation
            a[i] = [(p ^ polynomial_multiply(factor, q)) & mask for p, q in zip(a[i], a[r])]
            if factor & 1:
                left[i] = [p ^ q for p, q in zip(left[i], left[r])]
        for j in range(c+1, 16):
            assert not (a[r][j] & ((1 << valuation)-1))
            factor = a[r][j] >> valuation
            for i in range(16):
                a[i][j] ^= polynomial_multiply(factor, a[i][c]) & mask
            if factor & 1:
                for i in range(16):
                    right[i][j] ^= right[i][c]
    assert all(not (a[i][j] & 7) for i, j in product(range(2, 16), repeat=2))
    return left, right


def project(matrix, flag):
    left, right = flag
    temporary = [[sum(left[i][k]*matrix[k][j] for k in range(16)) % 2
                  for j in range(16)] for i in range(16)]
    return sum((sum(temporary[i][k]*right[k][j] for k in range(16)) % 2)
               << (14*(i-2)+j-2) for i, j in product(range(2, 16), repeat=2))


def jet(row, values, order, side):
    matrix = [[0]*16 for _ in range(16)]
    for j, p in enumerate(values):
        if side == 0:
            matrix[row][j] = (p >> (8+order-H[row])) & 1
        else:
            x, y = divmod(row, 4)
            v, u = divmod(j, 4)
            matrix[4*x+u][4*v+y] = (p >> (8+H[j]-order)) & 1
    return matrix
