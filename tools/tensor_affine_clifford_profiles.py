#!/usr/bin/env python3
"""Complete affine-CAR profiles through one auxiliary Gaussian mode.

An integral permutation conjugates G tensor I2 to a Gaussian operator.
Its 18-by-18 Clifford matrix determines the full 256-order active
profile. F4 polynomial one-mode maps and Laurent height gauges are
supported, with every determinant valuation and clearing charged.
"""
from functools import lru_cache
from native_matrix_smith_f4 import MatrixOracle
from tensor_literal_full_tail_core import images
from tensor_hot_compounds_f16 import mul as scalar_mul

ZERO, ONE = (0, 0), (1, 0)


def add(a, b):
    return a[0] ^ b[0], a[1] ^ b[1]


@lru_cache(maxsize=32768)
def binary_product(a, b):
    if not a or not b:
        return 0
    if a.bit_count() > b.bit_count():
        a, b = b, a
    out = 0
    while a:
        bit = a & -a
        out ^= b << (bit.bit_length()-1)
        a ^= bit
    return out


def multiply(a, b):
    if a == ZERO or b == ZERO:
        return ZERO
    if a == ONE:
        return b
    if b == ONE:
        return a
    first = binary_product(a[0], b[0])
    if not a[1] and not b[1]:
        return first, 0
    second = binary_product(a[1], b[1])
    return first ^ second, binary_product(a[0] ^ a[1], b[0] ^ b[1]) ^ first


def valuation(a):
    union = a[0] | a[1]
    if not union:
        raise ValueError('Valuation of zero is undefined here.')
    return (union & -union).bit_length()-1


def matrix_product(A, B):
    out = [[ZERO]*len(B[0]) for _ in A]
    for i, row in enumerate(A):
        for k, coefficient in enumerate(row):
            if coefficient == ZERO:
                continue
            for j, value in enumerate(B[k]):
                if value != ZERO:
                    out[i][j] = add(out[i][j], multiply(coefficient, value))
    return out


def affine(mask, degree):
    if not 0 <= mask < 256 or degree < 1:
        raise ValueError('Expected an F4 two-by-two mask and positive degree.')
    f, out = 1 ^ (1 << degree), []
    for i in range(2):
        row = []
        for j in range(2):
            coefficient = (mask >> (2*(2*i+j))) & 3
            row.append(((f if coefficient & 1 else 0) ^ int(i == j),
                        f if coefficient & 2 else 0))
        out.append(row)
    return out


def nilpotent_masks():
    out = []
    for mask in range(256):
        a, b, c, d = [(mask >> (2*i)) & 3 for i in range(4)]
        if a == d and scalar_mul(a, d) == scalar_mul(b, c):
            out.append(mask)
    assert len(out) == 16
    return out


def weighted(matrix, row_heights=None, column_heights=None, reverse=False, transpose=False):
    if transpose:
        matrix = [list(row) for row in zip(*matrix)]
    row_heights = [0, 0] if row_heights is None else row_heights
    column_heights = [0, 0] if column_heights is None else column_heights
    terms = [[[] for _ in range(2)] for _ in range(2)]
    for i, row in enumerate(matrix):
        for j, planes in enumerate(row):
            for exponent in range(max(p.bit_length() for p in planes)):
                coefficient = ((planes[0] >> exponent) & 1) | (((planes[1] >> exponent) & 1) << 1)
                if coefficient:
                    power = (-exponent if reverse else exponent)+row_heights[i]+column_heights[j]
                    terms[i][j].append((power, coefficient))
    clearing = max(0, -min(e for row in terms for cell in row for e, _ in cell))
    out = [[ZERO]*2 for _ in range(2)]
    for i in range(2):
        for j in range(2):
            low, high = 0, 0
            for exponent, coefficient in terms[i][j]:
                if coefficient & 1:
                    low ^= 1 << (exponent+clearing)
                if coefficient & 2:
                    high ^= 1 << (exponent+clearing)
            out[i][j] = low, high
    return out, clearing


def frame_maps(L, P, cut, low_height):
    H, minus_h = [0, low_height], [0, -1]
    if cut == 'zero':
        left, right = weighted(L, H, minus_h), weighted(P)
        return [left, left, right, right]
    if cut == 'infinity':
        return [weighted(L, reverse=True), weighted(P, H, minus_h, True, True),
                weighted(P, minus_h, H, True), weighted(L, reverse=True, transpose=True)]
    raise ValueError('Unknown place.')


def local_numerator(matrix):
    a, b = matrix[0]
    c, d = matrix[1]
    aa, bb, cc, dd = [multiply(x, x) for x in (a, b, c, d)]
    ab, ac, ad, bc, bd, cd = [multiply(x, y) for x, y in
                             ((a, b), (a, c), (a, d), (b, c), (b, d), (c, d))]
    determinant = add(ad, bc)
    if determinant == ZERO:
        raise ValueError('The one-mode factor must be invertible over the rational function field.')
    return [[dd, cc, cd, cd], [bb, aa, ab, ab],
            [bd, ac, ad, bc], [bd, ac, bc, ad]], determinant


def apply_local(numerator, operator_pair, local, right):
    """Multiply a rational Clifford matrix by a transported rank-four update."""
    n = len(numerator)//2
    columns = []
    for operator in operator_pair:
        columns.append([(v, 0) for v in operator[0]]+[ZERO]
                       +[(v, 0) for v in operator[1]]+[ZERO])
    for position in (n-1, 2*n-1):
        columns.append([ONE if i == position else ZERO for i in range(2*n)])
    embedding = [list(row) for row in zip(*columns)]
    dual = [columns[j][n:]+columns[j][:n] for j in (1, 0, 3, 2)]
    update, determinant = local_numerator(local)
    for i in range(4):
        update[i][i] = add(update[i][i], determinant)
    if right:
        correction = matrix_product(matrix_product(matrix_product(numerator, embedding), update), dual)
    else:
        correction = matrix_product(matrix_product(embedding, update), matrix_product(dual, numerator))
    out = [[add(multiply(determinant, value), correction[i][j]) for j, value in enumerate(row)]
           for i, row in enumerate(numerator)]
    return out, valuation(determinant)


def clifford_matrix(M, L, P, cut, low_height=1):
    operators, diagonal, active = images(M, 1)
    diagonal = diagonal+[0]
    n, base_clearing = len(diagonal), max(diagonal)
    numerator = [[ZERO]*(2*n) for _ in range(2*n)]
    for i, value in enumerate(diagonal):
        numerator[i][i] = 1 << (base_clearing+value), 0
        numerator[n+i][n+i] = 1 << (base_clearing-value), 0
    denominator_order, determinant_order_sum, clearing = base_clearing, 0, 0
    for side, (local, local_clearing) in enumerate(frame_maps(L, P, cut, low_height)):
        numerator, determinant_order = apply_local(numerator, operators[side][0], local, side >= 2)
        denominator_order += determinant_order
        determinant_order_sum += determinant_order
        clearing += local_clearing
    mean_twice = sum(diagonal)+determinant_order_sum-2*clearing
    return numerator, denominator_order, mean_twice, clearing, active


def profile(M, L, P, cut, low_height=1):
    matrix, denominator_order, mean_twice, clearing, active = clifford_matrix(M, L, P, cut, low_height)
    degree = max(max(a.bit_length(), b.bit_length())-1 for row in matrix for a, b in row)
    precision = 128
    while precision <= degree:
        precision *= 2
    attempts = []
    while precision <= 8192:
        try:
            values = MatrixOracle(matrix, 18, precision=precision).profile()['orders']
            attempts.append([precision, True])
            break
        except ArithmeticError:
            attempts.append([precision, False])
            precision *= 2
    else:
        raise ArithmeticError('The complete affine Clifford profile exceeded its precision budget.')
    assert len(values) == 18 and all(value == round(value) for value in values)
    signed = sorted(round(value)-denominator_order for value in values)
    assert signed == [-value for value in signed[::-1]]
    gaps = signed[9:]
    assert gaps[0] == 0
    gaps = gaps[1:]
    assert (mean_twice-sum(gaps)) % 2 == 0
    minimum = (mean_twice-sum(gaps))//2
    orders = [minimum]
    for gap in gaps:
        orders += [value+gap for value in orders]
    assert len(orders) == 256
    return {'cut': cut, 'orders': sorted(orders), 'clearing_power': clearing,
            'mean_order_twice': mean_twice, 'minimum_order': minimum, 'binomial_gaps': gaps,
            'clifford_signed_orders': signed, 'clifford_matrix_order': 18,
            'common_denominator_order': denominator_order, 'numerator_degree': degree,
            'precision_attempts': attempts, 'active_full_mode_indices': active,
            'active_ambient_order': 256, 'active_generic_rank': 256,
            'auxiliary_smith_multiplicity': 2}
