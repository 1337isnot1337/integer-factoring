#!/usr/bin/env python3
"""Generic degree-six Smith profiles over F4(y), with exact row certificates.

An entry is a pair of nonnegative integers encoding the binary and alpha
parts of a bivariate polynomial.  Bit x+256*y represents x^x y^y;
alpha^2=alpha+1.  ``profile_matrix`` accepts square matrices through order
16.  ``profile_tensor`` accepts 256 homogeneous degree-six polynomials,
as dictionaries from exponent triples to coefficients 0,1,2,3, or as
binary monomial sets.  No physical parity restriction is imposed.

The native engine uses fraction-free elimination, never specialization.
The Python certificate verifier independently multiplies F*A, divides
each row by its claimed x order, and checks det(F(0,y)) and det(H(0,y)).
These two nonzero determinants prove the complete Smith profile over
F4(y)[[x]].  Input degree six bounds minors by degree 96 and Bareiss
products by 192, so stride 256 is sufficient throughout.
"""
from __future__ import annotations

import hashlib
import json
import subprocess

from cube_bivariate_smith import native_program
from cube_homogeneous_gaussian import axis_matrix

STRIDE = 256
ZERO = (0, 0)
ONE = (1, 0)
F4_MUL = ((0,0,0,0), (0,1,2,3), (0,2,3,1), (0,3,1,2))
F4_INV = (0,1,3,2)


def binary_multiply(a, b):
    if a.bit_count() < b.bit_count():
        a, b = b, a
    result = 0
    while b:
        bit = b & -b
        b ^= bit
        result ^= a << (bit.bit_length()-1)
    return result


def add(a, b):
    return a[0] ^ b[0], a[1] ^ b[1]


def multiply(a, b):
    # Four binary cross products, independent of native Karatsuba packing.
    aa = binary_multiply(a[0], b[0])
    bb = binary_multiply(a[1], b[1])
    ab = binary_multiply(a[0], b[1])
    ba = binary_multiply(a[1], b[0])
    return aa ^ bb, ab ^ ba ^ bb


def scale(a, coefficient):
    if coefficient == 0:
        return ZERO
    if coefficient == 1:
        return a
    if coefficient == 2:
        return a[1], a[0] ^ a[1]
    if coefficient == 3:
        return a[0] ^ a[1], a[0]
    raise ValueError('An F4 coefficient must be in 0,1,2,3.')


def degree(a):
    return max(a[0].bit_length(), a[1].bit_length()) - 1


def coefficient(a, at):
    return ((a[0] >> at) & 1) | (((a[1] >> at) & 1) << 1)


def divide_exact(a, b):
    if b == ZERO:
        raise ZeroDivisionError('Zero F4 polynomial divisor.')
    db = degree(b)
    inverse = F4_INV[coefficient(b, db)]
    multiples = [scale(b, c) for c in range(4)]
    result = [0, 0]
    while degree(a) >= db:
        at = degree(a)
        shift = at-db
        c = F4_MUL[coefficient(a, at)][inverse]
        if c & 1:
            result[0] ^= 1 << shift
        if c & 2:
            result[1] ^= 1 << shift
        term = multiples[c]
        a = a[0] ^ (term[0] << shift), a[1] ^ (term[1] << shift)
    if a != ZERO:
        raise ArithmeticError('F4 polynomial division was not exact.')
    return tuple(result)


def determinant(matrix):
    """Exact polynomial Bareiss determinant, valid for either variable packing."""
    work = [list(row) for row in matrix]
    n = len(work)
    if not n:
        return ONE
    previous = ONE
    for k in range(n-1):
        at = next((i for i in range(k,n) if work[i][k] != ZERO), None)
        if at is None:
            return ZERO
        work[k], work[at] = work[at], work[k]
        pivot = work[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator = add(multiply(pivot,work[i][j]),
                                multiply(work[i][k],work[k][j]))
                work[i][j] = divide_exact(numerator,previous)
            work[i][k] = ZERO
        previous = pivot
    return work[-1][-1]


def terms(value):
    bits = value[0] | value[1]
    while bits:
        bit = bits & -bits
        bits ^= bit
        at = bit.bit_length()-1
        yield at % STRIDE, at // STRIDE, coefficient(value, at)


def encode(polynomial, axis):
    if axis not in range(3):
        raise ValueError('Unknown cube axis.')
    other = next(i for i in range(3) if i != axis)
    items = polynomial.items() if hasattr(polynomial,'items') else ((e,1) for e in polynomial)
    result = [0,0]
    for exponent,c in items:
        if len(exponent) != 3 or min(exponent) < 0 or sum(exponent) != 6:
            raise ValueError('Every tensor monomial must have total degree six.')
        if c not in range(4):
            raise ValueError('An F4 coefficient must be in 0,1,2,3.')
        bit = 1 << (exponent[axis]+STRIDE*exponent[other])
        if c & 1:
            result[0] ^= bit
        if c & 2:
            result[1] ^= bit
    return tuple(result)


def validate_tensor(tensor, endpoint=True):
    if len(tensor) != 256:
        raise ValueError('Expected the complete eight-leg binary-physical tensor.')
    for mask,polynomial in enumerate(tensor):
        encode(polynomial,0)
        if endpoint:
            value = 0
            coefficients = polynomial.values() if hasattr(polynomial,'values') else (1 for _ in polynomial)
            for c in coefficients:
                value ^= c
            wanted = int(all((mask>>i&1)==(mask>>(i^7)&1) for i in range(4)))
            if value != wanted:
                raise ValueError('The literal antipodal crossing endpoint failed.')


def validate_matrix(matrix):
    n = len(matrix)
    if not 1 <= n <= 16 or any(len(row) != n for row in matrix):
        raise ValueError('Expected a square matrix of order one through sixteen.')
    result = []
    for row in matrix:
        values = []
        for value in row:
            if isinstance(value,int):
                value = value,0
            if len(value) != 2 or any(not isinstance(v,int) or v < 0 for v in value):
                raise ValueError('Expected two nonnegative binary polynomial encodings.')
            value = tuple(value)
            if any(x+y > 6 for x,y,c in terms(value)):
                raise ValueError('Matrix entry total degree exceeds six.')
            values.append(value)
        result.append(values)
    return result


def profile_matrix(matrix):
    matrix = validate_matrix(matrix)
    n = len(matrix)
    data = str(n)+'\n'+'\n'.join(f'{a:x} {b:x}' for row in matrix for a,b in row)+'\n'
    result = subprocess.run([str(native_program('cube_f4_bivariate_smith'))],
                            input=data,text=True,capture_output=True)
    if result.returncode:
        raise ArithmeticError(result.stderr.strip())
    record = json.loads(result.stdout)
    record['frames'] = [[tuple(int(value,16) for value in pair) for pair in row]
                        for row in record['frames']]
    return record


def tensor_matrix(tensor, axis, endpoint=True):
    validate_tensor(tensor,endpoint=endpoint)
    matrix = axis_matrix(tensor,axis)
    return [[encode(value,axis) for value in row] for row in matrix]


def profile_tensor(tensor, axis, endpoint=True):
    record = profile_matrix(tensor_matrix(tensor,axis,endpoint=endpoint))
    record['axis'] = axis
    record['labels'] = list(range(16))
    return record


def x_coefficient(value, power):
    result = [0,0]
    for x,y,c in terms(value):
        if x == power:
            if c & 1:
                result[0] ^= 1 << y
            if c & 2:
                result[1] ^= 1 << y
    return tuple(result)


def verify_matrix(matrix, record):
    """Certify F*A=diag(x**orders)*H and both frames regular invertible."""
    matrix = validate_matrix(matrix)
    n = len(matrix)
    if record.get('field') != 'F4' or record.get('stride') != STRIDE:
        raise AssertionError('Incorrect certificate field or monomial stride.')
    orders = record['orders']
    frame = [[tuple(value) for value in row] for row in record['frames']]
    assert len(orders) == n and orders == sorted(orders) and min(orders) >= 0
    assert len(frame) == n and all(len(row) == n for row in frame)
    assert sorted(record['rows']) == sorted(record['columns']) == list(range(n))
    assert record['prefixes'] == [sum(orders[:i+1]) for i in range(n)]
    assert all(x+y <= 6*(n-1) for row in frame for value in row for x,y,c in terms(value))
    product = [[ZERO for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                product[i][j] = add(product[i][j],multiply(frame[i][k],matrix[k][j]))
            assert all(x >= orders[i] for x,y,c in terms(product[i][j]))
            assert all(x+y <= 6*n for x,y,c in terms(product[i][j]))
    constant = [[x_coefficient(value,0) for value in row] for row in frame]
    quotient = [[x_coefficient(value,orders[i]) for value in row]
                for i,row in enumerate(product)]
    first,second = determinant(constant),determinant(quotient)
    assert first != ZERO and second != ZERO
    # Shifting is exact: no x exponent is below the row's claimed order.
    right = [[(value[0]>>orders[i],value[1]>>orders[i]) for value in row]
             for i,row in enumerate(product)]
    payload = json.dumps([[[format(a,'x'),format(b,'x')] for a,b in row]
                          for row in right],separators=(',',':')).encode()
    return dict(polynomial_product_entries=n*n,complete_orders=n,
                unit_determinants=[[format(a,'x'),format(b,'x')] for a,b in (first,second)],
                right_frame_sha256=hashlib.sha256(payload).hexdigest())


def verify_tensor(tensor, axis, record, endpoint=True):
    assert record.get('axis') == axis and record.get('labels') == list(range(16))
    return verify_matrix(tensor_matrix(tensor,axis,endpoint=endpoint),record)


profile = profile_tensor
