#!/usr/bin/env python3
"""Implicit complete inverse contractions with regular bivariate circuits."""

from math import comb

from interval_unit_jets import QuadraticJets
from modular_interval_jets import unit_progression_jet
from pole_cauchy_polynomial_solve import _validate, normalized_interpolation_series


def compile_regular_odd_kernel(circuit, precision):
    """Compile F(a,b) modulo 2^precision on all odd pairs.

    Gates match the univariate compiler except that ('var',0) and
    ('var',1) select the two arguments. Total degree is truncated.
    """
    if not isinstance(precision, int) or precision < 1 or not circuit:
        raise ValueError('Use positive precision and a nonempty circuit.')
    modulus, zero = 1 << precision, (0, 0)
    counts = {'multiplications': 0, 'coefficient_products': 0,
              'unit_inversions': 0, 'binary_power_bits': 0}

    def constant(value):
        value %= modulus
        return {zero: value} if value else {}

    def add(left, right, sign=1):
        result = dict(left)
        for key, value in right.items():
            result[key] = (result.get(key, 0) + sign * value) % modulus
            if not result[key]:
                del result[key]
        return result

    def multiply(left, right):
        counts['multiplications'] += 1
        result = {}
        for (i, j), value in left.items():
            for (r, s), other in right.items():
                if i + j + r + s >= precision:
                    continue
                key = i + r, j + s
                result[key] = result.get(key, 0) + value * other
                counts['coefficient_products'] += 1
        return {key: value % modulus for key, value in result.items() if value % modulus}

    def invert(value):
        if not value.get(zero, 0) & 1:
            raise ValueError('Kernel inverse gates must be odd at (1,1).')
        counts['unit_inversions'] += 1
        result, degree = constant(pow(value[zero], -1, modulus)), 1
        while degree < precision:
            result = multiply(result, add(constant(2), multiply(value, result), -1))
            degree *= 2
        return result

    def power(value, exponent):
        if not isinstance(exponent, int):
            raise ValueError('Use integer binary-encoded exponents.')
        counts['binary_power_bits'] += abs(exponent).bit_length()
        if exponent < 0:
            value, exponent = invert(value), -exponent
        result = constant(1)
        while exponent:
            if exponent & 1:
                result = multiply(result, value)
            exponent >>= 1
            if exponent:
                value = multiply(value, value)
        return result

    values = []
    for node in circuit:
        if not isinstance(node, (tuple, list)) or not node:
            raise ValueError('A gate must be a nonempty tuple or list.')

        def operand(index):
            if not isinstance(index, int) or not 0 <= index < len(values):
                raise ValueError('Operands must refer to earlier gates.')
            return values[index]

        op = node[0]
        if op == 'var' and len(node) == 2 and node[1] in (0, 1):
            value = constant(1)
            if precision > 1:
                value[(1, 0) if node[1] == 0 else (0, 1)] = 2
        elif op == 'const' and len(node) == 2 and isinstance(node[1], int):
            value = constant(node[1])
        elif op in ('add', 'sub', 'mul') and len(node) == 3:
            left, right = operand(node[1]), operand(node[2])
            value = multiply(left, right) if op == 'mul' else add(left, right, 1 if op == 'add' else -1)
        elif op == 'inverse' and len(node) == 2:
            value = invert(operand(node[1]))
        elif op == 'pow' and len(node) == 3:
            value = power(operand(node[1]), node[2])
        else:
            raise ValueError('Unknown gate or invalid arguments.')
        if any(coefficient % (1 << sum(key)) for key, coefficient in value.items()):
            raise ArithmeticError('The total-degree coefficient filtration failed.')
        values.append(value)

    polynomial = {}
    for (i, j), coefficient in values[-1].items():
        weight = coefficient >> (i + j)
        for d in range(i + 1):
            for e in range(j + 1):
                value = weight * comb(i, d) * comb(j, e)
                if (i + j - d - e) & 1:
                    value = -value
                key = d, e
                polynomial[key] = (polynomial.get(key, 0) + value) % modulus
    terms = [(d, e, value) for (d, e), value in sorted(polynomial.items()) if value]
    return {'terms': terms, 'precision': precision, 'circuit_gates': len(circuit),
            'total_degree': max((d + e for d, e, _ in terms), default=0),
            'centered_nonzero_terms': len(values[-1]), 'operation_counts': counts,
            'scope': 'Pointwise agreement on every odd pair; no singular denominator inversion.'}


def inverse_polynomial_kernel_contraction(n, left, right, k, terms, precision):
    """Compute sum_ab 2^k (M^-1)[b,a] K(a,b) for an explicit polynomial K.

    Degrees, rather than only their binary lengths, are charged. The
    diagonal-prefix recurrence evaluates all requested moments in O(D E)
    scalar ring operations after the two interpolation jets are known.
    """
    coefficients = {}
    for d, e, value in terms:
        if any(not isinstance(v, int) for v in (d, e, value)) or min(d, e) < 0:
            raise ValueError('Terms must be nonnegative integer degrees and integer coefficients.')
        coefficients[d, e] = coefficients.get((d, e), 0) + value
    row_degree = max((d for d, _ in coefficients), default=0)
    column_degree = max((e for _, e in coefficients), default=0)
    size = _validate(n, left, right, right[0][0], k, precision)
    origin = normalized_interpolation_series(n, left, right, k, precision, row_degree)
    ring = QuadraticJets(precision, column_degree)
    y = tuple([(0, 0), (1, 0)][:ring.size] + [(0, 0)] * max(0, ring.size - 2))
    z = ring.constant(n, 1 << k)
    zy = ring.multiply(z, y)
    numerator, denominator = ring.constant(1), ring.constant(1)
    for first, length in left:
        factor = unit_progression_jet(ring, (0, length),
                                      ring.add(ring.constant(first), ring.scale(zy, -1)), ring.constant(2))
        factor = ring.multiply(factor, ring.inverse(ring.constant(*factor[0])))
        numerator = ring.multiply(numerator, factor)
    for first, length in right:
        factor = unit_progression_jet(ring, (0, length),
                                      ring.add(ring.constant(1), ring.scale(y, -first)), ring.scale(y, -2))
        denominator = ring.multiply(denominator, factor)
    infinity = ring.multiply(numerator, ring.inverse(denominator))
    scalar = QuadraticJets(precision, 0)
    z = scalar.constant(n, 1 << k)
    origin_constant = scalar.constant(1)
    for domain in (left, right):
        for first, length in domain:
            origin_constant = scalar.multiply(origin_constant, unit_progression_jet(
                scalar, (0, length), scalar.constant(first), scalar.constant(2)))
    origin_constant = scalar.multiply(origin_constant, scalar.inverse(scalar.power(z, size)))
    origin = tuple(scalar.pair_product(origin_constant[0], value) for value in origin)
    result, factor = scalar.constant(0), scalar.multiply(z, scalar.inverse(scalar.constant(0, 1)))
    previous = [(0, 0)] * (column_degree + 1)
    for d in range(row_degree + 1):
        current = []
        for e in range(column_degree + 1):
            convolution = scalar.pair_product(origin[d], infinity[e])
            if d and e:
                convolution = ((convolution[0] + previous[e - 1][0]) % scalar.modulus,
                               (convolution[1] + previous[e - 1][1]) % scalar.modulus)
            current.append(convolution)
            coefficient = coefficients.get((d, e), 0)
            if coefficient:
                moment = scalar.constant(int(d == e) - convolution[0], -convolution[1])
                result = scalar.add(result, scalar.scale(scalar.multiply(factor, moment), coefficient))
        previous = current
        factor = scalar.multiply(factor, z)
    return {'residue': list(result[0]), 'precision': precision, 'scale_exponent': -k,
            'row_degree': row_degree, 'column_degree': column_degree,
            'moment_recurrence_cells': (row_degree + 1) * (column_degree + 1)}


def inverse_analytic_kernel_contraction(n, left, right, k, circuit, precision,
                                        *, denominator_bits=0):
    if len(left) != 1 or len(right) != 1 or left[0][1] != right[0][1]:
        raise ValueError('Use two single odd progressions of equal length.')
    if not isinstance(denominator_bits, int) or denominator_bits < 0:
        raise ValueError('The external dyadic denominator exponent must be nonnegative.')
    compiled = compile_regular_odd_kernel(circuit, precision)
    result = inverse_polynomial_kernel_contraction(n, left, right, k, compiled['terms'], precision)
    result['scale_exponent'] -= denominator_bits
    result.update(absolute_precision=precision + result['scale_exponent'],
                  compiler=compiled, denominator_bits=denominator_bits,
                  enumerated_matrix_nodes=0,
                  scope='A complete regular bivariate-kernel inverse contraction, with all precision losses charged.')
    return result

