#!/usr/bin/env python3
"""Compile regular 2-adic weight circuits for implicit Cauchy bilinear forms.

Supported gates are integer constants, an odd input, addition, subtraction,
multiplication, binary-encoded integer powers, and inverses of units on the
odd residue class. No circuit node is evaluated at every matrix node.
"""

from pole_cauchy_polynomial_solve import inverse_polynomial_bilinear_form


def compile_regular_odd_circuit(circuit, precision):
    """Return an integer polynomial agreeing with the circuit on every odd x.

    Agreement is modulo 2^precision. The input expansion x=1+2X has
    coefficient i divisible by 2^i; every supported gate preserves this.
    Thus all degrees at least precision vanish when X is any integer.
    """
    if not isinstance(precision, int) or precision < 1 or not circuit:
        raise ValueError('Use positive precision and a nonempty circuit.')
    modulus, limit = 1 << precision, precision
    counts = {'multiplications': 0, 'integer_coefficient_products': 0,
              'unit_inversions': 0, 'binary_power_bits': 0}

    def trim(value):
        while len(value) > 1 and value[-1] == 0:
            value.pop()
        return value

    def add(left, right, sign=1):
        size = max(len(left), len(right))
        result = [0] * size
        for i in range(size):
            result[i] = ((left[i] if i < len(left) else 0)
                         + sign * (right[i] if i < len(right) else 0)) % modulus
        return trim(result)

    def multiply(left, right):
        counts['multiplications'] += 1
        result = [0] * min(limit, len(left) + len(right) - 1)
        for i, a in enumerate(left):
            if not a:
                continue
            for j, b in enumerate(right[:len(result) - i]):
                if b:
                    result[i + j] += a * b
                    counts['integer_coefficient_products'] += 1
        return trim([value % modulus for value in result])

    def invert(value):
        if not value[0] & 1:
            raise ValueError('An inverse gate must have an odd value at input one.')
        counts['unit_inversions'] += 1
        result = [pow(value[0], -1, modulus)]
        for i in range(1, limit):
            coefficient = sum(value[j] * result[i - j]
                              for j in range(1, min(i + 1, len(value))))
            result.append((-result[0] * coefficient) % modulus)
        return trim(result)

    def power(value, exponent):
        if not isinstance(exponent, int):
            raise ValueError('Power exponents must be integers encoded in binary.')
        counts['binary_power_bits'] += abs(exponent).bit_length()
        if exponent < 0:
            value, exponent = invert(value), -exponent
        result = [1]
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
            raise ValueError('Every circuit gate must be a nonempty tuple or list.')
        op = node[0]

        def operand(index):
            if not isinstance(index, int) or not 0 <= index < len(values):
                raise ValueError('Circuit operands must refer to earlier gates.')
            return values[index]

        if op == 'var' and len(node) == 1:
            value = [1, 2 % modulus][:limit]
        elif op == 'const' and len(node) == 2 and isinstance(node[1], int):
            value = [node[1] % modulus]
        elif op in ('add', 'sub', 'mul') and len(node) == 3:
            left, right = operand(node[1]), operand(node[2])
            value = multiply(left, right) if op == 'mul' else add(left, right, 1 if op == 'add' else -1)
        elif op == 'pow' and len(node) == 3:
            value = power(operand(node[1]), node[2])
        elif op == 'inverse' and len(node) == 2:
            value = invert(operand(node[1]))
        else:
            raise ValueError('Unknown gate or invalid gate arguments.')
        value = trim(value)
        if any(coefficient % (1 << i) for i, coefficient in enumerate(value)):
            raise ArithmeticError('The regular odd-input coefficient filtration failed.')
        values.append(value)

    centered = values[-1]
    # c_i X^i = (c_i/2^i)(x-1)^i. A selected lift of c_i suffices:
    # changing it by 2^precision changes the value on odd x by that modulus.
    result, basis = [0] * len(centered), [1]
    for i, coefficient in enumerate(centered):
        weight = coefficient >> i
        for j, entry in enumerate(basis):
            result[j] = (result[j] + weight * entry) % modulus
        basis = [(-basis[0]) % modulus] + [
            (basis[j - 1] - basis[j]) % modulus for j in range(1, len(basis))] + [basis[-1]]
    result = trim(result)
    return {'coefficients': tuple(result), 'precision': precision,
            'degree': len(result) - 1, 'centered_degree': len(centered) - 1,
            'circuit_gates': len(circuit), 'operation_counts': counts,
            'scope': 'Agreement modulo 2^precision on every odd integer.'}


def inverse_analytic_bilinear_form(n, left, right, k, row_circuit, column_circuit,
                                  precision, *, row_denominator_bits=0,
                                  column_denominator_bits=0):
    """Compute a complete inverse bilinear form with compressed regular weights.

    Each domain must be a single progression of positive odd nodes, with equal
    lengths. This gives the proved dimension-independent inverse-entry bound.
    External factors 2^-row_denominator_bits and 2^-column_denominator_bits
    are charged in the returned scale and absolute precision.
    """
    if len(left) != 1 or len(right) != 1 or left[0][1] != right[0][1]:
        raise ValueError('Use two single odd progressions of equal length.')
    for denominator in (row_denominator_bits, column_denominator_bits):
        if not isinstance(denominator, int) or denominator < 0:
            raise ValueError('External dyadic denominator exponents must be nonnegative.')
    row = compile_regular_odd_circuit(row_circuit, precision)
    column = compile_regular_odd_circuit(column_circuit, precision)
    result = inverse_polynomial_bilinear_form(
        n, left, right, k, row['coefficients'], column['coefficients'], precision)
    result['scale_exponent'] -= row_denominator_bits + column_denominator_bits
    result['absolute_precision'] = precision + result['scale_exponent']
    result.update(row_compiler=row, column_compiler=column,
                  row_denominator_bits=row_denominator_bits,
                  column_denominator_bits=column_denominator_bits,
                  enumerated_matrix_nodes=0,
                  scope='Regular analytic inverse moments; the singular factor-count kernel is not evaluated.')
    return result

