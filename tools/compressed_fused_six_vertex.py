#!/usr/bin/env python3
"""Polynomial-size, exactly checked fusion in the repository's F16 convention.

Every rectangular row and every growing output symmetric block is checked
through its two coefficient identities. No generic root-of-unity fusion or
regular endpoint assumption is used. Joint endpoint normalization is explicit.
"""
from __future__ import annotations

from itertools import product

from f16_identity_operator_profiles import MUL, INV, add, multiply, evaluate_one, inverse_matrix
from fused_six_vertex_crossings import power, pscale, pdiv, gcd


def local_transition(first, second, output_first, theta, factor):
    """Return (output_second, polynomial) for one specified first output bit."""
    c = theta ^ INV[theta]
    if first == second:
        if output_first != first:
            return None
        return second, {0: INV[theta], 1: MUL[theta][factor]}
    if output_first == first:
        value = {1: MUL[c][factor]} if first == 0 else {0: c}
        return second, {k: v for k, v in value.items() if v}
    return first, {0: 1, 1: factor}


def rectangular_row(m, theta, shift):
    """Checked bit x Sym_m -> Sym_m x bit map, after scalar gcd removal.

    Keys are (output symmetric weight, output bit, input bit, input weight).
    ``shift`` is the literal fusion row index i; gate j uses theta^(2(i-j))*z.
    The induction uses only Sym_n,k = theta^-k Sym_(n-1),k 0
    + Sym_(n-1),k-1 1, so it works even when quantum binomials vanish.
    """
    if m < 1 or not 0 < theta < 16:
        raise ValueError('Use a positive m and a nonzero F16 parameter.')
    previous = {(0, 0, 0, 0): {0: 1}, (0, 1, 1, 0): {0: 1}}
    checks = 0
    for n in range(1, m + 1):
        result = {}
        for a, b in product(range(2), range(n + 1)):
            expanded = {}
            for input_bit in range(2):
                input_weight = b - input_bit
                if not 0 <= input_weight < n:
                    continue
                scalar = power(theta, -b) if input_bit == 0 else 1
                for carrier in range(2):
                    old_weight = a + input_weight - carrier
                    polynomial = previous.get((old_weight, carrier, a, input_weight), {})
                    if not polynomial:
                        continue
                    for output_bit in range(2):
                        transition = local_transition(carrier, input_bit, output_bit,
                                                      theta, power(theta, 2 * (shift - n + 1)))
                        if transition is None:
                            continue
                        d, value = transition
                        key = old_weight + output_bit, d, output_bit
                        expanded[key] = add(expanded.get(key, {}),
                                            multiply(polynomial, value), scalar)
            for d in range(2):
                c = a + b - d
                if not 0 <= c <= n:
                    continue
                value = expanded.get((c, d, int(c > 0)), {})
                if 0 < c < n:
                    expected = pscale(value, power(theta, -c))
                    if expanded.get((c, d, 0), {}) != expected:
                        raise ArithmeticError(('Rectangular symmetry failed', n, shift, a, b, c, d))
                    checks += 1
                result[c, d, a, b] = value
        previous = result
    common = {}
    for value in previous.values():
        if value:
            common = gcd(common, value)
    if not common:
        raise ArithmeticError('The complete rectangular fusion row vanished.')
    for key, value in previous.items():
        previous[key], remainder = pdiv(value, common)
        if remainder:
            raise ArithmeticError('The rectangular scalar division was not exact.')
    return {'matrix': previous, 'common_factor': common,
            'symmetry_checks': checks,
            'degree': max(max(value, default=0) for value in previous.values())}


def joint_endpoint_normalize(matrix, q):
    """Right multiply by F(1)^-1 within its integer-charge blocks.

    This is a joint operation on the two input legs. It preserves the first
    crossed flattening's Smith profiles but need not preserve the second.
    """
    output = [[{} for _ in range(q * q)] for _ in range(q * q)]
    endpoint_inverse = [[0 for _ in range(q * q)] for _ in range(q * q)]
    for charge in range(2 * q - 1):
        indices = [q * a + b for a, b in product(range(q), repeat=2) if a + b == charge]
        endpoint = [[evaluate_one(matrix[i][j]) for j in indices] for i in indices]
        inverse = inverse_matrix(endpoint)
        for i, row_index in enumerate(indices):
            for j, column_index in enumerate(indices):
                endpoint_inverse[row_index][column_index] = inverse[i][j]
                value = {}
                for k, inner_index in enumerate(indices):
                    if inverse[k][j]:
                        value = add(value, matrix[row_index][inner_index], inverse[k][j])
                output[row_index][column_index] = value
    if not all(evaluate_one(value) == int(i == j)
               for i, row in enumerate(output) for j, value in enumerate(row)):
        raise ArithmeticError('The joint endpoint normalization failed.')
    return output, endpoint_inverse


def construct_compressed(m, theta, normalize_joint=False):
    """Construct and certify the finite fused matrix without full bit spaces."""
    q = m + 1
    rows = [rectangular_row(m, theta, i) for i in range(m)]
    matrix = [[{} for _ in range(q * q)] for _ in range(q * q)]
    checks = 0
    for a, b in product(range(q), repeat=2):
        values = {(a, b, 0): {0: 1}}
        for processed in range(m):
            remaining = m - 1 - processed
            expanded = {}
            for (left, middle, right), polynomial in values.items():
                for input_bit in range(2):
                    new_left = left - input_bit
                    if not 0 <= new_left <= remaining:
                        continue
                    scalar = power(theta, -left) if input_bit == 0 else 1
                    for output_bit in range(2):
                        new_middle = middle + input_bit - output_bit
                        if not 0 <= new_middle <= m:
                            continue
                        coefficient = rows[remaining]['matrix'].get(
                            (new_middle, output_bit, input_bit, middle), {})
                        key = new_left, new_middle, right + output_bit, output_bit
                        expanded[key] = add(expanded.get(key, {}),
                                            multiply(polynomial, coefficient), scalar)
            values = {}
            for left, middle, right in sorted({key[:3] for key in expanded}):
                canonical_bit = int(right > processed)
                value = expanded.get((left, middle, right, canonical_bit), {})
                if 0 < right <= processed:
                    expected = pscale(value, power(theta, -(processed + 1 - right)))
                    if expanded.get((left, middle, right, 1), {}) != expected:
                        raise ArithmeticError(('Output symmetry failed', a, b, processed, left, middle, right))
                    checks += 1
                if value:
                    values[left, middle, right] = value
        for (left, c, d), value in values.items():
            if left:
                raise ArithmeticError('The left input block did not finish.')
            matrix[c * q + d][a * q + b] = value
    common = {}
    for row in matrix:
        for value in row:
            if value:
                common = gcd(common, value)
    if not common:
        raise ArithmeticError('The complete fused matrix vanished.')
    for row in matrix:
        for j, value in enumerate(row):
            row[j], remainder = pdiv(value, common)
            if remainder:
                raise ArithmeticError('The final scalar division was not exact.')
    raw_endpoint = [[evaluate_one(value) for value in row] for row in matrix]
    scalar = raw_endpoint[0][0]
    scalar_identity = bool(scalar) and all(value == (scalar if i == j else 0)
                                          for i, row in enumerate(raw_endpoint)
                                          for j, value in enumerate(row))
    inverse = None
    if scalar_identity:
        matrix = [[pscale(value, INV[scalar]) for value in row] for row in matrix]
    elif normalize_joint:
        matrix, inverse = joint_endpoint_normalize(matrix, q)
    return {'m': m, 'theta': theta, 'alphabet': q, 'invariant': True,
            'rectangular_symmetry_checks': sum(row['symmetry_checks'] for row in rows),
            'output_symmetry_checks': checks,
            'row_degrees': [row['degree'] for row in rows],
            'row_common_degrees': [max(row['common_factor']) for row in rows],
            'common_degree': sum(max(row['common_factor']) for row in rows) + max(common),
            'degree': max(max(value, default=0) for row in matrix for value in row),
            'scalar_identity_endpoint_before_normalization': scalar_identity,
            'joint_endpoint_normalization': bool(inverse),
            'endpoint_identity': scalar_identity or bool(inverse),
            'endpoint_before_normalization': raw_endpoint,
            'endpoint_inverse': inverse,
            'matrix': matrix}
