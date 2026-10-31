#!/usr/bin/env python3
"""Exact secondary-height repairs of the four-by-three core over F4."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement, permutations, product
import json
from pathlib import Path
import random

from check_rectangular_f4_repair import (
    HX, HU, ZERO, add, scale, times_binary, multiply, binary_power,
    divide_binary, evaluate_one, c_order, ramify_binary, matrix_multiply,
    partial_repair, build, fmul, finv,
)


def height(index):
    return HX[index % 4]+HU[index // 4]


def validate_postcorrection(matrix, secondary_x, secondary_y, keep_background=True):
    if len(matrix) != 12 or any(len(row) != 12 for row in matrix):
        raise ValueError('A postcorrection must be a twelve-by-twelve matrix.')
    if len(secondary_x) != 4 or len(secondary_y) != 3:
        raise ValueError('The secondary height has four X values and three Y values.')
    keys = [(height(i), secondary_x[i % 4]+secondary_y[i // 4]) for i in range(12)]
    for i, row in enumerate(matrix):
        for j, coefficient in enumerate(row):
            if coefficient not in (0, 1, 2, 3):
                raise ValueError('Postcorrection coefficients must lie in F4.')
            if i == j:
                if coefficient != 1:
                    raise ValueError('The postcorrection diagonal must be one.')
            elif coefficient and not keys[i] < keys[j]:
                raise ValueError('A nonzero off-diagonal entry must strictly lower the lexicographic additive height.')
            if keep_background and i == 0 and j != 0 and coefficient:
                raise ValueError('This family fixes the Gaussian background row.')
    return keys


def single_shear(row, column, coefficient=1):
    if row == column or not 0 <= row < 12 or not 0 <= column < 12 or coefficient not in (1, 2, 3):
        raise ValueError('Supply distinct core indices and a nonzero F4 coefficient.')
    if height(row) > height(column):
        raise ValueError('The actual height may not increase.')
    wx, wy = [0]*4, [0]*3
    if height(row) == height(column):
        if row % 4 != column % 4:
            wx[row % 4], wx[column % 4] = -1, 1
        else:
            wy[row // 4], wy[column // 4] = -1, 1
    matrix = [[int(i == j) for j in range(12)] for i in range(12)]
    matrix[row][column] = coefficient
    validate_postcorrection(matrix, wx, wy, keep_background=row != 0)
    return matrix, wx, wy


def refined_repair(alpha, beta, postcorrection, secondary_x, secondary_y, keep_background=True):
    validate_postcorrection(postcorrection, secondary_x, secondary_y, keep_background)
    return matrix_multiply(postcorrection, partial_repair(alpha, beta))


def conjugate_postcorrection(postcorrection, multiplier):
    if multiplier not in (1, 2, 3):
        raise ValueError('A parameter gauge needs a nonzero F4 multiplier.')
    charges = [(1, 0, 1, 2)[i % 4]+(1, 0, 1)[i // 4] for i in range(12)]
    powers = (1, multiplier, fmul(multiplier, multiplier))
    return [[fmul(value, powers[(charges[i]-charges[j]) % 3]) for j, value in enumerate(row)]
            for i, row in enumerate(postcorrection)]


def postcorrect_seed(seed, postcorrection, secondary_x, secondary_y, check_types=True, *, keep_background=True):
    keys = validate_postcorrection(postcorrection, secondary_x, secondary_y, keep_background)
    if 'input_correction' in seed:
        validate_postcorrection(seed['input_correction'], secondary_x, secondary_y, keep_background=False)
    A, B, bulk = seed['q_x'], seed['q_u'], seed['bulk_size']
    def position(x, y, u, v):
        return x+A*(y+B*(u+B*v))
    gate = [[(source % 4, source // 4, coefficient,
              binary_power(3, height(source)-height(target)))
             for source, coefficient in enumerate(row) if coefficient]
            for target, row in enumerate(postcorrection)]
    values = [ZERO]*len(seed['values'])
    for x, y, u, v in product(range(A), range(B), range(B), range(A)):
        answer = ZERO
        for vs, ys, coefficient, polynomial in gate[v//bulk+4*(y//bulk)]:
            source = seed['values'][position(x, ys*bulk+y % bulk, u, vs*bulk+v % bulk)]
            answer = add(answer, scale(times_binary(source, polynomial), coefficient))
        values[position(x, y, u, v)] = answer
    endpoint = matrix_multiply(postcorrection, seed['ramified_endpoint_matrix'])
    result = _verified_repair(seed, values, endpoint, keys, check_types)
    total = (matrix_multiply(postcorrection, seed['postcorrection'])
             if 'postcorrection' in seed else postcorrection)
    result.update({'repair_kind': 'mixed_refined_height_partial',
                   'postcorrection': total, 'secondary_x': secondary_x, 'secondary_y': secondary_y,
                   'repair_matrix': matrix_multiply(postcorrection, seed['repair_matrix'])})
    if not keep_background or 'output_background_fixed' in seed:
        result['output_background_fixed'] = not any(total[0][j] for j in range(1, 12))
    return result


def _verified_repair(seed, values, endpoint, keys, check_types):
    A, B, bulk = seed['q_x'], seed['q_u'], seed['bulk_size']
    def position(x, y, u, v):
        return x+A*(y+B*(u+B*v))
    degree = seed['denominator_degree']
    cancellation = min(degree, min(c_order(value) for value in values if value != ZERO))
    if cancellation:
        values = [divide_binary(value, binary_power(3, cancellation)) for value in values]
        degree -= cancellation
    for i, row in enumerate(endpoint):
        assert row[i] == 1
        assert all(not coefficient or i == j or keys[i] < keys[j] for j, coefficient in enumerate(row))
    hx, hu, top = seed['logical_heights_x_v'], seed['logical_heights_u_y'], seed['indices'][-1]
    for x, y, u, v in product(range(A), range(B), range(B), range(A)):
        polynomial = tuple(ramify_binary(plane, top) for plane in values[position(x, y, u, v)])
        shift = hu[y]+hx[v]-hx[x]-hu[u]-degree*top
        if shift < 0:
            assert not ((polynomial[0] | polynomial[1]) & ((1 << -shift)-1))
            polynomial = tuple(plane >> -shift for plane in polynomial)
        else:
            polynomial = tuple(plane << shift for plane in polynomial)
        expected = endpoint[v//bulk+4*(y//bulk)][x//bulk+4*(u//bulk)]
        expected *= int(x % bulk == v % bulk and y % bulk == u % bulk)
        assert (polynomial[0] & 1) | ((polynomial[1] & 1) << 1) == expected
    checks = 0
    if check_types:
        denominator = binary_power(3, 2*degree)
        xtypes = [sorted(set(permutations(pair))) for pair in combinations_with_replacement(range(A), 2)]
        utypes = [sorted(set(permutations(pair))) for pair in combinations_with_replacement(range(B), 2)]
        for X, U in product(xtypes, utypes):
            for x, v, y, u in product(X, X, U, U):
                polynomial = multiply(values[position(x[0], y[0], u[0], v[0])], values[position(x[1], y[1], u[1], v[1])])
                polynomial = divide_binary(polynomial, denominator)
                assert evaluate_one(polynomial) == int(x == v and y == u)
                checks += 1
    result = dict(seed)
    result.update({'ramified_endpoint_matrix': endpoint, 'values': values,
                   'denominator_degree': degree, 'degree': max(max(value).bit_length()-1 for value in values),
                   'cancelled_c_order': 4-degree,
                   'ramified_endpoint_checks': len(values), 'two_copy_type_checks': checks,
                   'independent_exterior_checks': 0,
                   'nonzero_entries': sum(value != ZERO for value in values),
                   'c_order_histogram': dict(sorted(Counter(c_order(value) for value in values if value != ZERO).items())),
                   'coefficient_sha256': sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()})
    return result


def validate_precorrection(matrix, secondary_x, secondary_y, keep_background=True):
    keys = validate_postcorrection(matrix, secondary_x, secondary_y, keep_background=False)
    if keep_background and any(matrix[row][0] for row in range(1, 12)):
        raise ValueError('This input family fixes the Gaussian background column.')
    return keys


def precorrect_seed(seed, precorrection, secondary_x, secondary_y, check_types=True, *, keep_background=True):
    """Apply the dressed input map H^-1 W H, using one flag for V and W."""
    keys = validate_precorrection(precorrection, secondary_x, secondary_y, keep_background)
    if 'postcorrection' in seed:
        validate_postcorrection(seed['postcorrection'], secondary_x, secondary_y, keep_background=False)
    if 'input_correction' in seed:
        validate_precorrection(seed['input_correction'], secondary_x, secondary_y, keep_background=False)
    A, B, bulk = seed['q_x'], seed['q_u'], seed['bulk_size']
    def position(x, y, u, v):
        return x+A*(y+B*(u+B*v))
    gate = [[(source % 4, source // 4, precorrection[source][target],
              binary_power(3, height(target)-height(source)))
             for source in range(12) if precorrection[source][target]] for target in range(12)]
    values = [ZERO]*len(seed['values'])
    for x, y, u, v in product(range(A), range(B), range(B), range(A)):
        answer = ZERO
        for xs, us, coefficient, polynomial in gate[x//bulk+4*(u//bulk)]:
            source = seed['values'][position(xs*bulk+x % bulk, y, us*bulk+u % bulk, v)]
            answer = add(answer, scale(times_binary(source, polynomial), coefficient))
        values[position(x, y, u, v)] = answer
    endpoint = matrix_multiply(seed['ramified_endpoint_matrix'], precorrection)
    result = _verified_repair(seed, values, endpoint, keys, check_types)
    total = (matrix_multiply(seed['input_correction'], precorrection)
             if 'input_correction' in seed else precorrection)
    result.update({'repair_kind': 'mixed_refined_two_sided_partial',
                   'input_correction': total, 'secondary_x': secondary_x, 'secondary_y': secondary_y})
    if not keep_background or 'input_background_fixed' in seed:
        result['input_background_fixed'] = not any(total[i][0] for i in range(1, 12))
    return result


def parameter_controls():
    terms = [(1, 0), (2, 1), (2, 2)]
    bases = {k: build((0, 1, 2), 4, k, fmul(2, finv(k)), check_types=False, quartic_terms=terms)
             for k in (1, 2, 3)}
    conjugate_base = build((0, 1, 2), 4, 1, 3, check_types=False,
                           quartic_terms=[(fmul(a, a), e) for a, e in terms])
    A, B, bulk = 8, 6, 2
    qx = [(letter % bulk).bit_count()+(1, 0, 1, 2)[letter//bulk] for letter in range(A)]
    qu = [(letter % bulk).bit_count()+(1, 0, 1)[letter//bulk] for letter in range(B)]
    scaling, conjugation = 0, 0
    for row, column in ((3, 8), (8, 3), (10, 11), (1, 11)):
        matrix, wx, wy = single_shear(row, column, 2)
        base = postcorrect_seed(bases[1], matrix, wx, wy, check_types=False)
        for k in (1, 2, 3):
            transformed = postcorrect_seed(bases[k], conjugate_postcorrection(matrix, k), wx, wy, check_types=False)
            assert transformed['denominator_degree'] == base['denominator_degree']
            powers = (1, k, fmul(k, k))
            for x, y, u, v in product(range(A), range(B), range(B), range(A)):
                position = x+A*(y+B*(u+B*v))
                multiplier = powers[(qx[v]+qu[y]-qx[x]-qu[u]) % 3]
                assert transformed['values'][position] == scale(base['values'][position], multiplier)
                scaling += 1
        conjugate = postcorrect_seed(conjugate_base, [[fmul(a, a) for a in line] for line in matrix], wx, wy, check_types=False)
        assert conjugate['denominator_degree'] == base['denominator_degree']
        for original, transformed in zip(base['values'], conjugate['values']):
            assert transformed == (original[0] ^ original[1], original[1])
            conjugation += 1
    return {'quartic_terms': terms, 'scaling_entry_checks': scaling, 'frobenius_entry_checks': conjugation,
            'postcorrection_rule': 'V_rs -> k^(q(r)-q(s)) V_rs'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    base = build((0, 1), 2, 1, 2)
    summaries, counts = [], {'strict_height': 0, 'equal_height': 0}
    for row in range(1, 12):
        for column in range(12):
            if row == column or height(row) > height(column):
                continue
            counts['strict_height' if height(row) < height(column) else 'equal_height'] += 1
            for coefficient in (1, 2, 3):
                matrix, wx, wy = single_shear(row, column, coefficient)
                seed = postcorrect_seed(base, matrix, wx, wy)
                summaries.append({'row': row, 'column': column, 'coefficient': coefficient,
                                  'denominator_degree': seed['denominator_degree'],
                                  'ramified_endpoint_checks': seed['ramified_endpoint_checks'],
                                  'two_copy_type_checks': seed['two_copy_type_checks'],
                                  'coefficient_sha256': seed['coefficient_sha256']})
    rng = random.Random(77131)
    background = []
    for column in range(1, 12):
        for coefficient in (1, 2, 3):
            matrix, wx, wy = single_shear(0, column, coefficient)
            seed = postcorrect_seed(base, matrix, wx, wy, keep_background=False)
            assert not seed['output_background_fixed']
            background.append({'row': 0, 'column': column, 'coefficient': coefficient,
                               'denominator_degree': seed['denominator_degree'],
                               'ramified_endpoint_checks': seed['ramified_endpoint_checks'],
                               'two_copy_type_checks': seed['two_copy_type_checks'],
                               'coefficient_sha256': seed['coefficient_sha256']})
    dense = []
    for case in range(3):
        wx, wy = [rng.randrange(-50, 51) for _ in range(4)], [rng.randrange(-50, 51) for _ in range(3)]
        keys = [(height(i), wx[i % 4]+wy[i // 4]) for i in range(12)]
        matrix = [[int(i == j) if i == j else rng.randrange(4) if i and keys[i] < keys[j] else 0
                   for j in range(12)] for i in range(12)]
        seed = postcorrect_seed(build((0, 1, 2), 4, 1, 2, check_types=False), matrix, wx, wy)
        dense.append({key: value for key, value in seed.items() if key not in ('values', 'residual_matrix', 'repair_matrix', 'ramified_endpoint_matrix')})
    result = {'status': 'exact secondary-height polynomial repairs; no pressure claim',
              'single_shear_directions_with_row_zero_fixed': counts,
              'single_shear_cases': summaries, 'dense_flag_cases': dense,
              'background_changing_output_cases': background,
              'parameter_gauges': parameter_controls()}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
