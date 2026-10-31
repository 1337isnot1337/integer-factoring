#!/usr/bin/env python3
"""Exact controls for the affine family of polynomial frozen-core crossings."""
from __future__ import annotations

import argparse
from itertools import combinations_with_replacement, permutations, product
import json
from pathlib import Path
import random

from check_rectangular_f4_repair import (
    ZERO, add, multiply, binary_power, divide_binary, evaluate_one, fmul, scale,
    ramify_binary, full_inverse_repair, residual, build,
)
from rectangular_frozen_core_rows import row_components, combine_row


def ramified_control(seed):
    bulk, top = seed['bulk_size'], seed['indices'][-1]
    heights = seed['logical_heights']
    for x, y, u, v in product(range(bulk), repeat=4):
        value = seed['values'][x+bulk*(y+bulk*(u+bulk*v))]
        ramified = tuple(ramify_binary(plane, top) for plane in value)
        shift = heights[y]+heights[v]-heights[x]-heights[u]
        if shift < 0:
            assert not ((ramified[0] | ramified[1]) & ((1 << -shift)-1))
            ramified = tuple(plane >> -shift for plane in ramified)
        else:
            ramified = tuple(plane << shift for plane in ramified)
        actual = (ramified[0] & 1) | ((ramified[1] & 1) << 1)
        expected = seed['ramified_endpoint_scalar'] if x == v and y == u else 0
        assert actual == expected
    return bulk**4


def type_control(seed):
    bulk = seed['bulk_size']
    types = [sorted(set(permutations(pair)))
             for pair in combinations_with_replacement(range(bulk), 2)]
    checks = 0
    for xtype, utype in product(types, repeat=2):
        for x, v, y, u in product(xtype, xtype, utype, utype):
            a = seed['values'][x[0]+bulk*(y[0]+bulk*(u[0]+bulk*v[0]))]
            b = seed['values'][x[1]+bulk*(y[1]+bulk*(u[1]+bulk*v[1]))]
            assert evaluate_one(multiply(a, b)) == int(x == v and y == u)
            checks += 1
    return checks


def original_tensor_controls():
    data = {beta: row_components(beta=beta) for beta in (1, 2)}
    controls = []
    for alpha, beta, repair in ((1, 2, 'partial'), (1, 1, 'inverse'), (1, 2, 'inverse')):
        g = ([1]+[0]*11) if repair == 'partial' else full_inverse_repair(alpha, beta)[0]
        frozen = combine_row(data[beta], g)
        original = build(tuple(range(4)), 4, alpha, beta, check_types=False, repair=repair)
        A, B, bulk = original['q_x'], original['q_u'], original['bulk_size']
        denominator = binary_power(3, original['denominator_degree'])
        for x, y, u, v in product(range(bulk), repeat=4):
            expected = divide_binary(original['values'][x+A*(y+B*(u+B*v))], denominator)
            assert frozen['values'][x+bulk*(y+bulk*(u+bulk*v))] == expected
        controls.append({'alpha': alpha, 'beta': beta, 'repair': repair,
                         'row_coefficients': g, 'polynomial_entries_checked': bulk**4})
    return controls


def cubic_control():
    g = [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0]
    values = combine_row(row_components(beta=1), g)['values']
    assert values[0] == (4, 0)  # w^2
    # Bits 0,1,2 denote x_1,x_2,y_1 in the direct bulk coefficient order.
    residual_value = multiply(values[0], values[7])
    for pair, single in ((3, 4), (5, 2), (6, 1)):
        residual_value = add(residual_value, multiply(values[pair], values[single]))
    expected = (binary_power(3, 14) << 6, 0)
    assert residual_value == expected and residual_value != ZERO
    gaussian = combine_row(row_components(beta=2), [1]+[0]*11)['values']
    gaussian_residual = multiply(gaussian[0], gaussian[7])
    for pair, single in ((3, 4), (5, 2), (6, 1)):
        gaussian_residual = add(gaussian_residual, multiply(gaussian[pair], gaussian[single]))
    assert gaussian_residual == ZERO
    return {'M': 4, 'beta': 1, 'row_coefficients': g,
            'vacuum_coefficient': 'w^2', 'cubic_variables': ['x_1', 'x_2', 'y_1'],
            'cubic_relation_residual': 'w^6*(1+w)^14',
            'exact_residual_bitplanes': list(residual_value),
            'gaussian_control_residual_bitplanes': list(gaussian_residual),
            'conclusion': 'outside the scalar times product(1+l_i*z_i) product(1+q_ij*z_i*z_j) model'}


def parameter_orbit_controls():
    rng = random.Random(310829)
    checks = 0
    for indices, ambient in (((0, 1, 2), 4), ((0, 1, 2, 3), 4), ((0, 2, 5), 8)):
        data = {beta: row_components(indices, ambient, beta) for beta in (1, 2, 3)}
        bulk = data[1]['bulk_size']
        for gamma in (2, 3):
            powers = (1, gamma, fmul(gamma, gamma))
            for _ in range(3):
                g = [rng.randrange(4) for _ in range(12)]
                g[10] = 1 ^ g[0]
                changed = [fmul(c, powers[(2-(1, 0, 1, 2)[r % 4]-(1, 0, 1)[r//4]) % 3])
                           for r, c in enumerate(g)]
                actual = combine_row(data[gamma], g)['values']
                transformed = combine_row(data[1], changed)['values']
                for x, y, u, v in product(range(bulk), repeat=4):
                    position = x+bulk*(y+bulk*(u+bulk*v))
                    delta = x.bit_count()+u.bit_count()-v.bit_count()-y.bit_count()
                    assert actual[position] == scale(transformed[position], powers[delta % 3])
                    checks += 1
    return {'coefficient_checks': checks, 'nonzero_beta_representative': 1,
            'row_map': "g'_r=gamma^(2-q_r)*g_r",
            'external_map': 'gamma^(|x|+|u|-|v|-|y|)'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    residual_checks = 0
    for alpha, beta in product(range(4), repeat=2):
        column = [row[0] for row in residual(alpha, beta)]
        assert column == [int(i in (0, 10)) for i in range(12)]
        residual_checks += 12

    rng = random.Random(310826)
    cases = []
    for indices, ambient in (((0, 1), 2), ((0, 1, 2), 4),
                             ((0, 1, 2, 3), 4), ((0, 2, 5), 8)):
        for beta in (0, 1, 2):
            data = row_components(indices, ambient, beta)
            if beta == 1:
                inverse_row = [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0]
                vacuum = combine_row(data, inverse_row)['values'][0]
                assert vacuum == (1 ^ (5 << (2*(ambient-1-indices[-1]))), 0)
            basis_checks = 0
            for source in range(12):
                g = [0]*12
                g[source] = 1
                basis_checks += ramified_control(combine_row(data, g, require_crossing=False))
            typed_checks = 0
            for _ in range(3):
                g = [rng.randrange(4) for _ in range(12)]
                g[10] = 1 ^ g[0]
                seed = combine_row(data, g)
                ramified_control(seed)
                typed_checks += type_control(seed)
            cases.append({'indices': list(indices), 'ambient': ambient, 'beta': beta,
                          'basis_ramified_entries_checked': basis_checks,
                          'independent_two_copy_type_checks': typed_checks})

    try:
        combine_row(row_components((0, 1), 2), [0]*12)
    except ValueError:
        rejected = True
    else:
        raise AssertionError('A zero endpoint row was accepted.')

    result = {'status': 'exact finite controls of an all-size polynomial crossing family',
              'residual_column_entries_checked': residual_checks,
              'mode_set_controls': cases,
              'original_tensor_controls': original_tensor_controls(),
              'non_gaussian_control': cubic_control(),
              'parameter_orbit_controls': parameter_orbit_controls(),
              'incorrect_endpoint_rejected': rejected,
              'scope': 'construction and endpoint proof only; no pressure or factoring-exponent improvement'}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
