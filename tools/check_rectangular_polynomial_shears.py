#!/usr/bin/env python3
"""Endpoint-identity Laurent polynomial core shears and independent controls.

For P(1)=0, I+P(w)E_rs preserves the ramified endpoint for every r!=s.
Actual height may increase. The finite controls audit an ordinary all-length
proof; they do not certify any floating pressure or factoring exponent.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from check_rectangular_f4_repair import (
    ZERO, add, multiply, times_binary, binary_power, build,
)
from check_rectangular_identity_transports import IDENTITY, product_letters
from check_rectangular_refined_repairs import height, _verified_repair
from rectangular_coordinate_transforms import (
    transform_coordinate_matrix, linear_combination_coordinate_matrices,
)


def validated_terms(row, column, terms):
    if row == column or not 0 <= row < 12 or not 0 <= column < 12:
        raise ValueError('A shear needs distinct core indices between zero and eleven.')
    coefficients = {}
    for coefficient, exponent in terms:
        if coefficient not in (0, 1, 2, 3) or not isinstance(exponent, int):
            raise ValueError('A Laurent term needs an F4 coefficient and an integer exponent.')
        coefficients[exponent] = coefficients.get(exponent, 0) ^ coefficient
    answer = [(coefficient, exponent) for exponent, coefficient in sorted(coefficients.items()) if coefficient]
    value_at_one = 0
    for coefficient, _ in answer:
        value_at_one ^= coefficient
    if value_at_one:
        raise ValueError('The polynomial shear must have P(1)=0.')
    return answer


def polynomial_shear_seed(seed, row, column, terms, side='output', check_types=True):
    """Apply H^-1(I+P E_rs)H directly to a polynomial numerator tensor.

    Negative Laurent powers are cleared by one global w-power, recorded in
    global_w_power. The input convention is source row, target column.
    """
    terms = validated_terms(row, column, terms)
    if side not in ('output', 'input'):
        raise ValueError('Choose output or input for the physical core shear.')
    if not terms:
        return dict(seed)
    dh = height(column)-height(row)
    padding = max(0, -dh)
    extra_w = max(0, -min(exponent for _, exponent in terms))
    polynomial = tuple(sum(((coefficient >> plane)&1) << (exponent+extra_w)
                           for coefficient, exponent in terms) for plane in (0, 1))
    factor = times_binary(polynomial, binary_power(3, padding+dh))
    identity = binary_power(3, padding) << extra_w
    values = [times_binary(value, identity) for value in seed['values']]
    A, B, bulk = seed['q_x'], seed['q_u'], seed['bulk_size']
    def position(x, y, u, v):
        return x+A*(y+B*(u+B*v))
    for x, y, u, v in product_letters(A, B):
        if side == 'output':
            if v//bulk+4*(y//bulk) != row:
                continue
            source = position(x, (column//4)*bulk+y % bulk, u, (column % 4)*bulk+v % bulk)
        else:
            if x//bulk+4*(u//bulk) != column:
                continue
            source = position((row % 4)*bulk+x % bulk, y, (row//4)*bulk+u % bulk, v)
        target = position(x, y, u, v)
        values[target] = add(values[target], multiply(seed['values'][source], factor))
    uncancelled = seed['denominator_degree']+padding
    provisional = dict(seed, denominator_degree=uncancelled)
    endpoint = seed['ramified_endpoint_matrix']
    result = _verified_repair(provisional, values, endpoint,
                              [(height(i), 0) for i in range(12)], check_types)
    specification = {'side': side, 'row': row, 'column': column, 'terms': terms}
    result.update({'repair_kind': 'mixed_polynomial_identity_shear',
                   'polynomial_shear': specification,
                   'polynomial_shear_chain': seed.get('polynomial_shear_chain', [])+[specification],
                   'constant_base_repair_matrix': seed.get('constant_base_repair_matrix', seed.get('repair_matrix')),
                   'global_w_power': seed.get('global_w_power', 0)+extra_w,
                   'uncancelled_denominator_degree': uncancelled,
                   'cancelled_c_order': uncancelled-result['denominator_degree']})
    result.pop('repair_matrix', None)
    return result


def polynomial_coordinate_shear(compiled, data, row, column, terms,
                                side='output', unit_step=None):
    """Apply an endpoint-identity shear through two exact cached matrices.

    terms use the global variable w. A supplied unit_step must be the cached
    result for I+E_rs with the same base and side. This function preserves
    the original endpoint; the caller must certify that original endpoint.
    """
    terms = validated_terms(row, column, terms)
    if side not in ('output', 'input'):
        raise ValueError('Choose output or input for the physical core shear.')
    if not terms:
        return compiled
    if unit_step is None:
        gate = [line[:] for line in IDENTITY]
        gate[row][column] = 1
        unit_step = transform_coordinate_matrix(compiled, data,
                                                V=gate if side == 'output' else None,
                                                W=gate if side == 'input' else None)
    sign = -1 if data['place'] == 'infinity' else 1
    local = [(coefficient, sign*exponent) for coefficient, exponent in terms]
    result = linear_combination_coordinate_matrices([
        (compiled, [(1, 0)]+local), (unit_step, local)])
    result[2]['endpoint_preserving_polynomial_shear'] = {
        'side': side, 'row': row, 'column': column, 'terms': terms,
        'value_at_one': 0, 'base_endpoint_still_requires_certificate': True}
    return result


def compare_inverse(base, restored):
    delta = restored.get('global_w_power', 0)-base.get('global_w_power', 0)
    assert delta >= 0
    left = binary_power(3, base['denominator_degree'])
    right = binary_power(3, restored['denominator_degree'])
    for original, actual in zip(base['values'], restored['values']):
        assert times_binary(actual, left) == times_binary(tuple(plane << delta for plane in original), right)
    return len(base['values'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--dense', action='store_true')
    args = parser.parse_args()
    base = build((0, 1), 2, 1, 2)
    controls = []
    inverse_entries = 0
    for row in range(12):
        for column in range(12):
            if row == column:
                continue
            for coefficient, ell in ((1, -2), (2, 0), (3, 3)):
                terms = [(coefficient, ell), (coefficient, ell+1)]
                changed = polynomial_shear_seed(base, row, column, terms)
                restored = polynomial_shear_seed(changed, row, column, terms, check_types=False)
                inverse_entries += compare_inverse(base, restored)
                assert changed['denominator_degree'] <= 2
                controls.append({'row': row, 'column': column, 'coefficient': coefficient,
                                 'ell': ell, 'actual_height_change': height(row)-height(column),
                                 'denominator_degree': changed['denominator_degree'],
                                 'global_w_power': changed['global_w_power'],
                                 'ramified_endpoint_checks': changed['ramified_endpoint_checks'],
                                 'two_copy_type_checks': changed['two_copy_type_checks'],
                                 'coefficient_sha256': changed['coefficient_sha256']})
    input_controls = []
    for row, column, coefficient, ell in ((10, 1, 2, 0), (3, 8, 3, -2), (0, 11, 1, 3)):
        terms = [(coefficient, ell), (coefficient, ell+1)]
        changed = polynomial_shear_seed(base, row, column, terms, side='input')
        restored = polynomial_shear_seed(changed, row, column, terms, side='input', check_types=False)
        inverse_entries += compare_inverse(base, restored)
        input_controls.append({'row': row, 'column': column, 'coefficient': coefficient, 'ell': ell,
                               'denominator_degree': changed['denominator_degree'],
                               'ramified_endpoint_checks': changed['ramified_endpoint_checks'],
                               'two_copy_type_checks': changed['two_copy_type_checks']})
    dense = []
    if args.dense:
        from check_rectangular_defect_core_basis import data_for
        from check_rectangular_defect_coordinates import point_control, weighted_controls
        from rectangular_defect_coordinates import coordinate_matrix
        large_base = build(tuple(range(4)), 4, 1, 2, check_types=False)
        seeds = [(row, column, terms, polynomial_shear_seed(large_base, row, column, terms, check_types=False))
                 for row, column, terms in ((10, 1, [(2, 0), (2, 1)]),
                                            (3, 8, [(3, -2), (3, -1)]),
                                            (0, 11, [(1, 2), (1, 3)]))]
        for place in ('infinity', 'zero'):
            data = data_for(place, 4)
            compiled = coordinate_matrix(data)
            for row, column, terms, seed in seeds:
                transformed = polynomial_coordinate_shear(compiled, data, row, column, terms)
                power = seed['global_w_power']*(-1 if place == 'infinity' else 1)
                if power:
                    transformed = linear_combination_coordinate_matrices([(transformed, [(1, power)])])
                points = [point_control(data, transformed, z, seed) for z in (2, 3)]
                profiles = weighted_controls(data, transformed, 2, seed)
                dense.append({'place': place, 'M': 4, 'row': row, 'column': column,
                              'terms': terms, 'point_checks': points, 'weighted_profiles': profiles})
                print(json.dumps({'place': place, 'row': row, 'column': column,
                                  'entry_checks': sum(point['entries_checked'] for point in points)}), flush=True)
    result = {'status': 'exact endpoint-identity polynomial core shears; no pressure claim',
              'output_shear_cases': controls, 'inverse_entries_checked': inverse_entries,
              'input_shear_cases': input_controls,
              'dense_original_tensor_controls': dense}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
