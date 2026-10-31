#!/usr/bin/env python3
"""Exact two-cut row-lattice synthesis for occupation-preserving F4 kernels.

At zero, impose row orders in a regular Gaussian Smith basis of the
Weyl matrix. At infinity, impose row orders after an independent constant
binary row change. All Laurent coefficients, including off-diagonal
occupation transfers, are unknown. The endpoint must preserve the physical
Y energy filtration and have nonzero diagonal, so its constant inverse
repairs the endpoint without changing either profile.
"""
from __future__ import annotations

from functools import lru_cache

from coupled_identity_operator_filters import _inverse_matrix, states
from coupled_identity_weyl_profiles import reduced_components, weyl_data
from f4_identity_operator_filters import normalized_kernel
from f4_repaired_identity_filters import add
from weyl_lattice_filter_synthesis import (
    binary_image_basis, covering_pair, endpoint_image, matrix_multiply,
    nonzero_endpoint_lift, polynomial_degrees, reduce_row, regular_smith_basis,
)


def xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


@lru_cache(maxsize=16)
def operator_basis(n, ramification=1, modes=3, geometry='split:1:2'):
    """Every operator coefficient in a complete polynomial zero Smith basis."""
    if modes not in (2, 3):
        raise ValueError('This prototype supports two or three binary identity directions.')
    alphabet = states(2, modes)
    size = len(alphabet)
    data = weyl_data(2, n, modes, geometry)
    components = reduced_components(2, n, modes, geometry)
    sources = {}
    for index, entries in components.items():
        matrix = [[0 for _ in alphabet] for _ in alphabet]
        for (row, column), polynomial in entries.items():
            weight = ramification * sum(j * digit for j, digit in zip(data['selected'], alphabet[row]))
            matrix[row][column] = sum(1 << (ramification * degree + weight)
                                      for degree, coefficient in polynomial.items() if coefficient)
        sources[index] = matrix
    transformed = {index: {} for index in sources}
    gaussian_orders, occupations, blocks = [], [], []
    offset = 0
    for occupation in range(modes + 1):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        gaussian = [[xor(sources[state, state][i][j] for state in alphabet)
                     for j in indices] for i in indices]
        orders, left, right, diagonal = regular_smith_basis(gaussian)
        for index, matrix in sources.items():
            value = matrix_multiply(matrix_multiply(left, [[matrix[i][j] for j in indices]
                                                          for i in indices]), right)
            transformed[index].update({(offset+i, j): coefficient
                                       for i, row in enumerate(value)
                                       for j, coefficient in enumerate(row) if coefficient})
        gaussian_orders.extend(orders)
        occupations.extend([occupation] * len(indices))
        blocks.append({'occupation': occupation, 'indices': tuple(indices),
                       'row_offset': offset, 'orders': orders})
        offset += len(indices)
    assert offset == size
    infinity = tuple(ramification * sum((2*j+1-2*n)*digit
                                        for j, digit in zip(data['selected'], state))
                     for state in alphabet)
    return {'states': alphabet, 'selected': data['selected'], 'retained': data['retained'],
            'coefficients': transformed, 'gaussian_orders': tuple(gaussian_orders),
            'row_occupations': tuple(occupations), 'infinity_orders': infinity,
            'blocks': tuple(blocks)}


def _row_basis(data, basis):
    size = len(data['states'])
    if basis is None:
        basis = tuple(tuple(int(i == j) for j in range(size)) for i in range(size))
    basis = tuple(tuple(row) for row in basis)
    if len(basis) != size or any(len(row) != size for row in basis):
        raise ValueError('The infinity row basis has the wrong matrix size.')
    if any(value not in (0, 1) or (value and sum(data['states'][i]) != sum(data['states'][j]))
           for i, row in enumerate(basis) for j, value in enumerate(row)):
        raise ValueError('Use a binary row basis preserving total occupation.')
    return basis, tuple(tuple(row) for row in _inverse_matrix(basis, 2))


def coefficient_system(n, ramification, modes, geometry, lower_power,
                       zero_targets, infinity_targets, infinity_basis=None, upper_power=None):
    data = operator_basis(n, ramification, modes, geometry)
    alphabet = data['states']
    size = len(alphabet)
    if not isinstance(lower_power, int):
        raise ValueError('Use one integer lower Laurent power.')
    if upper_power is not None and (not isinstance(upper_power, int) or upper_power < lower_power):
        raise ValueError('An optional upper Laurent power must be an integer at least the lower bound.')
    if len(zero_targets) != size or len(infinity_targets) != size:
        raise ValueError('Use one target per zero Smith row and per infinity row.')
    if any(not isinstance(value, int) for value in tuple(zero_targets)+tuple(infinity_targets)):
        raise ValueError('The row targets must be integers.')
    basis, inverse = _row_basis(data, infinity_basis)
    sources, caps = {}, {}
    for r, row in enumerate(alphabet):
        for c, column in enumerate(alphabet):
            if sum(row) != sum(column):
                continue
            caps[r, c] = (data['infinity_orders'][c]-infinity_targets[r])//2
            if upper_power is not None:
                caps[r, c] = min(caps[r, c], upper_power)
            coefficients = {}
            for a, original_row in enumerate(alphabet):
                if inverse[a][r]:
                    for entry, coefficient in data['coefficients'][original_row, column].items():
                        coefficients[entry] = coefficients.get(entry, 0) ^ coefficient
            sources[r, c] = {entry: value for entry, value in coefficients.items() if value}
    variables = tuple((r, c, degree) for (r, c), cap in caps.items()
                      for degree in range(lower_power, cap+1))
    endpoint_masks = [0] * size
    forbidden = {}
    equations = {}
    energies = tuple(sum(j*digit for j, digit in zip(data['selected'], state)) for state in alphabet)
    for variable, (r, c, degree) in enumerate(variables):
        bit = 1 << variable
        for a in range(size):
            if not inverse[a][r]:
                continue
            if a == c:
                endpoint_masks[a] ^= bit
            elif energies[a] <= energies[c]:
                forbidden[a, c] = forbidden.get((a, c), 0) ^ bit
        for (row, column), coefficient in sources[r, c].items():
            cutoff = zero_targets[row]-degree
            if cutoff <= 0:
                continue
            for power in polynomial_degrees(coefficient & ((1 << cutoff)-1)):
                key = row, column, power+degree
                equations[key] = equations.get(key, 0) ^ bit
    pivots = {}
    constraints = tuple(value for value in equations.values() if value) + tuple(
        value for value in forbidden.values() if value)
    for equation in constraints:
        while equation:
            pivot = equation.bit_length()-1
            if pivot in pivots:
                equation ^= pivots[pivot]
            else:
                pivots[pivot] = equation
                break
    return {'data': data, 'variables': variables, 'pivots': pivots,
            'endpoint_masks': tuple(endpoint_masks), 'constraints': len(constraints),
            'sources': sources, 'caps': caps, 'infinity_basis': basis,
            'infinity_inverse': inverse}


def synthesize(n, ramification, modes, geometry, lower_power,
               zero_targets, infinity_targets, infinity_basis=None, construct=True, upper_power=None):
    """Return an exact witness, or None for this prescribed pair of lattices."""
    zero_targets, infinity_targets = tuple(zero_targets), tuple(infinity_targets)
    system = coefficient_system(n, ramification, modes, geometry, lower_power,
                                zero_targets, infinity_targets, infinity_basis, upper_power)
    variables, pivots = system['variables'], system['pivots']
    masks, data = system['endpoint_masks'], system['data']
    remainders = tuple(reduce_row(mask, pivots) for mask in masks)
    if not all(remainders) or covering_pair(binary_image_basis(remainders), len(masks)) is None:
        return None
    if not construct:
        return True
    image = endpoint_image(variables, pivots, masks)
    lifted = nonzero_endpoint_lift(image, len(masks))
    assert lifted is not None
    endpoint, lift = lifted
    kernel, transformed = {}, {}
    inverse, alphabet = system['infinity_inverse'], data['states']
    for variable, coefficient in lift.items():
        r, c, degree = variables[variable]
        assert lower_power <= degree <= system['caps'][r, c]
        for a in range(len(alphabet)):
            if inverse[a][r]:
                index = alphabet[a], alphabet[c]
                kernel[index] = add(kernel.get(index, {}), {degree: coefficient})
        for entry, polynomial in system['sources'][r, c].items():
            value = {power+degree: coefficient for power in polynomial_degrees(polynomial)}
            transformed[entry] = add(transformed.get(entry, {}), value)
    kernel = {index: value for index, value in kernel.items() if value}
    kernel, endpoint_inverse = normalized_kernel(modes, kernel, data['selected'])
    assert all(not value or min(value) >= zero_targets[row]
               for (row, column), value in transformed.items())
    return {'n': n, 'ramification': ramification, 'modes': modes, 'geometry': geometry,
            'lower_power': lower_power, 'zero_targets': zero_targets,
            'upper_power': upper_power,
            'infinity_targets': infinity_targets, 'infinity_basis': system['infinity_basis'],
            'gaussian_zero_orders': data['gaussian_orders'],
            'gaussian_infinity_orders': data['infinity_orders'],
            'coefficient_variables': len(variables), 'constraints': system['constraints'],
            'homogeneous_rank': len(pivots), 'endpoint_image_dimension': len(image),
            'endpoint_values': endpoint, 'kernel': kernel,
            'verified_zero_entries': len(transformed)}


def bound_blocks(n, ramification, modes, geometry, zero_targets):
    from collections import defaultdict
    data = operator_basis(n, ramification, modes, geometry)
    result = defaultdict(list)
    for auxiliary in data['states']:
        shift = -2*ramification*sum(j*digit for j, digit in zip(data['selected'], auxiliary))
        for occupation, value in zip(data['row_occupations'], zero_targets):
            result[occupation-sum(auxiliary)].append(2*value+shift)
    return tuple({'charge': charge, 'orders': tuple(sorted(orders))}
                 for charge, orders in sorted(result.items()))
