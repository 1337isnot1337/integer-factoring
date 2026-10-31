#!/usr/bin/env python3
"""Exact terminal binary operator cores for n=2^k+1.

The old 2^k-mode reflected Pascal change leaves only r+1 XV modes.
The existing identity direction j=0 is included with identity weight,
so the zero core is square without a formal Gaussian or division.
All polynomial sources below are complete; no series cutoff is used.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product

import numpy as np

from coupled_identity_mode_filters import add as add2, multiply as multiply2, shift, algebra_multiply
from coupled_identity_operator_filters import _components, states
from coupled_identity_weyl_profiles import determinant, adjugate
from compiled_f4_identity_operators import _compile_matrices, _orders
from f4_identity_operator_filters import normalized_kernel
from f4_repaired_identity_filters import MUL, add, multiply, evaluate_one, laurent_orders
from search_coupled_identity_operator_filters import Pressure


@lru_cache(maxsize=16)
def geometry_data(n, modes):
    parent = n-1
    if parent < 2 or parent & (parent-1) or not 1 <= modes <= parent:
        raise ValueError('Use n=2^k+1, k>=1, and between one and n−1 terminal modes.')
    selected = tuple(range(n-modes, n))
    augmented = selected+(0,)
    retained = tuple(range(modes-1))+(parent-1, parent)
    # For a power-of-two parent, (1+tau)^(parent−1) has every coefficient one.
    old_unit = {degree: 1 for degree in range(parent)}
    fx, fv = [], []
    for i in retained:
        xr, vr = [], []
        for j in augmented:
            if i == parent:
                value = {parent if j == parent else 0: 1}
                xr.append(value)
                vr.append(value)
            else:
                present = i == parent-1-(j if j < parent else 0)
                xr.append({0: 1} if present else {})
                vr.append(dict(old_unit) if present else {})
        fx.append(tuple(xr))
        fv.append(tuple(vr))
    assert determinant(fx, 2).get(0) == determinant(fv, 2).get(0) == 1
    return {'n': n, 'parent': parent, 'modes': modes,
            'selected': selected, 'augmented': augmented, 'retained': retained,
            'augmented_x': tuple(fx), 'augmented_v': tuple(fv),
            'unit': {0: 1}}


@lru_cache(maxsize=16)
def weyl_data(n, modes):
    data = geometry_data(n, modes)
    selected, retained = data['augmented'], data['retained']
    fx, fv = data['augmented_x'], data['augmented_v']
    size = modes+1
    dx, dv = determinant(fx, 2), determinant(fv, 2)
    delta = multiply2({0: 1, 1: 1}, multiply2(dx, dv, 2), 2)
    ax, av = adjugate(fx, 2), adjugate(fv, 2)
    numerator = []
    for i in range(size):
        row = []
        for j in range(size):
            value = dict(delta) if i == j else {}
            for k in range(size):
                term = multiply2(ax[i][k], av[j][k], 2)
                value = add2(value, shift(term, selected[i]+retained[k]+1, 2), 2)
            assert value.get(0, 0) == int(i == j)
            if value:
                assert min(value) >= max(0, selected[i]-selected[j])
            row.append(value)
        numerator.append(tuple(row))
    return {'selected': selected, 'retained': retained,
            'denominator_unit': delta, 'inverse_parameter_numerator': tuple(numerator)}


@lru_cache(maxsize=8)
def reduced_components(n, modes):
    """Complete source Delta^(r+1) Psi(F tensor I_2) on 2^(r+1) states."""
    data = weyl_data(n, modes)
    size = modes+1
    alphabet = states(2, size)
    lookup = {state: i for i, state in enumerate(alphabet)}
    delta, numerator = data['denominator_unit'], data['inverse_parameter_numerator']
    delta_powers = [{0: 1}]
    for exponent in range(size):
        delta_powers.append(multiply2(delta_powers[-1], delta, 2))
    forms = []
    for j in range(size):
        form = {}
        for i in range(size):
            powers = [0]*size
            powers[i] = 1
            if numerator[j][i]:
                form[tuple(powers)] = numerator[j][i]
        forms.append(form)
    cache = {(0,)*size: {(0,)*size: {0: 1}}}

    def monomial(powers):
        if powers not in cache:
            j = next(i for i, digit in enumerate(powers) if digit)
            previous = list(powers)
            previous[j] = 0
            cache[powers] = algebra_multiply(forms[j], monomial(tuple(previous)), 2)
        return cache[powers]

    result = {}
    for b in states(2, modes):
        for a in states(2, modes):
            if sum(a) != sum(b):
                continue
            matrix = {}
            for extra in (0, 1):
                aa, bb = a+(extra,), b+(extra,)
                polynomial = {
                    powers: multiply2(value, delta_powers[sum(bb)], 2)
                    for powers, value in monomial(tuple(1-digit for digit in bb)).items()}
                derivatives = tuple(1-digit for digit in aa)
                for column, source in enumerate(alphabet):
                    for powers, value in polynomial.items():
                        total = tuple(x+y for x, y in zip(source, powers))
                        if any(digit >= 2 for digit in total):
                            continue
                        if any(digit < derivative for digit, derivative in zip(total, derivatives)):
                            continue
                        target = tuple(digit-derivative for digit, derivative in zip(total, derivatives))
                        entry = (lookup[target], column)
                        matrix[entry] = add2(matrix.get(entry, {}), value, 2)
            result[a, b] = {entry: value for entry, value in matrix.items() if value}
    return result


def source_matrix(n, ramification, modes, kernel):
    if ramification <= 0 or ramification % 2 == 0:
        raise ValueError('Use a positive odd ramification.')
    alphabet = states(2, modes+1)
    selected = geometry_data(n, modes)['augmented']
    components = reduced_components(n, modes)
    matrix = [[{} for column in alphabet] for row in alphabet]
    for index, scalar in kernel.items():
        if index not in components:
            raise ValueError('The kernel must preserve total binary occupation.')
        for (row, column), polynomial in components[index].items():
            weight = ramification*sum(j*digit for j, digit in zip(selected, alphabet[row]))
            value = {ramification*degree+weight: coefficient for degree, coefficient in polynomial.items()}
            matrix[row][column] = add(matrix[row][column], multiply(value, scalar))
    return matrix


def zero_profiles(n, ramification, modes, kernel):
    data = geometry_data(n, modes)
    kernel, inverse = normalized_kernel(modes, kernel, data['selected'])
    alphabet = states(2, modes+1)
    matrix = source_matrix(n, ramification, modes, kernel)
    result = defaultdict(list)
    for occupation in range(modes+2):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        block = [[matrix[i][j] for j in indices] for i in indices]
        orders = tuple(2*value for value in laurent_orders(block)[0])
        for auxiliary in alphabet:
            adjustment = -2*ramification*sum(j*digit for j, digit in zip(data['augmented'], auxiliary))
            result[occupation-sum(auxiliary)].extend(order+adjustment for order in orders)
    assert sum(map(len, result.values())) == 2**(2*(modes+1))
    return tuple({'charge': charge, 'orders': tuple(sorted(values))}
                 for charge, values in sorted(result.items()))


def infinity_profiles(n, ramification, modes, kernel):
    selected = geometry_data(n, modes)['selected']
    kernel, inverse = normalized_kernel(modes, kernel, selected)
    alphabet = states(2, modes)
    matrix = [[{ramification*sum((2*j+1-2*n)*digit for j, digit in zip(selected, column))-2*degree: value
                for degree, value in kernel.get((row, column), {}).items()}
               for column in alphabet] for row in alphabet]
    return laurent_orders(matrix)[0]


def compile_family(n, ramification, modes, constant, perturbations):
    """Compile complete F4 sources for a finite affine kernel family."""
    data = geometry_data(n, modes)
    selected, augmented = data['selected'], data['augmented']
    alphabet, augmented_alphabet = states(2, modes), states(2, modes+1)
    zero_sources, infinity_sources, endpoint_sources = [], [], []
    for kernel in (constant,)+tuple(perturbations):
        zero_sources.append(source_matrix(n, ramification, modes, kernel))
        infinity_sources.append([
            [{ramification*sum((2*j+1-2*n)*digit for j, digit in zip(selected, column))-2*degree: value
              for degree, value in kernel.get((row, column), {}).items()}
             for column in alphabet] for row in alphabet])
        endpoint_sources.append(tuple(evaluate_one(kernel.get((row, column), {}))
                                      for row in alphabet for column in alphabet))
    zero_compiled, infinity_compiled = [], []
    for occupation in range(modes+2):
        indices = [i for i, state in enumerate(augmented_alphabet) if sum(state) == occupation]
        zero_compiled.append((occupation, _compile_matrices(
            [[[matrix[i][j] for j in indices] for i in indices] for matrix in zero_sources])))
    for occupation in range(modes+1):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        infinity_compiled.append((occupation, _compile_matrices(
            [[[matrix[i][j] for j in indices] for i in indices] for matrix in infinity_sources])))
    auxiliary = [(sum(state), -2*ramification*sum(j*digit for j, digit in zip(augmented, state)))
                 for state in augmented_alphabet]

    @lru_cache(maxsize=4096)
    def validate_endpoint(values):
        endpoint = {(row, column): {0: values[i*len(alphabet)+j]}
                    for i, row in enumerate(alphabet) for j, column in enumerate(alphabet)
                    if values[i*len(alphabet)+j]}
        normalized_kernel(modes, endpoint, selected)

    def evaluate(colors):
        if len(colors) != len(perturbations) or any(color not in (0, 1, 2, 3) for color in colors):
            raise ValueError('Use one F4 color per affine perturbation.')
        endpoint = list(endpoint_sources[0])
        for source, color in zip(endpoint_sources[1:], colors):
            if color:
                for i, value in enumerate(source):
                    endpoint[i] ^= MUL[color][value]
        validate_endpoint(tuple(endpoint))
        zero = defaultdict(list)
        for occupation, compiled in zero_compiled:
            orders = tuple(2*value for value in _orders(compiled, colors))
            for auxiliary_occupation, adjustment in auxiliary:
                zero[occupation-auxiliary_occupation].extend(order+adjustment for order in orders)
        infinity = tuple(sorted(value for occupation, compiled in infinity_compiled
                                for value in _orders(compiled, colors)))
        blocks = tuple({'charge': charge, 'orders': tuple(sorted(values))}
                       for charge, values in sorted(zero.items()))
        return blocks, infinity

    return evaluate


def profile_data(n, ramification, charge, modes, kernel, blocks=None):
    if not isinstance(charge, (int, Fraction)):
        raise ValueError('Use an exact rational aggregate charge.')
    data = geometry_data(n, modes)
    if blocks is None:
        blocks = zero_profiles(n, ramification, modes, kernel)
    return {'prime': 2, 'coefficient_field': 'F4', 'n': n, 'ramification': ramification,
            'modes': modes, 'charge': charge, 'identity_indices': data['selected'],
            'active_zero_identity_indices': data['augmented'], 'retained_xv_indices': data['retained'],
            'zero_active': tuple(sorted(value+ramification*(1-charge)*block['charge']
                                        for block in blocks for value in block['orders'])),
            'zero_charge_blocks_at_J1': blocks,
            'zero_spectators': tuple(ramification*(2*j+1-charge) for j in range(n) if j not in data['augmented'])
                +tuple(ramification*(2*j+1+charge) for j in range(n) if j not in data['retained']),
            'infinity_active': infinity_profiles(n, ramification, modes, kernel),
            'infinity_spectators': tuple(ramification*(2*j+1-2*n) for j in range(2*n) if j not in data['selected']),
            'paired_physical_first': tuple(ramification*(4*j-2*(n-1)+charge) for j in range(n)),
            'paired_physical_second': tuple(ramification*(4*j-2*(n-1)-charge) for j in range(n))}


class AbovePowerPressure(Pressure):
    def __init__(self, n, ramification, modes, blocks, infinity):
        data = geometry_data(n, modes)
        self.n, self.ramification, self.modes = n, ramification, modes
        self.selected, self.retained = data['augmented'], data['retained']
        self.infinity_selected = data['selected']
        self.zero = np.array([value/ramification for block in blocks for value in block['orders']])
        self.charges = np.array([block['charge'] for block in blocks for value in block['orders']])
        self.infinity = np.array(infinity)/ramification

    def grid(self, scales, charges, relative=False):
        # The parent removes j=0 at both places. Restore its infinity factor:
        # this physical direction is included only in the zero active core.
        scales, charges = np.broadcast_arrays(np.asarray(scales), np.asarray(charges))
        return super().grid(scales, charges, relative)+np.logaddexp(0, (2*self.n-1)*scales)


@lru_cache(maxsize=4)
def original_components(n, modes):
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel
    selected = geometry_data(n, modes)['selected']
    matrix = pascal_kernel(range(n), 2*n, 2)
    gaussian = {powers: {degree: value for degree, value in enumerate(polynomial) if value}
                for powers, polynomial in direct_kernel(matrix, 2).items()}
    left, right = [], []
    for j in selected:
        first, second = {}, {}
        for index, polynomial in [(n+j, (1,))]+[(i, matrix[i][j]) for i in range(n)]:
            if polynomial:
                powers = [0]*(4*n)
                powers[index] = 1
                first[tuple(powers)] = {degree: value for degree, value in enumerate(polynomial) if value}
        for index, polynomial in [(2*n+j, (1,))]+[(3*n+i, matrix[n+j][n+i]) for i in range(n)]:
            if polynomial:
                powers = [0]*(4*n)
                powers[index] = 1
                second[tuple(powers)] = {degree: value for degree, value in enumerate(polynomial) if value}
        left.append(first)
        right.append(second)
    return _components(gaussian, left, right, 2), gaussian


def polynomial_seed(n, ramification, modes, kernel):
    selected = geometry_data(n, modes)['selected']
    kernel, inverse = normalized_kernel(modes, kernel, selected)
    components, gaussian = original_components(n, modes)
    clearing = max(0, -min(degree for value in kernel.values() for degree in value))
    tensor = {}
    for index, scalar in kernel.items():
        scalar = {degree+clearing: value for degree, value in scalar.items()}
        for powers, polynomial in components[index].items():
            value = {ramification*degree: coefficient for degree, coefficient in polynomial.items()}
            tensor[powers] = add(tensor.get(powers, {}), multiply(value, scalar))
    alphabet = states(2, modes)
    lookup = {state: i for i, state in enumerate(alphabet)}
    repaired = {}
    for powers, value in tensor.items():
        column = lookup[tuple(powers[n+j] for j in selected)]
        for row, state in enumerate(alphabet):
            if inverse[row][column]:
                target = list(powers)
                for j, digit in zip(selected, state):
                    target[n+j] = digit
                target = tuple(target)
                repaired[target] = add(repaired.get(target, {}), value, inverse[row][column])
    repaired = {powers: value for powers, value in repaired.items() if value}
    assert all(evaluate_one(repaired.get(powers, {})) == evaluate_one(gaussian.get(powers, {}))
               for powers in set(repaired)|set(gaussian))
    bound = ramification*(max(degree for value in gaussian.values() for degree in value)+2*n*modes)
    bound += max(degree+clearing for value in kernel.values() for degree in value)
    assert max(degree for value in repaired.values() for degree in value) <= bound
    return repaired, clearing, gaussian
