#!/usr/bin/env python3
"""Operator lattice synthesis with F4 infinity flags and fixed F16 lifts."""
from __future__ import annotations

import f4_lattice_linear as linear
from f4_identity_operator_filters import inverse_matrix
from weyl_lattice_filter_synthesis import polynomial_degrees
from weyl_operator_lattice_synthesis import operator_basis as general_operator_basis
from scalar_pascal_operator_basis import operator_basis as scalar_operator_basis


def operator_basis(n, ramification=1, modes=3, geometry='split:1:2'):
    if n > 0 and not n & (n-1):
        return scalar_operator_basis(n, ramification, modes, geometry)
    return general_operator_basis(n, ramification, modes, geometry)


def coefficient_system(n, ramification, modes, geometry, lower_power,
                       zero_targets, infinity_targets, infinity_basis,
                       upper_power=None):
    data = operator_basis(n, ramification, modes, geometry)
    alphabet, size = data['states'], len(data['states'])
    basis = tuple(tuple(row) for row in infinity_basis)
    if len(basis) != size or any(len(row) != size for row in basis):
        raise ValueError('Give one square infinity row basis on the occupation states.')
    if any(value not in (0, 1, 2, 3) or (value and sum(alphabet[i]) != sum(alphabet[j]))
           for i, row in enumerate(basis) for j, value in enumerate(row)):
        raise ValueError('The F4 infinity row basis must preserve total occupation.')
    if len(zero_targets) != size or len(infinity_targets) != size:
        raise ValueError('Give one integer target per zero/infinity row.')
    if any(not isinstance(value, int) for value in tuple(zero_targets)+tuple(infinity_targets)+(lower_power,)):
        raise ValueError('Laurent powers and row targets must be integers.')
    if upper_power is not None and (not isinstance(upper_power, int) or upper_power < lower_power):
        raise ValueError('The optional upper power must be at least the lower power.')
    inverse = inverse_matrix(basis)
    caps, sources = {}, {}
    for r, row in enumerate(alphabet):
        for c, column in enumerate(alphabet):
            if sum(row) != sum(column):
                continue
            cap = (data['infinity_orders'][c]-infinity_targets[r])//2
            caps[r, c] = cap if upper_power is None else min(cap, upper_power)
            coefficients = {}
            for a, original in enumerate(alphabet):
                scalar = inverse[a][r]
                if scalar:
                    for entry, polynomial in data['coefficients'][original, column].items():
                        coefficients[entry] = linear.add(coefficients.get(entry, (0, 0)),
                                                          linear.scale((polynomial, 0), scalar))
            sources[r, c] = {entry: value for entry, value in coefficients.items() if value[0] | value[1]}
    variables = tuple((r, c, degree) for (r, c), cap in caps.items()
                      for degree in range(lower_power, cap+1))
    endpoints, forbidden, equations = [(0, 0)] * size, {}, {}
    energies = tuple(sum(j*digit for j, digit in zip(data['selected'], state)) for state in alphabet)
    for variable, (r, c, degree) in enumerate(variables):
        bit = 1 << variable
        for a in range(size):
            scalar = inverse[a][r]
            if not scalar:
                continue
            value = linear.scale((bit, 0), scalar)
            if a == c:
                endpoints[a] = linear.add(endpoints[a], value)
            elif energies[a] <= energies[c]:
                forbidden[a, c] = linear.add(forbidden.get((a, c), (0, 0)), value)
        for (row, column), polynomial in sources[r, c].items():
            cutoff = zero_targets[row]-degree
            if cutoff <= 0:
                continue
            support = (polynomial[0] | polynomial[1]) & ((1 << cutoff)-1)
            for power in polynomial_degrees(support):
                value = linear.scale((bit, 0), linear.coefficient(polynomial, power))
                key = row, column, power+degree
                equations[key] = linear.add(equations.get(key, (0, 0)), value)
    rows = tuple(value for value in tuple(equations.values())+tuple(forbidden.values()) if value[0] | value[1])
    pivots = linear.eliminate(rows)
    remainders = tuple(linear.reduce_row(value, pivots) for value in endpoints)
    return {'data': data, 'variables': variables, 'pivots': pivots, 'equations': rows,
            'endpoint_rows': tuple(endpoints), 'endpoint_remainders': remainders,
            'sources': sources, 'caps': caps, 'infinity_basis': basis,
            'infinity_inverse': inverse}


def synthesize(n, ramification, modes, geometry, lower_power,
               zero_targets, infinity_targets, infinity_basis,
               construct=True, field_size=16, upper_power=None):
    """Return a complete F4/F16 kernel or None for the selected two lattices."""
    if field_size not in (4, 16):
        raise ValueError('The supported coefficient fields are F4 and F16.')
    system = coefficient_system(n, ramification, modes, geometry, lower_power,
                                zero_targets, infinity_targets, infinity_basis, upper_power)
    remainders = system['endpoint_remainders']
    if any(not (a | b) for a, b in remainders):
        return None
    if not construct and field_size >= len(remainders):
        return True
    lifted = linear.full_support_lift(remainders, system['pivots'], field_size)
    if lifted is None:
        return None
    if not construct:
        return True
    endpoint, assignment, dimension = lifted
    assert linear.verify_assignment(system['equations'], assignment)
    kernel, transformed = {}, {}
    alphabet = system['data']['states']
    for variable, value in assignment.items():
        r, c, degree = system['variables'][variable]
        for a, original in enumerate(alphabet):
            coefficient = linear.f16_scale_f4(value, system['infinity_inverse'][a][r])
            if coefficient:
                polynomial = kernel.setdefault((original, alphabet[c]), {})
                polynomial[degree] = polynomial.get(degree, 0) ^ coefficient
        for entry, polynomial in system['sources'][r, c].items():
            result = transformed.setdefault(entry, {})
            for power in polynomial_degrees(polynomial[0] | polynomial[1]):
                coefficient = linear.f16_scale_f4(value, linear.coefficient(polynomial, power))
                result[power+degree] = result.get(power+degree, 0) ^ coefficient
    kernel = {index: {degree: value for degree, value in polynomial.items() if value}
              for index, polynomial in kernel.items()}
    kernel = {index: polynomial for index, polynomial in kernel.items() if polynomial}
    assert all(not value or min(value) >= zero_targets[row]
               for (row, column), polynomial in transformed.items()
               if (value := {degree: coefficient for degree, coefficient in polynomial.items() if coefficient}))
    actual_endpoint = []
    energies = {state: sum(j*digit for j, digit in zip(system['data']['selected'], state)) for state in alphabet}
    for row in alphabet:
        diagonal = 0
        for coefficient in kernel.get((row, row), {}).values():
            diagonal ^= coefficient
        actual_endpoint.append(diagonal)
    assert tuple(actual_endpoint) == endpoint and all(endpoint)
    for (row, column), polynomial in kernel.items():
        if row != column and energies[row] <= energies[column]:
            value = 0
            for coefficient in polynomial.values():
                value ^= coefficient
            assert value == 0
    return {'n': n, 'ramification': ramification, 'modes': modes, 'geometry': geometry,
            'coefficient_field': f'F{field_size}', 'lower_power': lower_power,
            'upper_power': upper_power, 'zero_targets': tuple(zero_targets),
            'infinity_targets': tuple(infinity_targets), 'infinity_basis': system['infinity_basis'],
            'coefficient_variables': len(system['variables']), 'homogeneous_rank': len(system['pivots']),
            'coefficient_equations': len(system['equations']),
            'endpoint_image_dimension': dimension, 'endpoint_values': endpoint,
            'kernel': kernel}
