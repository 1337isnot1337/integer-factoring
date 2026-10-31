#!/usr/bin/env python3
"""Exact polynomial joint-filter synthesis with a sufficient zero-cost bound.

Split each Gaussian diagonal energy E as theta*E+(1-theta)*E. Requiring
the corresponding rescaled coefficient matrix to be regular is linear
in all scalar state-filter coefficients. Every minor then has at least
the Gaussian prefix valuation, plus the requested common shift.

This sufficient certificate can miss valid low-cost filters. An
infeasible result only excludes the supplied degree bounds and split.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from math import factorial
from pathlib import Path

import coupled_identity_extreme_filters as extreme
import coupled_identity_mode_filters as initial


@lru_cache(maxsize=16)
def occupation_basis(prime, n, ramification, modes=2, geometry='initial'):
    """Literal tensor coefficient of every independent occupation weight.

    Powers are converted from z to x=z^2. Each entry also records its
    row and column Gaussian energy at aggregate charge J=1.
    """
    if geometry == 'initial':
        selected = tuple(range(modes))
        retained = tuple(range(n - modes, n))
        kernels = initial.active_power_kernels(prime, n, ramification, modes)
    else:
        data = extreme.geometry_data(prime, n, modes, geometry)
        selected, retained = data['selected'], data['retained']
        kernels = extreme.active_power_kernels(prime, n, ramification, modes, geometry)
    states = tuple(product(range(prime), repeat=modes))
    basis = []
    for state in states:
        tensor = {}
        for powers in product(*(range(digit, prime) for digit in state)):
            coefficient = (-1) ** sum(a - b for a, b in zip(powers, state))
            for a, b in zip(powers, state):
                coefficient *= pow(factorial(b) * factorial(a - b), -1, prime)
            for variables, value in kernels[powers].items():
                tensor[variables] = initial.add(tensor.get(variables, {}), value, prime, coefficient)
        result = {}
        for variables, value in tensor.items():
            if not value:
                continue
            if any(power % 2 for power in value):
                raise ArithmeticError('The base-charge active tensor must have integral x=z^2 powers.')
            v, y = variables[:modes], variables[modes:2 * modes]
            u, x = variables[2 * modes:3 * modes], variables[3 * modes:]
            column = ramification * (sum(j * digit for j, digit in zip(selected, u))
                                       + sum((j + 1) * digit for j, digit in zip(retained, v)))
            row = ramification * (sum(j * digit for j, digit in zip(selected, y))
                                    + sum((j + 1) * digit for j, digit in zip(retained, x)))
            result[variables] = ({power // 2: coefficient for power, coefficient in value.items()},
                                 column, row)
        basis.append(result)
    return states, tuple(basis)


def coefficient_system(prime, n, ramification, modes, geometry,
                       common_zero_shift, degrees, balance):
    states, basis = occupation_basis(prime, n, ramification, modes, geometry)
    variables = tuple((state, power) for state, bound in enumerate(degrees)
                      for power in range(bound + 1))
    equations = {}
    for column, (state, power) in enumerate(variables):
        for physical, (polynomial, column_energy, row_energy) in basis[state].items():
            target = (common_zero_shift + balance * row_energy
                      + (1 - balance) * column_energy)
            for degree, value in polynomial.items():
                if degree + power >= target:
                    continue
                key = physical, degree + power
                equation = equations.setdefault(key, {})
                equation[column] = (equation.get(column, 0) + value) % prime
                if not equation[column]:
                    del equation[column]
    return states, variables, tuple(value for value in equations.values() if value)


def sparse_nullspace(equations, columns, prime):
    pivots = {}
    for equation in equations:
        row = dict(equation)
        while row:
            column = min(row)
            coefficient = row[column]
            if column in pivots:
                row = initial.add(row, pivots[column], prime, -coefficient)
            else:
                inverse = pow(coefficient, -1, prime)
                pivots[column] = {j: value * inverse % prime for j, value in row.items()}
                break
    nullspace = []
    for free in range(columns):
        if free in pivots:
            continue
        vector = {free: 1}
        for column, row in sorted(pivots.items(), reverse=True):
            value = -sum(coefficient * vector.get(j, 0)
                         for j, coefficient in row.items() if j != column) % prime
            if value:
                vector[column] = value
        nullspace.append(vector)
    return tuple(nullspace), len(pivots)


def endpoint_witness(nullspace, variables, states, prime):
    image = []
    for lift in nullspace:
        lift = dict(lift)
        endpoint = [0] * len(states)
        for column, coefficient in lift.items():
            endpoint[variables[column][0]] = (endpoint[variables[column][0]] + coefficient) % prime
        for pivot, vector, old_lift in image:
            coefficient = endpoint[pivot]
            if coefficient:
                endpoint = [(a - coefficient * b) % prime for a, b in zip(endpoint, vector)]
                lift = initial.add(lift, old_lift, prime, -coefficient)
        pivot = next((j for j, value in enumerate(endpoint) if value), None)
        if pivot is None:
            continue
        inverse = pow(endpoint[pivot], -1, prime)
        image.append((pivot, tuple(value * inverse % prime for value in endpoint),
                      {j: value * inverse % prime for j, value in lift.items()}))
    target = [1] * len(states)
    witness = {}
    for pivot, vector, lift in image:
        coefficient = target[pivot]
        if coefficient:
            target = [(a - coefficient * b) % prime for a, b in zip(target, vector)]
            witness = initial.add(witness, lift, prime, coefficient)
    return (None if any(target) else witness), len(image)


def verify_witness(prime, n, ramification, modes, geometry,
                   common_zero_shift, degrees, balance, alpha):
    states, basis = occupation_basis(prime, n, ramification, modes, geometry)
    tensor = {}
    targets = {}
    for state, values in zip(states, basis):
        for variables, (polynomial, column, row) in values.items():
            targets[variables] = common_zero_shift + balance * row + (1 - balance) * column
            tensor[variables] = initial.add(tensor.get(variables, {}),
                                            initial.multiply(polynomial, alpha[state], prime), prime)
    assert all(value and max(value) <= bound and sum(value.values()) % prime == 1
               for state, bound in zip(states, degrees) for value in (alpha[state],))
    assert all(not value or min(value) >= targets[variables]
               for variables, value in tensor.items())
    return {'nonzero_literal_entries': sum(bool(value) for value in tensor.values()),
            'verified_regular_entry_targets': len(targets)}


def synthesize(prime, n, ramification, modes, geometry, common_zero_shift,
               degrees, balance=Fraction(1, 2)):
    degrees = tuple(degrees)
    balance = Fraction(balance)
    if (len(degrees) != prime ** modes
            or any(not isinstance(bound, int) or bound < 0 for bound in degrees)):
        raise ValueError('Supply one nonnegative integer degree bound per joint occupation state.')
    if not isinstance(common_zero_shift, int) or not 0 <= balance <= 1:
        raise ValueError('Use an integer common shift and a rational balance between zero and one.')
    states, variables, equations = coefficient_system(
        prime, n, ramification, modes, geometry, common_zero_shift, degrees, balance)
    nullspace, rank = sparse_nullspace(equations, len(variables), prime)
    witness, dimension = endpoint_witness(nullspace, variables, states, prime)
    result = {
        'prime': prime, 'n': n, 'ramification': ramification, 'modes': modes,
        'geometry': geometry, 'common_zero_shift': common_zero_shift,
        'degree_bounds': degrees, 'balance': str(balance),
        'coefficient_variables': len(variables), 'coefficient_constraints': len(equations),
        'homogeneous_rank': rank, 'homogeneous_nullity': len(nullspace),
        'endpoint_image_dimension': dimension, 'feasible': witness is not None,
    }
    if witness is None:
        return result
    alpha = {state: {} for state in states}
    for column, coefficient in witness.items():
        state, degree = variables[column]
        alpha[states[state]][degree] = coefficient
    result['witness_verification'] = verify_witness(
        prime, n, ramification, modes, geometry, common_zero_shift, degrees, balance, alpha)
    result['actual_degrees'] = tuple(max(alpha[state]) for state in states)
    result['alpha'] = alpha
    return result


def json_ready(result):
    result = dict(result)
    if 'alpha' in result:
        result['alpha'] = [{'state': state, 'polynomial': sorted(value.items())}
                           for state, value in result['alpha'].items()]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=2)
    parser.add_argument('--n', type=int, default=8)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--modes', type=int, default=2)
    parser.add_argument('--geometry', choices=('initial', 'terminal', 'mixed'), default='mixed')
    parser.add_argument('--zero-shift', type=int, default=1)
    parser.add_argument('--degrees', default='1,1,1,1')
    parser.add_argument('--balance', type=Fraction, default=Fraction(1, 2))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = synthesize(args.prime, args.n, args.ramification, args.modes, args.geometry,
                        args.zero_shift, tuple(map(int, args.degrees.split(','))), args.balance)
    encoded = json.dumps(json_ready(result), indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
