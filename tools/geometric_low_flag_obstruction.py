#!/usr/bin/env python3
"""Construct degree-three identities excluding every rank-one low-flag vector.

Polynomials use ordinary exponent tuples over F2. No Boolean relations or
finite-field equations are imposed, so the identities hold in every
characteristic-two field.
"""

from itertools import combinations_with_replacement
import json

from check_pure_two_copy_flag_floor import annihilator, profile_and_flags


def monomials(variables, degree):
    result = []
    for size in range(degree + 1):
        for indices in combinations_with_replacement(range(variables), size):
            exponents = [0] * variables
            for index in indices:
                exponents[index] += 1
            result.append(tuple(exponents))
    return result


def chart_equations(forms, first_left, first_right):
    variables = [('a', i) for i in range(first_left + 1, 5)] + [('b', i) for i in range(first_right + 1, 5)]
    indices = {value: i for i, value in enumerate(variables)}
    result = []
    for form in forms:
        polynomial = set()
        for i in range(first_left, 5):
            for j in range(first_right, 5):
                if not form >> (5 * i + j) & 1:
                    continue
                monomial = [0] * len(variables)
                if i > first_left:
                    monomial[indices['a', i]] += 1
                if j > first_right:
                    monomial[indices['b', j]] += 1
                term = tuple(monomial)
                polynomial.symmetric_difference_update({term})
        result.append(polynomial)
    return variables, result


def find_identity(forms, first_left, first_right, degree):
    variables, equations = chart_equations(forms, first_left, first_right)
    universe = monomials(len(variables), degree)
    coordinate = {monomial: i for i, monomial in enumerate(universe)}
    pivots, labels = {}, []
    for equation, polynomial in enumerate(equations):
        if not polynomial:
            continue
        maximum = max(map(sum, polynomial))
        if maximum > degree:
            continue
        for multiplier in monomials(len(variables), degree - maximum):
            shifted = [tuple(a + b for a, b in zip(multiplier, term)) for term in polynomial]
            row = sum(1 << coordinate[term] for term in shifted)
            trace = 1 << len(labels)
            labels.append((equation, multiplier))
            while row:
                pivot = row.bit_length() - 1
                if pivot not in pivots:
                    pivots[pivot] = row, trace
                    break
                previous, history = pivots[pivot]
                row ^= previous
                trace ^= history
            if 0 in pivots:
                assert pivots[0][0] == 1
                trace = pivots[0][1]
                multipliers = [[] for _ in equations]
                while trace:
                    bit = trace & -trace
                    trace ^= bit
                    equation, multiplier = labels[bit.bit_length() - 1]
                    multipliers[equation].append(multiplier)
                return {'first_left': first_left, 'first_right': first_right,
                        'variables': variables, 'identity_degree_bound': degree,
                        'equations': [sorted(polynomial) for polynomial in equations],
                        'multipliers': multipliers}
    return None


def construct():
    data = profile_and_flags()
    forms = annihilator(data['common_flag_basis'])
    charts = []
    for first_left in range(5):
        for first_right in range(5):
            for degree in range(4):
                identity = find_identity(forms, first_left, first_right, degree)
                if identity is not None:
                    charts.append(identity)
                    break
            else:
                raise AssertionError('The displayed degree-three geometric certificate was not recovered.')
    return {'construction': 'Exact ordinary-polynomial identities over F2 on every projective rank-one chart.',
            'reference_digest': data['reference_digest'],
            'relative_order_cutoff': 9, 'flag_dimension': 11,
            'physical_image_flag_basis': data['common_flag_basis'],
            'annihilating_forms': forms, 'charts': charts,
            'projective_chart_count': len(charts),
            'maximum_identity_degree': max(chart['identity_degree_bound'] for chart in charts),
            'total_multiplier_monomials': sum(sum(map(len, chart['multipliers'])) for chart in charts),
            'coefficient_scope': 'Every field of characteristic two, including its algebraic closure.',
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(construct(), indent=2))
