#!/usr/bin/env python3
"""Exact four-copy criterion and eighth-power gauges for binary typed endpoints."""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import gcd

from check_rectangular_defect_support import pdivmod
from weyl_lattice_filter_synthesis import polynomial_multiply, polynomial_degrees


DIRECTIONS = tuple(a for a in product((-1, 0, 1), repeat=2) if a != (0, 0))


def determinant(a, b):
    return a[0]*b[1]-a[1]*b[0]


def binary_charge_circuits():
    result = []
    for indices in combinations(range(len(DIRECTIONS)), 2):
        a, b = (DIRECTIONS[i] for i in indices)
        if tuple(x+y for x, y in zip(a, b)) == (0, 0):
            result.append((indices, (1, 1)))
    for indices in combinations(range(len(DIRECTIONS)), 3):
        a, b, c = (DIRECTIONS[i] for i in indices)
        weights = [determinant(b, c), determinant(c, a), determinant(a, b)]
        if all(value < 0 for value in weights):
            weights = [-value for value in weights]
        if not all(value > 0 for value in weights):
            continue
        divisor = gcd(gcd(weights[0], weights[1]), weights[2])
        result.append((indices, tuple(value//divisor for value in weights)))
    return tuple(result)


CIRCUITS = binary_charge_circuits()


def circuit_certificate(direction_valuations):
    values = dict(direction_valuations)
    if any(direction not in DIRECTIONS or not isinstance(value, int) for direction, value in values.items()):
        raise ValueError('Use integer valuations on nonzero binary paired-charge directions.')
    rows = []
    for indices, weights in CIRCUITS:
        directions = tuple(DIRECTIONS[i] for i in indices)
        if all(direction in values for direction in directions):
            total = sum(weight*values[direction] for weight, direction in zip(weights, directions))
            rows.append({'directions': directions, 'multiplicities': weights,
                         'word_length': sum(weights), 'valuation': total})
    return {'valid': all(row['valuation'] > 0 for row in rows), 'circuits': rows}


def eighth_power_potential(direction_valuations):
    """Return integer paired potentials after w -> z^8, or an invalid circuit.

    Valid circuit means are at least1/4. A feasible quarter-margin polygon
    has an eighth-integral vertex after intersection with a large integer box.
    Box sizes increase until such a vertex is found; exact arithmetic only.
    """
    values = dict(direction_valuations)
    certificate = circuit_certificate(values)
    if not certificate['valid']:
        return {'valid': False, 'failed_circuit': next(row for row in certificate['circuits'] if row['valuation'] <= 0)}
    constraints = [(a, b, Fraction(1, 4)-value) for (a, b), value in values.items()]
    bound = 1
    while True:
        bounded = constraints+[(1, 0, -bound), (-1, 0, -bound), (0, 1, -bound), (0, -1, -bound)]
        candidates = []
        for (a, b, c), (d, e, f) in combinations(bounded, 2):
            det = a*e-b*d
            if not det:
                continue
            x = Fraction(c*e-b*f, det); y = Fraction(a*f-c*d, det)
            if all(u*x+v*y >= rhs for u, v, rhs in bounded):
                candidates.append((x, y))
        if candidates:
            x, y = min(candidates, key=lambda row: (abs(row[0])+abs(row[1]), row))
            assert (8*x).denominator == (8*y).denominator == 1
            potential = int(8*x), int(8*y)
            margins = [8*value+a*potential[0]+b*potential[1] for (a, b), value in values.items()]
            assert all(value >= 2 for value in margins)
            return {'valid': True, 'ramification': 8, 'integer_potential': potential,
                    'minimum_transformed_margin': min(margins, default=None),
                    'box_bound': bound, 'circuits': certificate['circuits']}
        bound *= 2


def chi_power(power):
    if power < 0:
        raise ValueError('Use a nonnegative endpoint polynomial exponent.')
    result = 0; subset = power
    while True:
        result |= 1 << subset
        if subset == 0:
            return result
        subset = (subset-1) & power


def chi_order(polynomial):
    if not polynomial:
        return None
    order = 0
    while True:
        quotient, remainder = pdivmod(polynomial, 3)
        if remainder:
            return order
        polynomial = quotient; order += 1


def labels(entry):
    return tuple((entry >> bit) & 1 for bit in range(4))


def numerator_direction_valuations(values, pole):
    if len(values) != 16 or not isinstance(pole, int) or pole < 0:
        raise ValueError('Use sixteen F4 polynomial numerators and one common endpoint pole.')
    result = {}
    for entry, pair in enumerate(values):
        x, y, u, v = labels(entry)
        direction = x-v, y-u
        orders = [chi_order(value) for value in pair if value]
        if direction == (0, 0):
            evaluated = []
            for value in pair:
                quotient, remainder = pdivmod(value, chi_power(pole))
                if remainder:
                    raise ValueError('A matched coefficient has an endpoint pole.')
                evaluated.append(quotient.bit_count() & 1)
            if evaluated != [1, 0]:
                raise ValueError('Every matched coefficient must have endpoint one.')
        elif orders:
            valuation = min(orders)-pole
            result[direction] = min(result.get(direction, valuation), valuation)
    return result


def strict_laurent_gauge(values, pole):
    """Produce a literal strict polynomial numerator in the w-shift-zero case.

    A common Laurent w power can be carried separately and is multiplied
    by eight. The endpoint denominator is removed exactly, coefficient by
    coefficient. Field coefficients are unchanged by the substitution.
    """
    valuations = numerator_direction_valuations(values, pole)
    certificate = eighth_power_potential(valuations)
    if not certificate['valid']:
        return certificate
    f, g = certificate['integer_potential']
    output = []
    for entry, pair in enumerate(values):
        x, y, u, v = labels(entry)
        gauge = f*(x-v)+g*(y-u)
        remaining = 8*pole-gauge
        new_pair = []
        for polynomial in pair:
            substituted = sum(1 << (8*power) for power in polynomial_degrees(polynomial))
            if remaining >= 0:
                transformed, remainder = pdivmod(substituted, chi_power(remaining))
                assert remainder == 0
            else:
                transformed = polynomial_multiply(substituted, chi_power(-remaining))
            new_pair.append(transformed)
        assert ((new_pair[0].bit_count() & 1) | ((new_pair[1].bit_count() & 1) << 1)) == int(x == v and y == u)
        output.append(tuple(new_pair))
    return certificate | {'values': tuple(output), 'degree': max(max(pair).bit_length()-1 for pair in output)}
