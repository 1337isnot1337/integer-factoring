#!/usr/bin/env python3
"""Check the all-prime one-hole block formula and exact minimum tradeoffs."""
from __future__ import annotations

import argparse
from itertools import product
import json
from math import factorial
from pathlib import Path

from coupled_identity_operator_filters import identity_kernel, states
from coupled_identity_weyl_profiles import reduced_matrix


def universal_entry(prime, modes, row, column, source):
    """Apply the universal Weyl differential operator; return target,c^-degree,coefficient."""
    top = (prime-1,)*modes
    power = tuple(t-b for t, b in zip(top, column))
    derivative = tuple(t-a for t, a in zip(top, row))
    multiplied = tuple(s+b for s, b in zip(source, power))
    if any(value >= prime or value < order for value, order in zip(multiplied, derivative)):
        return None
    coefficient = (-1)**(modes+sum(power))
    for a, value, order in zip(row, multiplied, derivative):
        coefficient *= factorial(a)*factorial(value)//factorial(value-order)
    target = tuple(value-order for value, order in zip(multiplied, derivative))
    return target, sum(power), coefficient % prime


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    symbolic, lower = 0, 0
    for prime, modes in product((2, 3, 5, 7, 11), range(1, 5)):
        top = (prime-1,)*modes
        holes = tuple(tuple(value-int(i == j) for i, value in enumerate(top)) for j in range(modes))
        for source in holes:
            assert universal_entry(prime, modes, top, top, source) == (source, 0, 1)
            symbolic += 1
        for i, row in enumerate(holes):
            for j, column in enumerate(holes):
                for k, source in enumerate(holes):
                    expected = (holes[i], 1, prime-1) if j == k else None
                    assert universal_entry(prime, modes, row, column, source) == expected
                    symbolic += 1
        deficits = [value for value in product(range(min(3, prime)), repeat=modes) if sum(value) == 2]
        labels = [tuple(t-d for t, d in zip(top, deficit)) for deficit in deficits]
        for row, column, source in product(labels, labels, holes):
            assert universal_entry(prime, modes, row, column, source) is None
            lower += 1
    profiles = []
    for prime, n, P, selected in ((2, 4, 1, (0, 3)), (2, 8, 3, (1, 4, 7)),
                                   (3, 3, 1, (0, 2)), (5, 5, 1, (2,))):
        modes, geometry = len(selected), 'selected:'+','.join(map(str, selected))
        top = (prime-1,)*modes
        holes = tuple(tuple(value-int(i == j) for i, value in enumerate(top)) for j in range(modes))
        weights = tuple(P*sum(j*digit for j, digit in zip(selected, state)) for state in holes)
        kappas = tuple(P*sum((2*j+1-2*n)*digit for j, digit in zip(selected, state)) for state in holes)
        for gain in (1, 2, 3):
            value = identity_kernel(prime, modes)
            for row, kappa in zip(holes, kappas):
                degree = (kappa-min(kappas))//2-gain
                value[row, row] = {degree: 1}
            for i, row in enumerate(holes):
                for j, column in enumerate(holes):
                    if i != j:
                        degree = (kappas[j]-min(kappas))//2-gain
                        value[row, column] = {degree: 1, degree-1: prime-1}
            infinity_minimum = min(kappa-2*max(value[row, column])
                                   for row in holes for column, kappa in zip(holes, kappas)
                                   if value.get((row, column)))
            assert infinity_minimum == min(kappas)+2*gain
            matrix = reduced_matrix(prime, n, P, modes, value, geometry)
            alphabet = states(prime, modes)
            indices = [alphabet.index(state) for state in holes]
            zero_minimum = min(degree for i in indices for j in indices for degree in matrix[i][j])
            assert zero_minimum <= min(weights)-gain
            assert min(weights)+n*P-zero_minimum >= n*P+gain
            profiles.append({'prime': prime, 'n': n, 'ramification': P, 'selected': selected,
                             'infinity_gain': 2*gain, 'zero_q_minimum': zero_minimum,
                             'Gaussian_zero_q_minimum': min(weights)+n*P})
    result = {'status': 'Exact all-prime highest-block identities and minimum tradeoff controls; not a bound for every operator objective.',
              'top_and_one_hole_matrix_unit_actions': symbolic,
              'lower_occupation_vanishing_actions': lower,
              'exact_polynomial_controls': profiles,
              'strict_positive_gain_required': True}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'top_and_one_hole_matrix_unit_actions': symbolic,
                      'lower_occupation_vanishing_actions': lower,
                      'exact_polynomial_controls': len(profiles)}))


if __name__ == '__main__':
    main()
