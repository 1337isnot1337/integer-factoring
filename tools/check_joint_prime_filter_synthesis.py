#!/usr/bin/env python3
"""Independent full-coefficient and componentwise checks for joint prime filters."""
from itertools import product
import json
from pathlib import Path
import random

from joint_prime_filter_synthesis import synthesize
from polynomial_state_filter_synthesis import _rref, difference_matrix as single_matrix, transform as single_transform


def independent_matrix(prime, modes):
    states = tuple(product(range(prime), repeat=modes)); one = single_matrix(prime)
    result = []
    for row in states:
        values = []
        for column in states:
            value = {0: 1}
            for i, j in zip(row, column):
                updated = {}
                for a, x in value.items():
                    for b, y in one[i][j].items():
                        updated[a + b] = (updated.get(a + b, 0) + x * y) % prime
                value = {degree: coefficient for degree, coefficient in updated.items() if coefficient}
            values.append(value)
        result.append(tuple(values))
    return tuple(result)


def full_oracle(prime, modes, m, degrees, orders):
    matrix = independent_matrix(prime, modes)
    variables = tuple((state, degree) for state, bound in enumerate(degrees) for degree in range(bound + 1))
    constraints = []
    for row, order in enumerate(orders):
        maximum = max((degree + m * power for state, degree in variables for power in matrix[row][state]), default=-1)
        for at in range(min(order, maximum + 1)):
            constraints.append([matrix[row][state].get((at - degree) // m, 0)
                                if (at - degree) % m == 0 else 0 for state, degree in variables])
    reduced, pivots = _rref(constraints, len(variables), prime)
    support = [False] * prime ** modes
    for free in range(len(variables)):
        if free in pivots:
            continue
        lift = [0] * len(variables); lift[free] = 1
        for row, pivot in zip(reduced, pivots):lift[pivot] = -row[free] % prime
        endpoint = [0] * len(support)
        for value, (state, degree) in zip(lift, variables):endpoint[state] = (endpoint[state] + value) % prime
        support = [old or bool(value) for old, value in zip(support, endpoint)]
    return all(support)


def component_transform(prime, modes, m, alpha):
    states = tuple(product(range(prime), repeat=modes))
    result = {state: {} for state in states}
    for component in range(modes + 1):
        values = {state: {degree: (coefficient // prime ** component) % prime
                          for degree, coefficient in polynomial.items()
                          if (coefficient // prime ** component) % prime}
                  for state, polynomial in alpha.items()}
        for axis in range(modes):
            updated = {}
            for other in product(range(prime), repeat=modes - 1):
                indices = tuple(other[:axis] + (digit,) + other[axis:] for digit in range(prime))
                transformed = single_transform(tuple(values[state] for state in indices), m, prime)
                for state, polynomial in zip(indices, transformed):updated[state] = polynomial
            values = updated
        for state, polynomial in values.items():
            for degree, coefficient in polynomial.items():
                result[state][degree] = result[state].get(degree, 0) + coefficient * prime ** component
    return result


def main():
    rng = random.Random(0xC0F27)
    comparisons = witnesses = 0
    for prime, modes, m in ((3, 1, 7), (3, 2, 1), (3, 2, 2), (3, 2, 7), (5, 1, 7)):
        states = tuple(product(range(prime), repeat=modes))
        for trial in range(20):
            degrees = tuple(rng.choice((0, 1, m, m + 1)) for state in states)
            orders = tuple(max(0, m * sum(state) + rng.randrange(-1, 2)) for state in states)
            expected = full_oracle(prime, modes, m, degrees, orders)
            feasible = synthesize(prime, modes, m, degrees, orders, construct=False)
            assert bool(feasible) == expected
            comparisons += 1
            if feasible:
                witness = synthesize(prime, modes, m, degrees, orders)
                assert witness['beta'] == component_transform(prime, modes, m, witness['alpha'])
                witnesses += 1
    result = {'status': 'All exact joint-prime residue controls passed.',
              'full_coefficient_feasibility_comparisons': comparisons,
              'independent_componentwise_tensor_transform_replays': witnesses,
              'extension_fields': [9, 27, 25]}
    path = Path(__file__).resolve().parents[1] / 'evidence' / 'joint-prime-filter-synthesis.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
