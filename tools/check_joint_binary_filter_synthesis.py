#!/usr/bin/env python3
"""Compare compressed joint synthesis with full-coefficient exact equations."""
from itertools import product
import json
from pathlib import Path
import random

from joint_binary_filter_synthesis import synthesize
from weyl_lattice_filter_synthesis import polynomial_multiply, reduce_row


def independent_matrix(modes):
    states = tuple(product(range(2), repeat=modes))
    one = ((0, 1), (1, 3))  # [0,1;1,1+X] in F2[X].
    result = []
    for row in states:
        values = []
        for column in states:
            value = 1
            for i, j in zip(row, column):
                value = polynomial_multiply(value, one[i][j])
            values.append(value)
        result.append(tuple(values))
    return tuple(result)


def full_oracle(modes, m, degrees, orders):
    transform = independent_matrix(modes)
    variables = tuple((state, degree) for state, bound in enumerate(degrees)
                      for degree in range(bound + 1))
    equations = {}
    endpoints = [0] * (2 ** modes)
    for column, (state, degree) in enumerate(variables):
        bit = 1 << column
        endpoints[state] |= bit
        for row, values in enumerate(transform):
            polynomial = values[state]
            while polynomial:
                low = polynomial & -polynomial
                exponent = degree + m * (low.bit_length() - 1)
                polynomial ^= low
                if exponent < orders[row]:
                    equations[row, exponent] = equations.get((row, exponent), 0) ^ bit
    pivots = {}
    for equation in equations.values():
        while equation:
            pivot = equation.bit_length() - 1
            if pivot in pivots:
                equation ^= pivots[pivot]
            else:
                pivots[pivot] = equation
                break
    return all(reduce_row(endpoint, pivots) for endpoint in endpoints)


def main():
    rng = random.Random(0xB37A2026)
    comparisons = witnesses = 0
    examples = []
    for modes, m in product((2, 3, 4), (1, 2, 7, 32)):
        size = 2 ** modes
        transform = independent_matrix(modes)
        for trial in range(24):
            degrees = tuple(rng.choice((0, 1, m, m + 1, 2 * m)) for _ in range(size))
            if trial % 3:
                orders = tuple(max(0, m * k.bit_count() + rng.randrange(-2, 3)) for k in range(size))
            else:
                orders = tuple(rng.randrange(m * modes + 3) for _ in range(size))
            independent = full_oracle(modes, m, degrees, orders)
            result = synthesize(modes, m, degrees, orders, construct=False)
            assert bool(result) == independent
            comparisons += 1
            if not result:
                continue
            witness = synthesize(modes, m, degrees, orders)
            alpha = tuple(witness['alpha'].values())
            beta = []
            for row in transform:
                value = {}
                for polynomial, coefficient in zip(alpha, row):
                    while coefficient:
                        bit = coefficient & -coefficient
                        shift = m * (bit.bit_length() - 1)
                        coefficient ^= bit
                        for degree, scalar in polynomial.items():
                            value[degree + shift] = value.get(degree + shift, 0) ^ scalar
                beta.append({degree: scalar for degree, scalar in value.items() if scalar})
            assert tuple(beta) == tuple(witness['beta'].values())
            assert all(max(value) <= bound for value, bound in zip(alpha, degrees))
            assert all(min(value) >= bound for value, bound in zip(beta, orders))
            witnesses += 1
            if len(examples) < 12:
                examples.append({key: value for key, value in witness.items() if key not in ('alpha', 'beta')})
    result = {'status': 'All complete-coefficient versus residue-compressed controls passed.',
              'exact_feasibility_comparisons': comparisons,
              'literal_polynomial_witness_replays': witnesses,
              'modes_tested': [2, 3, 4], 'residue_scales_tested': [1, 2, 7, 32],
              'sample_witnesses': examples}
    path = Path(__file__).resolve().parents[1] / 'evidence' / 'joint-binary-filter-synthesis.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'sample_witnesses'}))


if __name__ == '__main__':
    main()
