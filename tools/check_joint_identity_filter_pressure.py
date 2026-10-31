#!/usr/bin/env python3
"""Replay fast all-index profile and pressure formulas against complete sums."""
from fractions import Fraction
from itertools import combinations, permutations, product
import json
from math import log, sqrt
from pathlib import Path
import random

import numpy as np
from scipy.special import logsumexp

from joint_identity_filter_pressure import compiled_binary_beta_orders, index_gradient, pressure, relative_log_cost
from prime_power_identity_filters import beta_orders, profile_data
from search_all_index_joint_filters import subset_weights
from search_f4_repaired_identity_filters import candidate_exponents


def main():
    rng = random.Random(0xA11B37A)
    controls = gradient_controls = 0
    discrepancy = 0.0
    for n, modes, field in ((4, 2, 4), (8, 3, 16), (8, 4, 32)):
        states = tuple(product(range(2), repeat=modes))
        for trial in range(20):
            P = (1, 3)[trial % 2]
            selected = tuple(sorted(rng.sample(range(n), modes)))
            exponents = tuple(rng.randrange(-n * P, n * P + 1) for state in states)
            colors = tuple(rng.randrange(1, field) for state in states)
            alpha = {state: {-exponent: color} for state, exponent, color in zip(states, exponents, colors)}
            orders = compiled_binary_beta_orders(n, P, modes, exponents)(colors)
            assert orders == beta_orders(2, n, P, alpha, field)
            degrees = tuple(-value for value in exponents)
            J = Fraction(rng.randrange(-3, 4), 2)
            t = rng.uniform(.03, .3)
            data = profile_data(2, n, P, J, selected, alpha, field)
            ell = t / P
            complete = 0.0
            for cut in ('zero', 'infinity'):
                complete += float(logsumexp(-ell * np.asarray(data[cut + '_active'], dtype=float)))
                complete += sum(float(np.logaddexp(0, -ell * value)) for value in data[cut + '_spectators'])
            physical = data['paired_physical_first'] + data['paired_physical_second']
            complete -= sqrt(2) * sum(float(np.logaddexp(0, -ell * value / sqrt(2))) for value in physical)
            fast = pressure(2, n, P, float(J), t, selected, orders, degrees)
            discrepancy = max(discrepancy, abs(complete - fast))
            assert abs(complete - fast) < 1e-10
            wz, wi = subset_weights(n, modes, (selected,), t, float(J))
            occupations = np.asarray([sum(state) for state in states])
            az = 2 * t * (n * occupations - np.asarray(orders) / P)
            ai = 2 * t * np.asarray(degrees[::-1]) / P
            batched = log(float(wz[0] @ np.exp(az - az.max()))) + float(az.max())
            batched += log(float(wi[0] @ np.exp(ai - ai.max()))) + float(ai.max())
            assert abs(batched - relative_log_cost(2, n, P, float(J), t, selected, orders, degrees)) < 1e-10
            value, derivative = index_gradient(2, n, P, float(J), t, selected, orders, degrees)
            for index in range(modes):
                delta = np.zeros(modes); delta[index] = 1e-5
                finite = (relative_log_cost(2, n, P, float(J), t, np.asarray(selected) + delta, orders, degrees)
                          - relative_log_cost(2, n, P, float(J), t, np.asarray(selected) - delta, orders, degrees)) / 2e-5
                assert abs(finite - derivative[index]) < 1e-7
                gradient_controls += 1
            controls += 1
    closure = []
    states = tuple(product(range(2), repeat=3))
    lookup = {state: index for index, state in enumerate(states)}
    for n in (16, 32):
        family = set(candidate_exponents(n, 1, 3))
        for permutation in permutations(range(3)):
            indices = tuple(lookup[tuple(state[j] for j in permutation)] for state in states)
            assert {tuple(exponents[j] for j in indices) for exponents in family} == family
        closure.append({'n': n, 'exponent_tables': len(family), 'complete_mode_permutations': 6})
    result = {'status': 'Exact compiled beta and full partition controls passed.',
              'complete_partition_comparisons': controls,
              'compiled_beta_comparisons': controls, 'index_derivative_controls': gradient_controls,
              'maximum_pressure_discrepancy': discrepancy,
              'three_mode_exponent_family_permutation_closure': closure}
    path = Path(__file__).resolve().parents[1] / 'evidence' / 'joint-identity-filter-pressure.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
