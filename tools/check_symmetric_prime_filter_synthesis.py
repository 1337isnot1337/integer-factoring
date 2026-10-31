#!/usr/bin/env python3
"""Replay symmetric prime filters in the original digit-state transform."""

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
import random

from joint_prime_filter_synthesis import extension_add, extension_scale, transform as original_transform
from prime_power_identity_filters import finite_difference_terms
from symmetric_prime_filter_synthesis import synthesize, transform


def endpoint(polynomial, prime):
    result = 0
    for value in polynomial.values():
        result = extension_add(result, value, prime)
    return result


def original_symmetric_transform(prime, modes, m, values):
    last = modes * (prime - 1)
    alpha = {state: values[last - sum(state)] for state in product(range(prime), repeat=modes)}
    return original_transform(prime, modes, m, alpha)


def uncompressed_feasible(prime, modes, m, degrees, orders):
    last = modes * (prime - 1)
    variables = tuple((weight, degree) for weight, bound in enumerate(degrees)
                      for degree in range(bound + 1))
    equations = {}
    terms = finite_difference_terms(prime, modes)
    for column, (weight, degree) in enumerate(variables):
        for state, row in terms.items():
            for source, power, scalar in row:
                if last - sum(source) != weight:
                    continue
                for shift in range(power + 1):
                    coefficient = scalar * comb(power, shift) * (-1) ** shift % prime
                    at = degree + m * shift
                    if coefficient and at < orders[sum(state)]:
                        key = state, at
                        equations.setdefault(key, {})
                        equations[key][column] = (equations[key].get(column, 0) + coefficient) % prime
    pivots = {}

    def reduce(vector):
        vector = {key: value % prime for key, value in vector.items() if value % prime}
        while vector:
            pivot = max(vector)
            if pivot not in pivots:
                break
            coefficient = vector[pivot]
            for column, value in pivots[pivot].items():
                vector[column] = (vector.get(column, 0) - coefficient * value) % prime
                if not vector[column]:
                    del vector[column]
        return vector

    for row in equations.values():
        row = reduce(row)
        if row:
            pivot = max(row)
            inverse = pow(row[pivot], -1, prime)
            pivots[pivot] = {column: value * inverse % prime for column, value in row.items()}
    for weight in range(last + 1):
        vector = {column: 1 for column, (source, degree) in enumerate(variables) if source == weight}
        # Unlike insertion, row-space membership must skip absent pivots
        # while continuing to reduce all lower pivot coordinates.
        for pivot in sorted(pivots, reverse=True):
            coefficient = vector.get(pivot, 0)
            if coefficient:
                for column, value in pivots[pivot].items():
                    vector[column] = (vector.get(column, 0) - coefficient * value) % prime
                    if not vector[column]:
                        del vector[column]
        if not vector:
            return False
    return True


def prime_product(first, second, prime):
    result = {}
    for a, c in first.items():
        for b, d in second.items():
            result[a + b] = (result.get(a + b, 0) + c * d) % prime
    return {degree: value for degree, value in result.items() if value}


def exact_ratio(prime, n, energies, beta, values, base, charge):
    positive = negative = normal_zero = normal_infinity = Fraction(0)
    for state in product(range(prime), repeat=len(energies)):
        weight = sum(state)
        energy = sum(a * b for a, b in zip(energies, state))
        zero, infinity = base ** (energy + charge * weight), base ** energy
        normal_zero += zero
        normal_infinity += infinity
        positive += zero * base ** (2 * (min(beta[state]) - n * weight))
        negative += infinity * base ** (-2 * max(values[state]))
    return positive * negative / (normal_zero * normal_infinity)


def tensor_controls():
    prime, n, m = 3, 3, 3
    first = ({0: 1, 3: 1, 6: 2}, {0: 1}, {0: 1, 1: 1, 4: 2})
    second = ({0: 1}, {0: 1, 3: 2, 7: 1}, {0: 1})
    beta_first = original_symmetric_transform(prime, 1, m, first)
    beta_second = original_symmetric_transform(prime, 1, m, second)
    values = {(a, b): prime_product(first[a], second[b], prime)
              for a, b in product(range(prime), repeat=2)}
    alpha = {tuple(prime - 1 - digit for digit in state): polynomial
             for state, polynomial in values.items()}
    beta = original_transform(prime, 2, m, alpha)
    for a, b in values:
        assert beta[a, b] == prime_product(beta_first[a,], beta_second[b,], prime)
        assert max(values[a, b]) == max(first[a]) + max(second[b])
        assert min(beta[a, b]) == min(beta_first[a,]) + min(beta_second[b,])
        assert endpoint(values[a, b], prime) == 1
        assert endpoint(beta[a, b], prime) == (-1) ** (a + b) % prime
    count = 0
    for base in (Fraction(4, 5), Fraction(49, 50)):
        for charge in (-6, 0, 6):
            left = exact_ratio(prime, n, (5,), beta_first,
                               {(a,): first[a] for a in range(prime)}, base, charge)
            right = exact_ratio(prime, n, (1,), beta_second,
                                {(a,): second[a] for a in range(prime)}, base, charge)
            combined = exact_ratio(prime, n, (5, 1), beta, values, base, charge)
            assert combined == left * right
            count += 1
    return {'prime': prime, 'different_block_full_states': 9,
            'unequal_physical_mode_energies': [5, 1], 'exact_cost_factorizations': count}


def main():
    directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=directory.parent / 'evidence' /
                        'symmetric-prime-filter-synthesis.json')
    args = parser.parse_args()
    generator = random.Random(720493)
    baseline = originals = feasibility = witness_states = 0
    for prime, modes in ((2, 3), (3, 1), (3, 2), (3, 4), (5, 1), (5, 2)):
        last = modes * (prime - 1)
        for m in (1, 3, 9):
            values = ({0: 1},) * (last + 1)
            beta = transform(prime, modes, m, values)
            assert beta == tuple({m * weight: (-1) ** weight % prime} for weight in range(last + 1))
            assert all(polynomial == beta[sum(state)] for state, polynomial in
                       original_symmetric_transform(prime, modes, m, values).items())
            baseline += 1
        for case in range(5):
            m = 1 + case % 3
            values = []
            for weight in range(last + 1):
                polynomial = {generator.randrange(1, 7): generator.randrange(1, prime ** 2)
                              for _ in range(3)}
                polynomial[0] = extension_add(1, extension_scale(endpoint(polynomial, prime), -1, prime), prime)
                values.append({degree: value for degree, value in polynomial.items() if value})
            beta = transform(prime, modes, m, values)
            original = original_symmetric_transform(prime, modes, m, values)
            assert all(polynomial == beta[sum(state)] for state, polynomial in original.items())
            assert all(endpoint(value, prime) == (-1) ** weight % prime
                       for weight, value in enumerate(beta))
            originals += len(original)
        for case in range(8):
            m = 1 + case % 3
            degrees = tuple(generator.randrange(5) for _ in range(last + 1))
            orders = tuple(generator.randrange(2 * m + 5) for _ in range(last + 1))
            expected = uncompressed_feasible(prime, modes, m, degrees, orders)
            witness = synthesize(prime, modes, m, degrees, orders)
            assert bool(witness) == expected
            feasibility += 1
            if witness:
                original = original_symmetric_transform(prime, modes, m, witness['f'])
                assert all(polynomial == witness['beta'][sum(state)] for state, polynomial in original.items())
                witness_states += len(original)
    sources = ('symmetric_prime_filter_synthesis.py', Path(__file__).name,
               'joint_prime_filter_synthesis.py', 'prime_power_identity_filters.py',
               'polynomial_state_filter_synthesis.py')
    record = {
        'status': 'Prime-weight compression agrees with the original digit transform, endpoints, and exact feasibility.',
        'scope': 'Constructive prime-state symmetric multilevel filters; no factoring exponent improvement.',
        'exact_baseline_and_endpoint_sign_controls': baseline,
        'nonconstant_original_digit_state_controls': originals,
        'uncompressed_feasibility_comparisons': feasibility,
        'synthesized_original_state_controls': witness_states,
        'different_block_tensor_controls': tensor_controls(),
        'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in sources},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
