#!/usr/bin/env python3
"""Exact first-jet certificates and independent arbitrary-rank frame controls."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random

import numpy as np

from check_reflection_laurent_frames import exact_orders
from f16_identity_operator_profiles import inverse_matrix
from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import identity_split, transform_leg
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize, bell_restore, coupled_terms, digest
from square_zero_coupled_reflection_frames import (
    construct, first_jet_identity, matching_generator, matrix_terms,
    product, rank, square_zero_matrix, transpose,
)


def conjugated(generator, steps, rng):
    q = len(generator)
    change = identity_split(q)
    for _ in range(steps):
        i, j = rng.sample(range(q), 2)
        scalar = rng.choice((1, 2, 3))
        from physical_split_crossings import F4_PRODUCT
        change[i] = [a ^ F4_PRODUCT[scalar][b] for a, b in zip(change[i], change[j])]
    inverse = [list(row) for row in inverse_matrix(change)]
    return product(product(change, generator), inverse)


def sequential(seed, generator):
    """Independent dense physical-leg contractions, with one common shift."""
    q = seed['q_x']
    current = [(a << 4, b << 4) for a, b in seed['values']]
    transposed = transpose(generator)
    for leg, matrix in enumerate((generator, transposed, generator, transposed)):
        changed, dimensions = transform_leg(current, [q] * 4, leg, matrix)
        assert dimensions == [q] * 4
        assert not any((a | b) & 1 for a, b in changed)
        current = [(a ^ (c >> 1), b ^ (d >> 1)) for (a, b), (c, d) in zip(current, changed)]
    assert not any((a | b) & 15 for a, b in current)
    return [(a >> 4, b >> 4) for a, b in current]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(312131)
    cases, profile_controls = [], []
    for interval in (4, 8):
        base = tensor(interval)
        data = bell_normalize(base)
        bell, factors = data
        q = base['q_x']
        jet = first_jet_identity(base, data)
        assert jet['scalar'] == 0
        assert jet['operator'] == [[int(i == j) * (i.bit_count() % 2) for j in range(q)] for i in range(q)]
        residue_zero = sum(1 << exponent for exponent in range(0, base['degree'] + 1, 4))
        for index, value in enumerate(base['values']):
            x, v = index % q, index // (q ** 3)
            parity = (x.bit_count() + v.bit_count()) % 2
            for component in value:
                selected = component & residue_zero
                assert component == (selected ^ (selected << 2) if parity else selected)
        generators = []
        for mask in ((1, 2, 3) if interval == 4 else (1, q - 1)):
            generators.append(('xor-matching-' + str(mask), matching_generator([i ^ mask for i in range(q)])))
        permutation = list(range(q))
        for i in range(0, q // 2, 2):
            permutation[i], permutation[i + 1] = permutation[i + 1], permutation[i]
        partial = matching_generator(permutation)
        generators.append(('F4-conjugated-partial-matching', conjugated(partial, q + 3, rng)))
        records = []
        for index, (name, generator) in enumerate(generators):
            assert square_zero_matrix(generator)
            observed = construct(base, generator, restore=False, bell_data=data)
            assert observed['values'] == sequential(bell, generator)
            inverse_terms = matrix_terms(observed, generator)
            from square_zero_coupled_reflection_frames import from_terms
            assert from_terms(observed, inverse_terms)['values'] == bell['values']
            restored = bell_restore(observed, factors)
            original_generator = product(product(factors['first_factor'], generator), factors['inverse_factor'])
            independently_converted = construct(base, original_generator, coordinates='original', bell_data=data)
            assert independently_converted['values'] == restored['values']
            record = {'name': name, 'generator': generator, 'rank': rank(generator),
                      'degree': restored['degree'], 'content': observed['monomial_content'],
                      'matrix_sha256': digest(restored['values']),
                      'independent_coefficient_entries': q ** 4,
                      'involution_entries': q ** 4, 'coordinate_conversion_entries': q ** 4,
                      'endpoint_entries': q ** 4}
            records.append(record)
            if interval == 4 and index in (0, len(generators) - 1):
                for chosen_cut in (0, 1):
                    native = RectangularOracle(restored['values'], q, q)
                    heights = [[[Fraction(0)] + [Fraction(rng.randrange(-4, 5), 3) for _ in range(q - 1)]
                                for _ in range(4)] for _ in range(2)]
                    scale = Fraction(5, 4)
                    parameters = [float(value) for endpoint in heights for leg in endpoint for value in leg[1:]] + [float(scale)]
                    actual = np.asarray(native.prefixes(parameters, chosen_cut)['orders']).reshape(2, q * q)
                    maximum_error = 0.0
                    for side, cut, point in ((0, chosen_cut, None), (1, 1 - chosen_cut, (0, 1))):
                        exact = exact_orders(restored, cut, point, heights[side], scale)
                        if point is None:
                            exact = [value - restored['degree'] * scale for value in exact]
                        error = max(abs(float(value) - observed) for value, observed in zip(exact, actual[side]))
                        assert error < 2e-8
                        maximum_error = max(maximum_error, error)
                    profile_controls.append({'generator': name, 'cut_at_infinity': chosen_cut,
                                             'complete_profiles': 2, 'maximum_error': maximum_error})
        coordinate = [[0] * q for _ in range(q)]
        coordinate[0][q - 1] = 2
        assert matrix_terms(bell, coordinate) == coupled_terms(bell, 0, q - 1, 2)
        invalid_generator_rejected = False
        try:
            construct(base, identity_split(q), bell_data=data)
        except ValueError:
            invalid_generator_rejected = True
        assert invalid_generator_rejected
        altered_values = list(base['values'])
        a, b = altered_values[0]
        altered_values[0] = a ^ 6, b  # w+w² is zero at the crossing endpoint.
        invalid_jet_rejected = False
        try:
            first_jet_identity(dict(base, values=altered_values))
        except ArithmeticError:
            invalid_jet_rejected = True
        assert invalid_jet_rejected
        cases.append({'interval': interval, 'alphabet': q, 'degree': base['degree'],
                      'first_jet': jet, 'generators': records,
                      'full_polynomial_mod_four_entries': q ** 4,
                      'coordinate_rank_one_recovery': True,
                      'nonnilpotent_generator_rejected': invalid_generator_rejected,
                      'nonzero_linear_jet_rejected': invalid_jet_rejected})
    directory = Path(__file__).parent
    names = (Path(__file__).name, 'square_zero_coupled_reflection_frames.py',
             'reflection_laurent_physical_frames.py', 'reflection_dual_rail_gaussian.py',
             'physical_split_crossings.py', 'check_reflection_laurent_frames.py',
             'f16_weighted_local_oracle.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp')
    result = {'scope': 'Exact first-jet hypotheses plus arbitrary-rank coupled-frame coefficient and profile controls; no pressure improvement certified.',
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names},
              'cases': cases, 'profile_controls': profile_controls,
              'totals': {'first_jet_entries': sum(case['first_jet']['exact_coefficient_entries'] for case in cases),
                         'full_polynomial_mod_four_entries': sum(case['full_polynomial_mod_four_entries'] for case in cases),
                         'independent_coefficient_entries': sum(sum(row['independent_coefficient_entries'] for row in case['generators']) for case in cases),
                         'coordinate_conversion_entries': sum(sum(row['coordinate_conversion_entries'] for row in case['generators']) for case in cases),
                         'complete_exact_profiles': sum(row['complete_profiles'] for row in profile_controls)}}
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(json.dumps(result['totals']))


if __name__ == '__main__':
    main()
