#!/usr/bin/env python3
"""Independent Laurent expansion and complete-profile controls for coupled frames."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random

import numpy as np

from check_reflection_laurent_frames import exact_orders
from native_rectangular_smith_f4 import RectangularOracle
from reflection_laurent_physical_frames import (
    bell_normalize, bell_restore, coupled_admissible_poles, coupled_terms,
    coupled_transformed, digest, scale, tensor,
)


def sequential_laurent(seed, target, source, pole, scalar):
    """Independent four physical-leg actions with one common 4k shift."""
    q, shift = seed['q_x'], 4 * pole
    values = [(a << shift, b << shift) for a, b in seed['values']]
    low_mask = (1 << pole) - 1
    for leg, output, input_symbol in ((0, target, source), (1, source, target),
                                      (2, target, source), (3, source, target)):
        previous = values
        values = previous.copy()
        stride = q ** leg
        for index in range(len(previous)):
            if (index // stride) % q != output:
                continue
            source_index = index + (input_symbol - output) * stride
            a, b = scale(previous[source_index], scalar)
            assert not ((a | b) & low_mask)
            before = values[index]
            values[index] = before[0] ^ (a >> pole), before[1] ^ (b >> pole)
    mask = (1 << shift) - 1
    if any((a | b) & mask for a, b in values):
        raise ArithmeticError('The independent complete Laurent expansion retains a pole.')
    return [(a >> shift, b >> shift) for a, b in values]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    generator = random.Random(312128)
    cases, profile_records = [], []
    for interval, maximum in ((4, 8), (8, 4)):
        seed = tensor(interval)
        bell, factors = bell_normalize(seed)
        assert bell_restore(bell, factors)['values'] == seed['values']
        q = seed['q_x']
        catalog = []
        for target in range(q):
            for source in range(q):
                if target == source:
                    continue
                terms = coupled_terms(bell, target, source)
                allowed = coupled_admissible_poles(terms, maximum)
                if allowed:
                    catalog.append({'target': target, 'source': source, 'poles': allowed})
        assert len(catalog) == q * (q - 1)
        assert sum(len(row['poles']) for row in catalog) == (34 if interval == 4 else 688)
        selected = [catalog[i] for i in (0, len(catalog) // 3, 2 * len(catalog) // 3, len(catalog) - 1)]
        records = []
        for row in selected:
            target, source = row['target'], row['source']
            for scalar in (1, 2):
                terms = coupled_terms(bell, target, source, scalar)
                allowed = coupled_admissible_poles(terms, maximum)
                assert allowed
                pole = max(allowed) if scalar == 1 else min(allowed)
                changed = coupled_transformed(bell, target, source, -pole, scalar, terms)
                assert changed['values'] == sequential_laurent(bell, target, source, pole, scalar)
                inverse = coupled_transformed(changed, target, source, -pole, scalar)
                assert inverse['values'] == bell['values']
                restored = bell_restore(changed, factors)
                assert restored['degree'] <= seed['degree']
                record = {'target': target, 'source': source, 'pole': pole, 'scalar': scalar,
                          'degree': restored['degree'], 'content': changed['monomial_content'],
                          'bell_matrix_sha256': digest(changed['values']),
                          'restored_matrix_sha256': digest(restored['values']),
                          'independent_coefficient_entries': q ** 4,
                          'involution_entries': q ** 4, 'endpoint_entries': q ** 4}
                records.append(record)
                if interval == 4:
                    native = RectangularOracle(restored['values'], q, q)
                    height_sets = [[[Fraction(0)] + [Fraction(generator.randrange(-5, 6), 4) for _ in range(q - 1)]
                                    for _ in range(4)] for _ in range(2)]
                    parameter_scale = Fraction(3, 2)
                    params = [float(value) for heights in height_sets for leg in heights for value in leg[1:]]
                    params.append(float(parameter_scale))
                    observed = np.asarray(native.prefixes(params, 1)['orders']).reshape(2, q * q)
                    errors = []
                    for side, cut, point in ((0, 1, None), (1, 0, (0, 1))):
                        exact = exact_orders(restored, cut, point, height_sets[side], parameter_scale)
                        if point is None:
                            exact = [value - restored['degree'] * parameter_scale for value in exact]
                        error = max(abs(float(value) - actual) for value, actual in zip(exact, observed[side]))
                        assert error < 2e-8
                        errors.append(error)
                    profile_records.append({'frame': [target, source, -pole, scalar],
                                            'complete_profiles': 2, 'maximum_profile_error': max(errors)})
        # This value lies beyond the degree and cannot be admissible for a
        # nontrivial coefficient action. Both independent expansions reject it.
        rejected = []
        for function in (lambda: coupled_transformed(bell, 0, 1, -seed['degree'] - 1),
                         lambda: sequential_laurent(bell, 0, 1, seed['degree'] + 1, 1)):
            try:
                function()
            except ArithmeticError:
                rejected.append(True)
            else:
                raise AssertionError('An actual coupled pole was accepted.')
        representative_source = 1 if interval == 4 else 8
        representative = bell_restore(coupled_transformed(bell, 0, representative_source, -1), factors)
        odd_entries = [(index, value) for index, value in enumerate(representative['values'])
                       if index.bit_count() % 2 and value != (0, 0)]
        assert representative['values'][0] != (0, 0) and odd_entries
        first_odd_index, first_odd_value = odd_entries[0]
        cases.append({'interval': interval, 'alphabet': q, 'degree': seed['degree'],
                      'maximum_tested_pole': maximum, 'coefficient_scalar_for_catalogue': 1,
                      'factors': factors, 'catalogue': catalog,
                      'admissible_frames_in_stated_range': sum(len(row['poles']) for row in catalog),
                      'selected_full_coefficient_controls': records,
                      'both_pole_negative_controls_rejected': all(rejected),
                      'bell_normalization_entries': q ** 4, 'bell_restoration_entries': q ** 4,
                      'even_gaussian_locus_escape': {'frame': [0, representative_source, -1, 1],
                                                    'vacuum': representative['values'][0],
                                                    'odd_entry_index': first_odd_index,
                                                    'odd_entry': first_odd_value,
                                                    'nonzero_odd_entries': len(odd_entries)}})
    directory = Path(__file__).parent
    names = (Path(__file__).name, 'reflection_laurent_physical_frames.py',
             'check_reflection_laurent_frames.py', 'reflection_dual_rail_gaussian.py',
             'physical_split_crossings.py', 'f16_weighted_local_oracle.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    result = {'scope': 'Exact finite coupled rational-frame construction, bounded-pole catalogues, and independent weighted-profile controls; no pressure or factoring improvement certified.',
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names},
              'cases': cases, 'profile_controls': profile_records,
              'totals': {'complete_exact_profiles': sum(row['complete_profiles'] for row in profile_records),
                         'independent_coefficient_entries': sum(sum(row['independent_coefficient_entries'] for row in case['selected_full_coefficient_controls']) for case in cases),
                         'involution_entries': sum(sum(row['involution_entries'] for row in case['selected_full_coefficient_controls']) for case in cases),
                         'endpoint_entries': sum(sum(row['endpoint_entries'] for row in case['selected_full_coefficient_controls']) for case in cases)}}
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(json.dumps(result['totals']))


if __name__ == '__main__':
    main()
