#!/usr/bin/env python3
"""Independent coefficient, exact-profile and gradient controls for closed places."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random

import numpy as np
from scipy.special import logsumexp

from affine_place_profile_oracle import AffinePlaceOracle
from f4_closed_place_oracle import F4ClosedPlaceOracle, checked_factor
from physical_split_crossings import F4_PRODUCT
from search_finite_fiber_restrictions import restriction


def extension_multiply(left, right, factor):
    d = len(factor) - 1
    value = [0] * (2 * d - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            value[i + j] ^= F4_PRODUCT[a][b]
    for k in range(2 * d - 2, d - 1, -1):
        scalar = value[k]
        for j in range(d):
            value[k - d + j] ^= F4_PRODUCT[scalar][factor[j]]
    return value[:d]


def extension_evaluate(polynomial, value, factor):
    d = len(factor) - 1
    answer = [0] * d
    for k in range(max(polynomial, default=0), -1, -1):
        answer = extension_multiply(answer, value, factor)
        answer[0] ^= polynomial.get(k, 0)
    return answer


def f4_evaluate(pair, point):
    answer = 0
    for k in range(max(pair[0].bit_length(), pair[1].bit_length()) - 1, -1, -1):
        answer = F4_PRODUCT[answer][point]
        answer ^= ((pair[0] >> k) & 1) ^ (((pair[1] >> k) & 1) << 1)
    return answer


def coefficient_controls(oracle, factor, cut):
    expanded, d = oracle.local_matrix(cut, factor)
    assert d == len(factor) - 1
    alpha = [0] * d
    if d == 1:
        alpha[0] = factor[0]
    else:
        alpha[1] = 1
    controls = 0
    for i, row in enumerate(oracle.cuts[cut]):
        for j, polynomial in enumerate(row):
            # Independent Horner multiplication by alpha+t, without the
            # module's Lucas-binomial formula or cached alpha powers.
            series = [[0] * d for _ in range(oracle.degree + 1)]
            for k in range(oracle.degree, -1, -1):
                updated = [extension_multiply(value, alpha, factor) for value in series]
                for order in range(1, oracle.degree + 1):
                    updated[order] = [a ^ b for a, b in zip(updated[order], series[order - 1])]
                updated[0][0] ^= polynomial.get(k, 0)
                series = updated
            for order, expected in enumerate(series):
                for column in range(d):
                    basis = [int(k == column) for k in range(d)]
                    target = extension_multiply(expected, basis, factor)
                    actual = []
                    for output in range(d):
                        low, high = expanded[i * d + output][j * d + column]
                        actual.append(((low >> order) & 1) ^ (((high >> order) & 1) << 1))
                    assert actual == target
                    controls += d
    return controls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=Path('evidence/deflated-finite-fiber-input.json'))
    parser.add_argument('--catalogue', type=Path,
                        default=Path('evidence/finite-fiber-restriction-catalogue.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    input_bytes, catalogue_bytes = args.input.read_bytes(), args.catalogue.read_bytes()
    inputs, catalogue = json.loads(input_bytes), json.loads(catalogue_bytes)
    selected = []
    for degree in range(1, 8):
        for row in catalogue['records']:
            candidates = [p for p in row['closed_points']
                          if p['residue_degree'] == degree and p['factor'] != [0, 1]]
            if candidates:
                selected.append((row, tuple(candidates[0]['factor'])))
                break
        else:
            raise AssertionError('The catalogue no longer covers every tested residue degree.')
    for dimensions in ((2, 3), (3, 2), (3, 3)):
        row = next(r for r in catalogue['records']
                   if (len(r['subset_x']), len(r['subset_y'])) == dimensions
                   and any(p['residue_degree'] > 1 for p in r['closed_points']))
        point = next(tuple(p['factor']) for p in row['closed_points'] if p['residue_degree'] > 1)
        selected.append((row, point))
    rng = random.Random(313017)
    records = []
    coefficient_entries = profile_channels = gradient_coordinates = comparison_profiles = 0
    maximum_order_error = maximum_gradient_error = 0.0
    for row, point in selected:
        seed, _ = restriction(inputs['cases'][row['source_index']]['seed'],
                              row['subset_x'], row['subset_y'])
        oracle = F4ClosedPlaceOracle(seed)
        heights = [Fraction(rng.randrange(-100, 101), 137 + 2 * j)
                   for j in range(sum(oracle.alphabets))]
        scale = Fraction(13, 11)
        current = []
        for cut in (0, 1):
            coefficient_entries += coefficient_controls(oracle, point, cut)
            for place in (point, None):
                exact = oracle.profile(cut, place, heights, scale)
                value, gradient, dscale, native = oracle.log_profile_with_gradient(
                    cut, place, list(map(float, heights)), float(scale))
                error = max(abs(float(a) - b) for a, b in zip(exact, native))
                assert len(exact) == oracle.n and error < 3e-8
                assert abs(value - logsumexp([-float(x) for x in exact])) < 3e-8
                maximum_order_error = max(maximum_order_error, error)
                profile_channels += oracle.n
                for leg in range(4):
                    start = oracle.offsets[leg]
                    assert abs(sum(gradient[start:start + oracle.alphabets[leg]]) + 1) < 2e-9
                # The nonzero rational heights make these local evaluations
                # generic enough for a centered finite-difference control.
                variables = np.r_[list(map(float, heights)), float(scale)]
                analytic = np.r_[gradient, dscale]
                for coordinate in range(len(variables)):
                    step = 2e-6
                    plus, minus = variables.copy(), variables.copy()
                    plus[coordinate] += step
                    minus[coordinate] -= step
                    high = oracle.log_profile_with_gradient(cut, place, plus[:-1], plus[-1])[0]
                    low = oracle.log_profile_with_gradient(cut, place, minus[:-1], minus[-1])[0]
                    observed = (high - low) / (2 * step)
                    error = abs(observed - analytic[coordinate])
                    assert error < 4e-6, (row['subset_x'], row['subset_y'], point,
                                         cut, place, coordinate, observed, analytic[coordinate])
                    maximum_gradient_error = max(maximum_gradient_error, error)
                    gradient_coordinates += 1
                current.append({'cut': cut, 'place': place, 'orders': [str(v) for v in exact],
                                'native_error': max(abs(float(a) - b) for a, b in zip(exact, native))})
        records.append({'source_index': row['source_index'], 'subset_x': row['subset_x'],
                        'subset_y': row['subset_y'], 'degree': seed['degree'],
                        'factor': point, 'residue_degree': len(point) - 1,
                        'heights': [str(v) for v in heights], 'scale': str(scale), 'profiles': current})
    # Compare the old rational-place adapter without altering its source.
    full_seed = inputs['cases'][0]['seed']
    old, new = AffinePlaceOracle(full_seed), F4ClosedPlaceOracle(full_seed)
    heights = np.asarray([rng.uniform(-1, 1) for _ in range(16)])
    for cut in (0, 1):
        for place in ((0, 1), (2, 1), (3, 1), None):
            first = old.log_profile_with_gradient(cut, place, heights, .87)
            second = new.log_profile_with_gradient(cut, place, heights, .87)
            assert abs(first[0] - second[0]) < 2e-9
            assert np.max(np.abs(first[1] - second[1])) < 2e-8
            assert abs(first[2] - second[2]) < 2e-8
            assert np.max(np.abs(first[3] - second[3])) < 2e-8
            comparison_profiles += 1
    rejected = 0
    for factor in ((1,), (0, 2), (1, 0, 1), (0, 1, 1), (0, 0, 1), (4, 1)):
        try:
            checked_factor(factor)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('An invalid or reducible place was accepted.')
    names = (Path(__file__).name, 'f4_closed_place_oracle.py',
             'affine_place_profile_oracle.py', 'fused_physical_basis_oracle.py',
             'check_fused_closed_profiles.py', 'f16_weighted_local_oracle.py',
             'physical_split_crossings.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'search_finite_fiber_restrictions.py')
    result = {
        'scope': 'Exact finite-field restriction, complete rational-height profiles, and charged local gradients; no pressure threshold claim.',
        'input_sha256': sha256(input_bytes).hexdigest(),
        'catalogue_sha256': sha256(catalogue_bytes).hexdigest(),
        'records': records,
        'counts': {'independent_extension_coefficient_entries': coefficient_entries,
                   'complete_exact_profiles': sum(len(r['profiles']) for r in records),
                   'complete_exact_channels': profile_channels,
                   'finite_difference_gradient_coordinates': gradient_coordinates,
                   'old_rational_adapter_comparisons': comparison_profiles,
                   'rejected_invalid_places': rejected},
        'maximum_order_error': maximum_order_error,
        'maximum_gradient_error': maximum_gradient_error,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names}, 'passed': True,
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('counts', 'maximum_order_error',
                                                 'maximum_gradient_error', 'passed')}))


if __name__ == '__main__':
    main()
