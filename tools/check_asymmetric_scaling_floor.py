#!/usr/bin/env python3
"""Reproduce the degenerate-scale asymmetric support-value diagnostic.

For a rectangular A,B,B,A complete profile the zero-height limit is
cos(theta)+sin(theta), independent of the seed and the cut.  The recorded
q4 frames approach that surface from the physical side and the best tested
corner stays above one.  This is numerical evidence, not a pressure or
factoring certificate.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure
from dense_coupled_reflection_frames import digest, rank_one_terms, transformed_from_terms
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize
from search_dense_coupled_companions import translate_terms


def quadratic_fit(rows):
    xs = np.asarray([row['log_s'] ** 2 for row in rows[:3]], float)
    ys = np.asarray([row['support_value'] - 1 for row in rows[:3]], float)
    return float(xs @ ys / (xs @ xs))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=Path('evidence/asymmetric-scaling-floor.json'))
    parser.add_argument('--catalogue', type=Path, default=Path('evidence/dense-rank-one-coupled-frames.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    evidence = json.loads(args.evidence.read_text())
    catalogue = json.loads(args.catalogue.read_text())
    if sha256(args.catalogue.read_bytes()).hexdigest() != evidence['source_sha256']['evidence/dense-rank-one-coupled-frames.json']:
        raise SystemExit('Catalogue hash differs from the recorded diagnostic source.')
    base = tensor(4)
    _, factors = bell_normalize(base)
    records, max_zero_error, max_scaled_error = [], 0.0, 0.0
    for case in evidence['cases']:
        source = catalogue['records'][case['generator_index']]
        if source['u'] != case['u'] or source['v'] != case['v']:
            raise SystemExit('Generator vectors differ from the catalogue.')
        terms = rank_one_terms(base, case['u'], case['v'], factors)
        if case['constant']:
            terms = translate_terms(terms, case['constant'])
        changed = transformed_from_terms(base, terms, case['pole'])
        if digest(changed['values']) != case['seed_sha256']:
            raise SystemExit('Framed tensor hash differs from the recorded seed.')
        oracle = AsymmetricPressure(changed, 1)
        ratios = np.asarray(case['ratios'], float)
        zero_errors = []
        for row in evidence['zero_scale_tests']:
            if (row['generator_index'], row['constant'], row['pole']) != (
                    case['generator_index'], case['constant'], case['pole']):
                continue
            params = np.r_[ratios * row['log_s'], row['log_s'], row['angle']]
            value, _, _ = oracle.evaluate(params)
            error = abs(value - row['support_value'])
            predicted = math.cos(row['angle']) + math.sin(row['angle'])
            zero_errors.append(max(error, abs(value - predicted)))
        scaled_errors, values = [], []
        for row in case['scaled']:
            params = np.r_[ratios * row['log_s'], row['log_s'], row['angle']]
            value, _, _ = oracle.evaluate(params)
            scaled_errors.append(abs(value - row['support_value']))
            values.append(value)
            if value <= 1:
                raise SystemExit('The tested support value is not above one.')
        if any(values[i] <= values[i + 1] for i in range(len(values) - 1)):
            raise SystemExit('The tested scaled values are not decreasing with log(s).')
        delta = quadratic_fit(case['scaled'])
        if abs(delta - case['quadratic_coefficient']) > 1e-5 or delta <= 0:
            raise SystemExit('The recorded positive quadratic coefficient did not replay.')
        max_zero_error = max(max_zero_error, max(zero_errors, default=0.0))
        max_scaled_error = max(max_scaled_error, max(scaled_errors))
        records.append({'generator_index': case['generator_index'], 'constant': case['constant'],
                        'pole': case['pole'], 'seed_sha256': case['seed_sha256'],
                        'tested_log_s': [row['log_s'] for row in case['scaled']],
                        'tested_support_values': values,
                        'quadratic_coefficient': delta,
                        'max_scaled_replay_error': max(scaled_errors),
                        'max_zero_scale_identity_error': max(zero_errors, default=0.0)})
    if max_zero_error >= 2e-6 or max_scaled_error >= 5e-12:
        raise SystemExit('A replay error exceeded its declared tolerance.')
    names = ('check_asymmetric_scaling_floor.py', 'asymmetric_crossing_pressure.py',
             'dense_coupled_reflection_frames.py', 'search_dense_coupled_companions.py',
             'native_rectangular_smith_f4.py')
    result = {'status': 'pass', 'scope': 'numerical degenerate-scale diagnostic; no rate certificate',
              'records': records, 'max_zero_scale_identity_error': max_zero_error,
              'max_scaled_replay_error': max_scaled_error,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass', 'cases': len(records),
                      'max_zero_scale_identity_error': max_zero_error,
                      'max_scaled_replay_error': max_scaled_error,
                      'quadratic_coefficients': [record['quadratic_coefficient'] for record in records]}))


if __name__ == '__main__':
    main()
