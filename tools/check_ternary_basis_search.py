#!/usr/bin/env python3
"""Validate and replay the best retained candidates in a bounded search archive."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from search_ternary_physical_bases import change_basis, inverse, oracles_for
import ternary_mixed_pressure as pressure


def check(archive):
    totals = {'stages': 0, 'bases': 0, 'completed_optimizations': 0,
              'exact_weighted_profile_controls': 0, 'finite_difference_gradient_controls': 0}
    stages = []
    for stage in archive['stages']:
        data = stage['data']
        assert len(data['records']) > 0 and data['elapsed_seconds'] > 0
        runs = profile_controls = gradient_controls = 0
        for index, record in enumerate(data['records']):
            assert record['index'] == index
            left, right = np.array(record['basis_a']), np.array(record['basis_b'])
            assert np.array_equal(left @ inverse(left) % 3, np.eye(3, dtype=int))
            assert np.array_equal(right @ inverse(right) % 3, np.eye(3, dtype=int))
            assert [item['run'] for item in record['runs']] == list(range(len(record['runs'])))
            assert record['best_rho'] <= data['baseline_rho'] + 1e-10
            assert record['best_rho'] <= min((item['rho'] for item in record['runs']), default=2) + 1e-10
            assert record['best_rho'] >= data['best_rho'] - 1e-10
            runs += len(record['runs'])
            profile_controls += record['exact_weighted_profile_controls']
            gradient_controls += record['finite_difference_gradient_controls']
        left, right = np.array(data['best_basis_a']), np.array(data['best_basis_b'])
        pressure.ORACLES, minor_count = oracles_for(change_basis(left, right))
        replay = pressure.objective(np.array(data['best_vector']))
        assert abs(replay - data['best_rho']) < 1e-11
        totals['stages'] += 1
        totals['bases'] += len(data['records'])
        totals['completed_optimizations'] += runs
        totals['exact_weighted_profile_controls'] += profile_controls
        totals['finite_difference_gradient_controls'] += gradient_controls
        stages.append({'name': stage['name'], 'recorded_elapsed_seconds': data['elapsed_seconds'],
                       'complete': data.get('complete', False), 'bases': len(data['records']),
                       'completed_optimizations': runs, 'recorded_best': data['best_rho'],
                       'replayed_best': replay, 'replayed_nonzero_minor_valuations': minor_count})
    return {'status': 'pass', 'counts': totals, 'stages': stages,
            'scope': 'candidate replay and complete recorded-trial accounting; no global optimum claim'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('evidence/ternary-physical-basis-search.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    result = check(json.loads(raw))
    result['archive_sha256'] = sha256(raw).hexdigest()
    result['checker_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
