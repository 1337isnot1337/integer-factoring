#!/usr/bin/env python3
"""Replay the completed65536-word q4 Laurent-series search checkpoint."""
from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure
from check_dense_coupled_reflection_frames import apply_matrix, coefficients, matrix_product
from laurent_series_coupled_frames import construct
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize, digest
from search_affine_coupled_reflection_frames import balanced
from search_reflection4_laurent_series import decode


def independent_coefficients(base, factors, matrices):
    values = [coefficients(value) for value in base['values']]
    values = apply_matrix(values, factors['inverse_factor'], 0, 4)
    values = apply_matrix(values, list(map(list, zip(*factors['first_factor']))), 3, 4)
    for first in matrices.values():
        for second in matrices.values():
            assert matrix_product(first, second) == [[0] * 4 for _ in range(4)]
    # All products G_i G_j vanish.  Hence the sequential factors on a leg
    # equal I+sum_j w^-j G_j, including the transpose on the paired legs.
    for leg in range(4):
        for pole, matrix in sorted(matrices.items()):
            local = matrix if leg % 2 == 0 else list(map(list, zip(*matrix)))
            values = apply_matrix(values, local, leg, 4, -pole, identity=True)
    values = apply_matrix(values, factors['first_factor'], 0, 4)
    values = apply_matrix(values, list(map(list, zip(*factors['inverse_factor']))), 3, 4)
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path('evidence/reflection4-laurent-series-search.json.gz'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.archive.read_bytes()
    archive = json.loads(gzip.decompress(raw))
    for name, expected in archive['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected
    screens = archive['screen_records']
    assert archive['completed_full_family']
    assert archive['completed_screens'] == len(screens) == 65536
    assert [row[0] for row in screens] == list(range(65536))
    assert all(row[1] is None or math.isfinite(row[1]) for row in screens)
    assert sum(row[1] is None for row in screens) <= archive['unresolved_calls']
    refinements = archive['refinements']
    assert archive['completed_refinements'] == len(refinements) == 32
    assert all(row['optimizer_completed'] for row in refinements)
    assert [row['index'] for row in refinements] == [row['index'] for row in archive['retained']]
    base = tensor(4)
    bell_data = bell_normalize(base)
    seed_cache, coefficient_entries = {}, 0
    selected = sorted({row['index'] for row in refinements} | {1, 255, 256, 65535})
    for index in selected:
        matrices = decode(index)
        seed = base if index == 0 else construct(base, matrices, bell_data=bell_data)
        reference = independent_coefficients(base, bell_data[1], matrices)
        assert all(all(exponent >= 0 for exponent in entry) for entry in reference)
        assert reference == [coefficients(value) for value in seed['values']]
        coefficient_entries += len(reference)
        seed_cache[index] = seed
    replays, maximum_error = [], 0.0
    for row in refinements:
        seed = seed_cache[row['index']]
        assert digest(seed['values']) == row['tensor_sha256']
        value, _, details = balanced(AsymmetricPressure(seed, 1), np.asarray(row['best_parameters']))
        error = abs(value - row['best_rho'])
        assert error < 2e-8
        maximum_error = max(maximum_error, error)
        replays.append({'index': row['index'], 'recorded_rho': row['best_rho'],
                        'replayed_rho': float(value), 'tensor_sha256': row['tensor_sha256']})
    best = min(refinements, key=lambda row: row['best_rho'])
    nonzero = min((row for row in refinements if row['index']), key=lambda row: row['best_rho'])
    names = (Path(__file__).name, 'search_reflection4_laurent_series.py',
             'check_dense_coupled_reflection_frames.py', 'laurent_series_coupled_frames.py',
             'asymmetric_crossing_pressure.py', 'search_affine_coupled_reflection_frames.py')
    report = {'status': 'passed',
              'scope': 'Complete finite family accounting, independent full-coefficient replay, and fresh floating winner costs; no global pressure optimum or factoring bound.',
              'archive_sha256': sha256(raw).hexdigest(), 'completed_screens': len(screens),
              'completed_refinements': len(refinements), 'requested_seconds': archive['arguments']['seconds'],
              'actual_elapsed_seconds': archive['elapsed_seconds'],
              'unresolved_calls': archive['unresolved_calls'], 'independent_tensor_replays': len(selected),
              'independent_coefficient_entries': coefficient_entries, 'fresh_pressure_replays': len(replays),
              'maximum_replay_error': maximum_error, 'baseline_rho': archive['baseline_rho'],
              'best_word': best['index'], 'best_rho': best['best_rho'],
              'best_nonzero_word': nonzero['index'], 'best_nonzero_rho': nonzero['best_rho'],
              'replays': replays,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key not in ('replays', 'source_sha256')}))


if __name__ == '__main__':
    main()
