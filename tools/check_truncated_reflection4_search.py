#!/usr/bin/env python3
"""Reconstruct and replay every saved minimum in the complete q4 cap screen."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import digest
from search_truncated_reflection4 import balanced
from truncate_reflection_crossings import truncate_even_crossing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('evidence/truncated-reflection4-pressure-search.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    data = json.loads(args.source.read_text())
    assert data['completed'] and data['queue_exhausted'] and len(data['cases']) == 24
    assert {(r['degree_bound'], r['cut_at_infinity']) for r in data['cases']} == {
        (degree, cut) for degree in range(2, 25, 2) for cut in (0, 1)}
    for name, expected in data['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected
    base = tensor(4)
    errors = []
    for row in data['cases']:
        seed = truncate_even_crossing(base, row['degree_bound'])
        assert seed['degree'] == row['actual_degree']
        assert digest(seed['values']) == row['tensor_sha256']
        oracle = AsymmetricPressure(seed, row['cut_at_infinity'])
        for mode in ('balanced', 'asymmetric'):
            saved = row[mode]
            params = np.asarray(saved['params'])
            actual = (balanced(oracle, params) if mode == 'balanced' else oracle.evaluate(params))[0]
            errors.append(abs(actual - saved['value']))
            assert errors[-1] < 2e-10
    result = {'status': 'passed', 'scope': 'complete exact tensor reconstruction and fresh numerical minimum replays; no exact pressure certificate',
              'source_record_sha256': sha256(args.source.read_bytes()).hexdigest(),
              'reconstructed_tensors': 24, 'endpoint_entries': 24 * 256,
              'fresh_complete_profile_minima': len(errors), 'maximum_replay_error': max(errors),
              'best_rho': data['best_rho'], 'best_support_value': data['best_support_value']}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
