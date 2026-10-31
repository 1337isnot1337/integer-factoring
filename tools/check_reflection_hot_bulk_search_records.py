#!/usr/bin/env python3
"""Replay saved best points and schedule counts from five finite bulk searches."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from reflection_hot_leg_bulk_gauges import LegBulkHotStatePressure
from reflection_hot_all_bulk_gauges import AllBulkHotStatePressure
from search_reflection_hot_core_filters import factor_from_row
from search_reflection_hot_all_bulk_gauges import derivative_controls


FILES = (
    'reflection-hot-leg-bulk-search-M64.json',
    'reflection-hot-leg-bulk-search-M128.json',
    'reflection-hot-all-bulk-search-M64.json',
    'reflection-hot-all-bulk-search-M128.json',
    'reflection-hot-scaled-bulk-search-M128.json',
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    summary = {
        'status': 'Complete saved-point replay, not a new optimization or global minimum certificate.',
        'chain_rule_controls': derivative_controls(),
        'records': [],
    }
    for name in FILES:
        path = root / 'evidence' / name
        saved = json.loads(path.read_text())
        assert saved['complete']
        assert len(saved['records']) == saved['planned_candidates'] == saved['completed_candidates']
        interval = saved['interval']
        count, starts, error = 0, 0, 0.0
        for record in saved['records']:
            model_class = LegBulkHotStatePressure if '-leg-' in name else AllBulkHotStatePressure
            model = model_class(interval, factor_from_row(tuple(record['row']), 1),
                                profile_precision=max(512, 8 * interval))
            points = [record['best']]
            if 'cones' in record:
                assert len(record['cones']) == 16
                assert len({tuple(c['larger_leg_choices']) for c in record['cones']}) == 16
                points += [c['best'] for c in record['cones']]
                starts += sum(len(c['starts']) for c in record['cones'])
            else:
                starts += len(record['starts'])
            for point in points:
                offsets = np.zeros((2, 4, 2))
                offsets[:, :, 1] = point['relative_hot_heights']
                bulk = np.array(point['physical_leg_bulks'])
                value, _, detail = model.evaluate(point['log_s'], offsets, bulk, True)
                error = max(error, abs(value - point['rho']))
                assert abs(value - point['rho']) < 2e-11, (name, record['row'])
                for key in ('active_orders', 'spectator_orders', 'type_x_core',
                            'type_y_core', 'cold_x_energies', 'cold_y_energies'):
                    assert np.allclose(np.array(detail[key]), np.array(point[key]),
                                       rtol=0, atol=2e-10), (name, key)
                count += 1
        row = {
            'file': name, 'sha256': sha256(path.read_bytes()).hexdigest(),
            'candidates': len(saved['records']), 'starts': starts,
            'original_queries': sum(r['queries'] for r in saved['records']),
            'replayed_saved_points': count, 'maximum_pressure_error': error,
        }
        summary['records'].append(row)
        print(json.dumps(row), flush=True)
    if args.output:
        args.output.write_text(json.dumps(summary, indent=2) + '\n')


if __name__ == '__main__':
    main()
