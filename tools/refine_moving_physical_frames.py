#!/usr/bin/env python3
"""Bounded full physical-gauge refinements of declared moving frame pairs.

All tensor entries and endpoints are exact. Gauge weights, types and
optimized pressures are numerical discovery evidence, without a theorem.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np

from check_coherent_hot_boundary import logical, natural_profiles
from moving_physical_frames import compress_frame
from restrict_implicit_reflection import restricted_tensor
from search_moving_interior_crossings import pressure
from search_physical_code_gauges import optimize, start_params

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maxiter', type=int, default=80)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    pairs = [(128, 128), (128, 0), (0, 128), (1, 128), (128, 1),
             (170, 85), (85, 170), (240, 15), (15, 240), (204, 51), (51, 204), (255, 255)]
    raw = restricted_tensor(16, list(map(logical, [0, 1, 2, 3, 4, 5, 8, 6, 9])))
    energies = np.array([0, 4, 12, 16, 20, 24, 28, 32])
    shifts = [(32 - int(value)) // 4 for value in energies]
    sources = ['tools/refine_moving_physical_frames.py', 'tools/moving_physical_frames.py',
               'tools/moving_boundary_crossings.py', 'tools/search_moving_interior_crossings.py',
               'tools/search_physical_code_gauges.py', 'tools/native_rectangular_smith_f4.py',
               'tools/native_rectangular_smith_f4.cpp', 'tools/restrict_implicit_reflection.py']
    report = {'scope': __doc__, 'declared_pairs': pairs, 'height_flow': 1,
              'maximum_iterations_per_start': args.maxiter, 'starts_per_pair': 3,
              'completed': False, 'records': [],
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in sources}}
    began, rng = time.monotonic(), np.random.default_rng(19520014)
    for am, bm in pairs:
        a = [(int(am >> i & 1), 0) for i in range(8)]
        b = [(int(bm >> i & 1), 0) for i in range(8)]
        seed = compress_frame(raw, a, b, shifts, [-value for value in shifts])
        natural = pressure(natural_profiles(seed, energies.tolist()), energies.tolist())
        initial = start_params(energies * natural['scale'], natural['scale'])
        starts = [('natural', initial)]
        for amplitude in (.002, .02):
            changed = initial.copy()
            changed[:-1] += rng.uniform(-amplitude, amplitude, len(changed) - 1)
            starts.append((str(amplitude), changed))
        record = {'a_mask': am, 'b_mask': bm, 'natural_numerical_pressure': natural}
        try:
            best, attempts = optimize(seed, starts, args.maxiter)
            record.update(best=best, attempts=attempts)
            assert abs(attempts[0]['initial_rho'] - natural['rho']) < 1e-8
        except (ArithmeticError, ValueError) as error:
            record.update(error_type=type(error).__name__, error=str(error))
        report['records'].append(record)
        report['elapsed_seconds'] = time.monotonic() - began
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
        print(json.dumps({'a_mask': am, 'b_mask': bm,
                          'natural_rho': natural['rho'],
                          'refined_rho': record.get('best', {}).get('rho'),
                          'error': record.get('error'), 'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat())
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
