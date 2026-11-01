#!/usr/bin/env python3
"""Complete F4 two-coordinate frames that can change a protected low letter.

For each interior index 0,...,4, support A and B on that index and index7.
All 4^4 coefficient choices are tested with reciprocal height flow one.
The 1280 complete profiles and numerical fits do not prove a new bound.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time

from check_coherent_hot_boundary import logical, natural_profiles
from moving_physical_frames import compress_frame
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from search_moving_interior_crossings import pressure

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a fresh path to preserve prior finite results.')
    letters = [0, 1, 2, 3, 4, 5, 8, 6, 9]
    raw = restricted_tensor(16, list(map(logical, letters)))
    energies = [0, 4, 12, 16, 20, 24, 28, 32]
    shifts = [(32 - value) // 4 for value in energies]
    names = ['tools/search_f4_interior_frames.py', 'tools/moving_physical_frames.py',
             'tools/moving_boundary_crossings.py', 'tools/search_moving_interior_crossings.py',
             'tools/check_coherent_hot_boundary.py', 'tools/restrict_implicit_reflection.py',
             'tools/native_matrix_smith_f4.py', 'tools/native_matrix_smith_f4.cpp',
             'tools/native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'started_utc': datetime.now(timezone.utc).isoformat(),
              'planned_cases': 1280, 'completed': False, 'records': [], 'energies': energies,
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in names}}
    began = time.monotonic()
    for interior in range(5):
        for ac in range(16):
            av = [0] * 8
            av[interior], av[7] = ac % 4, ac // 4
            a = [(value & 1, value >> 1) for value in av]
            for bc in range(16):
                bv = [0] * 8
                bv[interior], bv[7] = bc % 4, bc // 4
                b = [(value & 1, value >> 1) for value in bv]
                seed = compress_frame(raw, a, b, shifts, [-value for value in shifts])
                sides = natural_profiles(seed, energies)
                fit = pressure(sides, energies)
                record = {'interior_index': interior, 'interior_letter': letters[interior],
                          'a_code': ac, 'b_code': bc, 'a_coefficients': av, 'b_coefficients': bv,
                          'orders_infinity_zero': sides, 'degree': seed['degree'],
                          'tensor_sha256': digest(seed['values']), 'numerical_pressure': fit}
                report['records'].append(record)
                if 'best' not in report or fit['rho'] < report['best']['numerical_pressure']['rho'] - 1e-12:
                    report['best'] = record
                    print(json.dumps({'completed_cases': len(report['records']), 'interior_index': interior,
                                      'a_code': ac, 'b_code': bc, 'numerical_rho': fit['rho']}), flush=True)
        report['elapsed_seconds'] = time.monotonic() - began
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
        temporary.replace(args.output)
        print(json.dumps({'completed_cases': len(report['records']), 'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    assert len(report['records']) == report['planned_cases']
    report.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat())
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
