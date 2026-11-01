#!/usr/bin/env python3
"""Reproducible binary polynomial codimension-one frame screen.

For every declared A,B in {0,...,255}, use a=e_7+w*A, b=B, or a
declared reciprocal monomial flow of A,B, in L=[I,a], P=[I+a*b^T;b^T].
Both complete literal weighted profiles are resolved. A numerical
pressure threshold is discovery evidence only.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from random import Random
import time

from check_coherent_hot_boundary import logical, natural_profiles
from moving_physical_frames import compress_frame
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from search_moving_interior_crossings import pressure

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--a-masks', type=int, nargs='+')
    parser.add_argument('--b-masks', type=int, nargs='+')
    parser.add_argument('--uniform-flow', type=int)
    parser.add_argument('--height-flow', type=int)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--records', type=Path, required=True)
    args = parser.parse_args()
    if args.uniform_flow is not None and args.height_flow is not None:
        parser.error('Choose at most one Laurent flow.')
    rng = Random(1952001)
    a_order, b_order = list(range(1, 256)), [i for i in range(256) if i != 128]
    rng.shuffle(a_order)
    rng.shuffle(b_order)
    a_order, b_order = args.a_masks or [0] + a_order, args.b_masks or [128] + b_order
    if any(not 0 <= value < 256 for value in a_order + b_order):
        parser.error('Masks must be integers from zero through255.')
    if len(set(a_order)) != len(a_order) or len(set(b_order)) != len(b_order):
        parser.error('A declared mask list must not contain duplicates.')
    if args.records.exists() or args.output.exists():
        parser.error('Use fresh output paths to preserve prior receipts.')
    raw = restricted_tensor(16, list(map(logical, [0, 1, 2, 3, 4, 5, 8, 6, 9])))
    energies = [0, 4, 12, 16, 20, 24, 28, 32]
    source_names = ['tools/search_moving_physical_frames.py', 'tools/moving_physical_frames.py',
                    'tools/moving_boundary_crossings.py', 'tools/search_moving_interior_crossings.py',
                    'tools/check_coherent_hot_boundary.py', 'tools/restrict_implicit_reflection.py',
                    'tools/native_matrix_smith_f4.py', 'tools/native_matrix_smith_f4.cpp',
                    'tools/native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'a_order': a_order, 'b_order': b_order,
              'planned_cases': len(a_order) * len(b_order), 'completed_cases': 0,
              'completed': False, 'energies': energies,
              'uniform_flow': args.uniform_flow, 'height_flow': args.height_flow,
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in source_names}}
    began = time.monotonic()

    def save():
        report['elapsed_seconds'] = time.monotonic() - began
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
        temporary.replace(args.output)

    with args.records.open('x') as records:
        for am in a_order:
            flow = args.uniform_flow is not None or args.height_flow is not None
            a = [(int(am >> i & 1), 0) if flow else
                 (int(i == 7) | (2 if am >> i & 1 else 0), 0) for i in range(8)]
            shifts = ([args.uniform_flow] * 8 if args.uniform_flow is not None else
                      [args.height_flow * (32 - energy) // 4 for energy in energies]
                      if args.height_flow is not None else [0] * 8)
            for bm in b_order:
                b = [(int(bm >> i & 1), 0) for i in range(8)]
                seed = compress_frame(raw, a, b, shifts, [-value for value in shifts])
                sides = natural_profiles(seed, energies)
                fit = pressure(sides, energies)
                record = {'a_mask': am, 'b_mask': bm, 'degree': seed['degree'],
                          'laurent_clearing_power': seed['laurent_clearing_power'],
                          'tensor_sha256': digest(seed['values']),
                          'orders_infinity_zero': sides, 'numerical_pressure': fit}
                records.write(json.dumps(record, sort_keys=True, separators=(',', ':')) + '\n')
                records.flush()
                report['completed_cases'] += 1
                if 'best' not in report or fit['rho'] < report['best']['numerical_pressure']['rho'] - 1e-12:
                    report['best'] = record
                    save()
                    print(json.dumps({'completed_cases': report['completed_cases'], 'a_mask': am,
                                      'b_mask': bm, 'numerical_rho': fit['rho']}), flush=True)
            save()
            print(json.dumps({'completed_cases': report['completed_cases'],
                              'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat(),
                  records_sha256=sha256(args.records.read_bytes()).hexdigest())
    save()


if __name__ == '__main__':
    main()
