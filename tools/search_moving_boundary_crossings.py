#!/usr/bin/env python3
"""Complete declared finite screens of polynomial reflection boundary maps.

Both Smith lists contain every promised pivot. Fits are numerical search
evidence, not exact pressure certificates or factoring improvements.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time

from check_coherent_hot_boundary import logical, natural_profiles
from hot_code_charge_profiles import hot_energy
from moving_boundary_crossings import compress
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from search_coherent_hot_boundaries import pressure

ROOT = Path(__file__).resolve().parents[1]


def polynomial(code, field, degree):
    planes = [0, 0]
    for power in range(degree + 1):
        coefficient = code % field
        code //= field
        planes[0] |= (coefficient & 1) << power
        planes[1] |= (coefficient >> 1) << power
    return tuple(planes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field', type=int, choices=(2, 4), default=4)
    parser.add_argument('--degree', type=int, choices=(0, 1, 2, 3), default=1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    lower, boundary = [0, 1, 2, 3, 4, 5, 8], [6, 9]
    raw = restricted_tensor(16, list(map(logical, lower + boundary)))
    energies = list(map(hot_energy, lower)) + [32]
    count = args.field ** (args.degree + 1)
    source_names = ['search_moving_boundary_crossings.py', 'moving_boundary_crossings.py',
                    'check_coherent_hot_boundary.py', 'restrict_implicit_reflection.py',
                    'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
                    'native_rectangular_smith_f4.cpp', 'search_coherent_hot_boundaries.py']
    sources = [str(Path(__file__).with_name(name)) for name in source_names]
    report = {'scope': __doc__, 'started_utc': datetime.now(timezone.utc).isoformat(),
              'field': args.field, 'parameter_degree_bound': args.degree,
              'planned_pairs': count * count, 'records': [], 'completed': False,
              'source_sha256': {str(Path(p).relative_to(ROOT)): sha256(Path(p).read_bytes()).hexdigest()
                                for p in sources},
              'lower_letters': lower, 'boundary_letters': boundary, 'energies': energies}
    began = time.monotonic()
    best = None
    for ai in range(count):
        for bi in range(count):
            a, b = polynomial(ai, args.field, args.degree), polynomial(bi, args.field, args.degree)
            seed = compress(raw, len(lower), a, b)
            sides = natural_profiles(seed, energies)
            fit = pressure(sides, energies)
            record = {'a': a, 'b': b, 'degree': seed['degree'], 'orders_infinity_zero': sides,
                      'tensor_sha256': digest(seed['values']), 'numerical_pressure': fit}
            report['records'].append(record)
            if best is None or fit['rho'] < best['numerical_pressure']['rho'] - 1e-12:
                best = record
                report['best'] = best
                print(json.dumps({'completed_pairs': len(report['records']), 'a': a, 'b': b,
                                  'degree': seed['degree'], 'numerical_rho': fit['rho']}), flush=True)
        report['elapsed_seconds'] = time.monotonic() - began
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
        temporary.replace(args.output)
    report.update(completed=True, completed_pairs=len(report['records']),
                  elapsed_seconds=time.monotonic() - began,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: report[k] for k in ('completed', 'completed_pairs', 'elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    main()
