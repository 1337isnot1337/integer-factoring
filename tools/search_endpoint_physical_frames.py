#!/usr/bin/env python3
"""Physical frames with only the necessary hot split identity imposed.

Use a=e7+w*A, r=w*(A+e7), b=B. Then a(1)=r(1), while
L(w)P(w)=I+(1+w)*e7*B^T need not equal I away from the endpoint.
The literal two-cut numerical screen is not a factoring-bound certificate.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from check_coherent_hot_boundary import logical, natural_profiles
from moving_physical_frames import compress_frame
from restrict_implicit_reflection import restricted_tensor
from search_moving_interior_crossings import pressure

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    masks = [0, 1, 128, 255, 170, 85, 240, 15, 204, 51]
    raw = restricted_tensor(16, list(map(logical, [0, 1, 2, 3, 4, 5, 8, 6, 9])))
    energies = [0, 4, 12, 16, 20, 24, 28, 32]
    names = ['tools/search_endpoint_physical_frames.py', 'tools/moving_physical_frames.py',
             'tools/moving_boundary_crossings.py', 'tools/search_moving_interior_crossings.py',
             'tools/check_coherent_hot_boundary.py', 'tools/restrict_implicit_reflection.py',
             'tools/native_matrix_smith_f4.py', 'tools/native_matrix_smith_f4.cpp',
             'tools/native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'masks': masks, 'completed': False, 'records': [],
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in names}}
    began = time.monotonic()
    for am in masks:
        a = [(int(i == 7) | (2 if am >> i & 1 else 0), 0) for i in range(8)]
        r = [(2 if (am ^ 128) >> i & 1 else 0, 0) for i in range(8)]
        for bm in masks:
            b = [(int(bm >> i & 1), 0) for i in range(8)]
            seed = compress_frame(raw, a, b, right_a=r)
            sides = natural_profiles(seed, energies)
            fit = pressure(sides, energies)
            record = {'a_mask': am, 'b_mask': bm, 'degree': seed['degree'],
                      'orders_infinity_zero': sides, 'numerical_pressure': fit}
            report['records'].append(record)
            if 'best' not in report or fit['rho'] < report['best']['numerical_pressure']['rho'] - 1e-12:
                report['best'] = record
                print(json.dumps({'completed_cases': len(report['records']), 'a_mask': am,
                                  'b_mask': bm, 'numerical_rho': fit['rho']}), flush=True)
        report['elapsed_seconds'] = time.monotonic() - began
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    report['completed'] = True
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
