#!/usr/bin/env python3
"""Compute one explicit rank-three height-25 coherent q96 candidate.

The map input and source hashes are saved before the exact complete-profile
calculation. Numerical pressure is a discovery value, not a certificate.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time

from coherent_hot_compounds_fast import profiles
from search_coherent_hot_retractions import _atomic_write, _lift, pressure, q96_start


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--right', default='3,5,9')
    parser.add_argument('--left', default='13,11,8')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    right, left = [list(map(int, value.split(','))) for value in (args.right, args.left)]
    if len(right) != 3 or len(left) != 3:
        raise ValueError('Use exactly three paired boundary vectors.')
    base, boundary, _, _ = q96_start()
    embedding = base + [_lift(value, boundary) for value in right]
    retraction = base + [_lift(value, boundary) for value in left]
    names = (Path(__file__).name, 'coherent_hot_compounds_fast.py', 'coherent_hot_compounds_fast.cpp',
             'coherent_hot_compounds.py', 'search_coherent_hot_retractions.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'completed': False,
              'description': {'right': right, 'left': left, 'boundary_words': boundary,
                              'fixed_boundary_word': 31, 'modes': 12},
              'embedding': embedding, 'retraction': retraction,
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    _atomic_write(args.output, report)
    start = time.monotonic()
    result = profiles(embedding, retraction, 12)
    fit = pressure(result)
    report.update(completed=True, profile=result, pressure=fit,
                  elapsed_seconds=time.monotonic() - start,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    _atomic_write(args.output, report)
    print(json.dumps({'completed': True, 'pressure': fit, 'seconds': report['elapsed_seconds'],
                      'block_precisions': [[b['precision_t'] for b in c['blocks']]
                                            for c in result['cuts']]}, indent=2), flush=True)


if __name__ == '__main__':
    main()
