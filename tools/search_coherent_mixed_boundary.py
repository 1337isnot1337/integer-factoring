#!/usr/bin/env python3
"""All 496 binary q32 height-boundary retractions, allowing mixed occupation.

The 121 occupation-preserving maps reuse the already complete exact
catalogue. Every additional map is evaluated with both full parity-block
profiles. Numerical pressure is exploratory; all inputs and profiles persist.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time

from coherent_hot_compounds import height
from coherent_hot_mixed import profiles, validate_maps
from search_coherent_hot_retractions import _atomic_write, _lift, pressure


ROOT = Path(__file__).resolve().parents[1]
BOUNDARY = (24, 36, 66, 129, 15)


def boundary_maps(kernel, covector):
    if not 1 <= kernel < 32 or not 1 <= covector < 32 or (kernel & covector).bit_count() % 2 != 1:
        raise ValueError('Use nonzero five-dimensional masks with odd dot product.')
    lower = sorted((word for word in range(256) if height(word) < 16),
                   key=lambda word: (height(word), word))
    assert len(lower) == 28
    base = [1 << word for word in lower]
    pivot = (covector & -covector).bit_length() - 1
    keep = [i for i in range(5) if i != pivot]
    right = [(1 << i) ^ (((covector >> i) & 1) << pivot) for i in keep]
    left = [(1 << i) ^ (covector if kernel >> i & 1 else 0) for i in keep]
    embedding = base + [_lift(value, BOUNDARY) for value in right]
    retraction = base + [_lift(value, BOUNDARY) for value in left]
    return {'kernel': kernel, 'covector': covector, 'right': right, 'left': left,
            'boundary_words': list(BOUNDARY)}, embedding, retraction


def catalogue_keys():
    priority = [(23, 15), (7, 31), (31, 31), (31, 7), (15, 23)]
    keys = [(kernel, covector) for kernel in range(1, 32) for covector in range(1, 32)
            if (kernel & covector).bit_count() % 2 == 1]
    assert len(keys) == 496 and all(key in keys for key in priority)
    return priority + [key for key in keys if key not in priority]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=7200)
    parser.add_argument('--max-trials', type=int, default=496)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.seconds <= 0 or args.max_trials <= 0:
        raise ValueError('Use positive time and trial budgets.')
    old_path = ROOT / 'evidence/coherent-hot-q32-retractions.json'
    old = json.loads(old_path.read_text())
    assert old['completed'] and old['finite_list_exhausted'] and old['unresolved'] == 0
    reused = {}
    for record in old['records']:
        description = record['description']
        key = ((description['kernel'], description['covector']) if 'kernel' in description
               else (16, 16))
        reused[key] = record
    assert len(reused) == 121
    names = (Path(__file__).name, 'coherent_hot_mixed.py', 'coherent_hot_mixed.cpp',
             'coherent_hot_compounds.py', 'coherent_hot_compounds_fast.cpp',
             'search_coherent_hot_retractions.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'completed': False, 'total_maps': 496, 'records': [],
              'seconds_budget': args.seconds, 'trial_budget': args.max_trials,
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'previous_catalogue_sha256': sha256(old_path.read_bytes()).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'unresolved': 0, 'in_progress': None}
    started, best, exhausted = time.monotonic(), None, True
    for kernel, covector in catalogue_keys():
        if len(report['records']) >= args.max_trials or time.monotonic() - started >= args.seconds:
            exhausted = False
            break
        description, embedding, retraction = boundary_maps(kernel, covector)
        _, _, heights, checks = validate_maps(embedding, retraction, 8)
        report['in_progress'] = {'description': description, 'embedding': embedding,
                                 'retraction': retraction}
        _atomic_write(args.output, report)
        trial_start = time.monotonic()
        record = {'description': description, 'embedding': embedding, 'retraction': retraction,
                  'heights': list(heights), 'deletion_checks': checks, 'completed': False}
        try:
            if (kernel, covector) in reused:
                previous = reused[kernel, covector]
                assert sorted(heights) == sorted(previous['heights'])
                record.update(zero_orders_t=previous['zero_orders_t'],
                              infinity_orders_t=previous['infinity_orders_t'],
                              pressure=previous['pressure'], reused_occupation_profile=True)
            else:
                result = profiles(embedding, retraction, 8)
                record.update(zero_orders_t=result['zero_orders_t'],
                              infinity_orders_t=result['infinity_orders_t'],
                              pressure=pressure(result), reused_occupation_profile=False,
                              block_sizes=[b['size'] for b in result['cuts'][0]['blocks']],
                              precisions=[[b['precision_t'] for b in c['blocks']] for c in result['cuts']])
            record['completed'] = True
            if best is None or record['pressure']['rho'] < best['pressure']['rho']:
                best = record
                report['best_key'] = [kernel, covector]
                print(json.dumps({'best_key': [kernel, covector], 'rho': record['pressure']['rho'],
                                  'elapsed_seconds': time.monotonic() - started}), flush=True)
        except (ValueError, ArithmeticError) as error:
            record['error'] = f'{type(error).__name__}: {error}'
            report['unresolved'] += 1
        record['seconds'] = time.monotonic() - trial_start
        report['records'].append(record)
        report['in_progress'] = None
        report['elapsed_seconds'] = time.monotonic() - started
        _atomic_write(args.output, report)
    report.update(completed=True, finite_list_exhausted=exhausted,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    _atomic_write(args.output, report)


if __name__ == '__main__':
    main()
