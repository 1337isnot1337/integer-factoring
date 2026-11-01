#!/usr/bin/env python3
"""Finite complete search of weighted homogeneous binary inverse-system codes.

All nonzero top forms in one stated monomial cell are listed. Complete
two-cut profiles determine each numerical pressure; unfinished budgets
retain the exact in-progress form and never imply exhaustive coverage.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time

from coherent_hot_compounds import height, profiles
from coherent_inverse_systems import inverse_system
from search_coherent_hot_retractions import _atomic_write, pressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', type=int, default=8)
    parser.add_argument('--particles', type=int, default=4)
    parser.add_argument('--height', type=int, default=32)
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--max-trials', type=int, default=255)
    parser.add_argument('--start-mask', type=int, default=1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.modes <= 16 or not 0 <= args.particles <= args.modes or args.seconds <= 0:
        raise ValueError('Use a valid occupation cell and a positive time budget.')
    monomials = [w for w in range(1 << args.modes)
                 if w.bit_count() == args.particles and height(w) == args.height]
    if not monomials or len(monomials) > 24:
        raise ValueError('Use a nonempty explicit cell of at most twenty-four monomials.')
    names = (Path(__file__).name, 'coherent_inverse_systems.py', 'coherent_hot_compounds.py',
             'coherent_hot_compounds.cpp', 'search_coherent_hot_retractions.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'modes': args.modes, 'particles': args.particles,
              'height': args.height, 'top_monomials': monomials,
              'total_nonzero_forms': (1 << len(monomials)) - 1, 'start_mask': args.start_mask,
              'seconds_budget': args.seconds, 'trial_budget': args.max_trials,
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'records': [], 'completed': False,
              'unresolved': 0, 'in_progress': None}
    began, best = time.monotonic(), None
    next_mask = args.start_mask
    for mask in range(args.start_mask, 1 << len(monomials)):
        if len(report['records']) >= args.max_trials or time.monotonic() - began >= args.seconds:
            break
        top = sum(1 << word for i, word in enumerate(monomials) if mask >> i & 1)
        code = inverse_system([top], args.modes)
        report['in_progress'] = {'mask': mask, 'code': code}
        _atomic_write(args.output, report)
        started = time.monotonic()
        record = {'mask': mask, 'code': code, 'q': len(code['embedding']), 'completed': False}
        try:
            result = profiles(code['embedding'], code['retraction'], args.modes)
            fit = pressure(result)
            record.update(completed=True, pressure=fit, zero_orders_t=result['zero_orders_t'],
                          infinity_orders_t=result['infinity_orders_t'],
                          block_sizes=[b['size'] for b in result['cuts'][0]['blocks']],
                          precisions=[[b['precision_t'] for b in cut['blocks']] for cut in result['cuts']])
            if best is None or fit['rho'] < best['pressure']['rho']:
                best = record
                report['best_mask'] = mask
                print(json.dumps({'mask': mask, 'q': record['q'], 'rho': fit['rho'],
                                  'elapsed': time.monotonic() - began}), flush=True)
        except (ValueError, ArithmeticError) as error:
            record['error'] = f'{type(error).__name__}: {error}'
            report['unresolved'] += 1
        record['seconds'] = time.monotonic() - started
        report['records'].append(record)
        report['in_progress'] = None
        report['elapsed_seconds'] = time.monotonic() - began
        next_mask = mask + 1
        _atomic_write(args.output, report)
    report.update(completed=True, next_mask=next_mask,
                  finite_list_exhausted=next_mask == 1 << len(monomials),
                  elapsed_seconds=time.monotonic() - began,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    _atomic_write(args.output, report)


if __name__ == '__main__':
    main()
