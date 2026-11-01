#!/usr/bin/env python3
"""Complete-profile search of homogeneous nonlinear physical quotient maps."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from random import Random
import sys
import time

from coherent_graded_quotients import monomial_parameters, parameter_images, quotient_maps
from coherent_hot_mixed import profiles
from search_coherent_core_products import pressure_grid
from search_coherent_hot_retractions import _atomic_write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--virtual-modes', type=int, default=4)
    parser.add_argument('--physical-modes', type=int, default=8)
    parser.add_argument('--max-trials', type=int, default=100000)
    parser.add_argument('--seconds', type=float, default=7200)
    parser.add_argument('--seed', type=int, default=104729)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    parameters = monomial_parameters(args.virtual_modes, args.physical_modes)
    if len(parameters) > 20:
        raise ValueError('Use at most twenty explicit binary parameters for this finite screen.')
    total = 1 << len(parameters)
    priority = [0, total - 1] + [1 << j for j in range(len(parameters))]
    priority = list(dict.fromkeys(priority))
    remaining = [mask for mask in range(total) if mask not in set(priority)]
    Random(args.seed).shuffle(remaining)
    queue = priority + remaining
    names = (Path(__file__).name, 'coherent_graded_quotients.py',
             'coherent_hot_mixed.py', 'coherent_hot_mixed.cpp', 'coherent_hot_compounds.py',
             'coherent_hot_compounds_fast.cpp', 'search_coherent_core_products.py',
             'search_coherent_hot_retractions.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'completed': False, 'records': [],
              'virtual_modes': args.virtual_modes, 'physical_modes': args.physical_modes,
              'parameters': parameters, 'total_parameter_masks': total,
              'trial_budget': args.max_trials, 'seconds_budget': args.seconds, 'seed': args.seed,
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started, exhausted, best = time.monotonic(), True, None
    for mask in queue:
        if len(report['records']) >= args.max_trials or time.monotonic() - started >= args.seconds:
            exhausted = False
            break
        images = parameter_images(args.virtual_modes, args.physical_modes, mask)
        code = quotient_maps(images, args.virtual_modes)
        report['in_progress'] = {'mask': mask, 'code': code}
        _atomic_write(args.output, report)
        profile = profiles(code['embedding'], code['retraction'], args.physical_modes)
        fit = pressure_grid(profile)
        record = {'mask': mask, 'code': code, 'profile': profile, 'pressure': fit}
        report['records'].append(record)
        report['in_progress'] = None
        report['elapsed_seconds'] = time.monotonic() - started
        if best is None or fit['best']['rho'] < best['pressure']['best']['rho']:
            best = record
            report['best_mask'] = mask
            print(json.dumps({'mask': mask, 'q': code['q'], 'best': fit['best'],
                              'elapsed_seconds': report['elapsed_seconds']}), flush=True)
        _atomic_write(args.output, report)
    report.update(completed=True, finite_list_exhausted=exhausted,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    _atomic_write(args.output, report)


if __name__ == '__main__':
    main()
