#!/usr/bin/env python3
"""Declared monomial joint-retraction screen with complete two-cut profiles.

The finite scale scan is numerical. No all-scale exclusion follows from it.
Frobenius identifies scalar omega and omega+1 because the source is binary.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from math import sqrt
from pathlib import Path
import time

import numpy as np
from scipy.special import logsumexp

from joint_five_letter_splits import HEIGHTS
from joint_moving_retractions import coordinate_indices, moving_crossing, profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--min-power', type=int, default=-12)
    parser.add_argument('--max-power', type=int, default=12)
    parser.add_argument('--max-cases', type=int, default=14400)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.min_power > args.max_power or args.max_cases < 1:
        parser.error('Use an ordered power interval and a positive case budget.')
    names = (Path(__file__).name, 'joint_moving_retractions.py', 'joint_five_letter_splits.py',
             'joint_energy_conjugations.py', 'search_f4_recentered_crossings.py',
             'physical_split_crossings.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp', 'native_weighted_smith_f4.py')
    report = dict(scope='Finite complete-profile sample and finite numerical scale scan; no uniform pressure or factoring claim.',
                  min_power=args.min_power, max_power=args.max_power, scalars=[1, 2],
                  omitted_scalar_three='Frobenius of scalar two has the same binary-source Smith profiles.',
                  records=[], classes=[], source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names})
    scales = np.exp(np.linspace(-7, 1, 65))
    physical = 2*sqrt(2)*logsumexp(-sqrt(2)*scales[:, None]*np.asarray(HEIGHTS)[None, :], axis=1)
    report['scale_grid'] = scales.tolist()
    outside = [i for i in range(25) if i not in coordinate_indices()]
    known, started, best, exhausted = {}, time.monotonic(), None, True
    for inside, extra, side, power, scalar in product(range(16), outside, ('row', 'column'),
                                                    range(args.min_power, args.max_power+1), (1, 2)):
        if len(report['records']) >= args.max_cases:
            exhausted = False
            break
        code = moving_crossing(inside, extra, power, scalar, side)
        data = profiles(code)
        key = tuple(data['zero_orders']), tuple(data['infinity_orders'])
        if key not in known:
            known[key] = len(report['classes'])
            zero, infinity = (np.asarray(data[k]) for k in ('zero_orders', 'infinity_orders'))
            values = logsumexp(-scales[:, None]*zero[None, :], axis=1)+logsumexp(-scales[:, None]*infinity[None, :], axis=1)-physical
            j = int(np.argmin(values))
            target = values+physical
            lower, upper = np.full(65, 1e-7), np.full(65, 64.)
            for _ in range(44):
                middle = (lower+upper)/2
                pressure = target-4*middle*logsumexp(-scales[:, None]*np.asarray(HEIGHTS)[None, :]/middle[:, None], axis=1)
                lower, upper = np.where(pressure > 0, middle, lower), np.where(pressure > 0, upper, middle)
            r = int(np.argmin(upper))
            result = dict(**data, minimum_sampled_critical_pressure=float(values[j]), sampled_scale=float(scales[j]),
                          minimum_sampled_rho=float(upper[r]), rho_sampled_scale=float(scales[r]),
                          representative=dict(inside=inside, outside=extra, side=side, power=power, scalar=scalar))
            report['classes'].append(result)
            if best is None or upper[r] < best:
                best = float(upper[r])
                print(json.dumps(dict(tested=len(report['records'])+1, classes=len(known), best_sampled_rho=best,
                                      representative=result['representative'])), flush=True)
        report['records'].append([inside, extra, side, power, scalar, known[key]])
        if len(report['records']) % 1000 == 0:
            report.update(elapsed_seconds=time.monotonic()-started, finite_list_exhausted=False)
            args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    report.update(completed=True, finite_list_exhausted=exhausted, elapsed_seconds=time.monotonic()-started,
                  minimum_sampled_rho=best)
    args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps(dict(cases=len(report['records']), classes=len(known), exhausted=exhausted, best=best,
                          seconds=report['elapsed_seconds'])))


if __name__ == '__main__':
    main()
