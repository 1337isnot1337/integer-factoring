#!/usr/bin/env python3
"""Complete-prefix search of energy-preserving joint coordinate grids."""
import argparse
from hashlib import sha256
import json
from math import cos, exp, log, pi, sin
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import logsumexp

from joint_grid_compressions import coordinate_maps, crossing, energy_grids, profiles
from screen_unequal_coherent_hot_codes import support


def numerical_cost(hx, hy, data):
    zero, infinity = (np.asarray(data[key], dtype=float) for key in ('zero_orders', 'infinity_orders'))
    first, second = np.asarray(hx, dtype=float), np.asarray(hy, dtype=float)
    def objective(s):
        t = exp(s)
        target = float(logsumexp(-t*zero)+logsumexp(-t*infinity))
        return brentq(lambda rho: target-2*rho*float(logsumexp(-t*first/rho)+logsumexp(-t*second/rho)),
                      1e-7, 2, xtol=1e-12)
    balanced = minimize_scalar(objective, bounds=(log(.005), log(2)), method='bounded')
    sides, ex, ey = [data['infinity_orders'], data['zero_orders']], [2*h for h in hx], [2*h for h in hy]
    trials = []
    for angle in (pi/6, pi/4, pi/3):
        fit = minimize(lambda point: support(sides, ex, ey, exp(point[0]), cos(point[1]), sin(point[1]), 1),
                       [float(balanced.x)-log(2)/2, angle], method='L-BFGS-B',
                       bounds=((log(.005), log(2)), (.05, pi/2-.05)))
        trials.append(dict(support=float(fit.fun), scale=exp(float(fit.x[0])),
                           angle=float(fit.x[1]), completed=bool(fit.success)))
    return dict(balanced_rho=float(balanced.fun), balanced_scale=exp(float(balanced.x)),
                balanced_completed=bool(balanced.success), asymmetric_trials=trials,
                best_support=min(trial['support'] for trial in trials))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--x', type=int, required=True)
    parser.add_argument('--y', type=int, required=True)
    parser.add_argument('--max-cases', type=int, default=2000)
    parser.add_argument('--seconds', type=float, default=900)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.x < 2 or args.y < 2 or args.x*args.y > 16 or args.max_cases < 1 or args.seconds <= 0:
        parser.error('Use nontrivial grid dimensions with product at most sixteen and positive budgets.')
    names = (Path(__file__).name, 'joint_grid_compressions.py', 'joint_height_conjugations.py',
             'check_triangular_auxiliary_row_fiber.py', 'weighted_smith_oracle_f4.py',
             'check_joint_energy_conjugations.py', 'screen_unequal_coherent_hot_codes.py')
    report = dict(q_x=args.x, q_y=args.y, records=[], classes=[],
                  source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
                  scope='Complete exact profiles for the declared finite prefix; numerical scale/support fits only, no global optimum or factoring claim.')
    started, known, best, exhausted = time.monotonic(), {}, None, True
    for hx, hy in energy_grids(args.x, args.y):
        for mapping in coordinate_maps(hx, hy):
            if len(report['records']) >= args.max_cases or time.monotonic()-started >= args.seconds:
                exhausted = False
                break
            data = profiles(crossing(hx, hy, mapping))
            key = (hx, hy, tuple(data['zero_orders']), tuple(data['infinity_orders']))
            if key not in known:
                known[key] = len(report['classes'])
                cost = numerical_cost(hx, hy, data)
                report['classes'].append(dict(heights_x=list(hx), heights_y=list(hy), **data, **cost))
                if best is None or cost['best_support'] < best:
                    best = cost['best_support']
                    print(json.dumps(dict(tested=len(report['records'])+1, classes=len(known), best=best,
                                          balanced_rho=cost['balanced_rho'], hx=hx, hy=hy)), flush=True)
            report['records'].append(dict(mapping=list(mapping), profile_class=known[key]))
        if not exhausted:
            break
        report.update(elapsed_seconds=time.monotonic()-started, finite_list_exhausted=False)
        args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    report.update(completed=True, finite_list_exhausted=exhausted, elapsed_seconds=time.monotonic()-started,
                  best_numerical_support=best)
    args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps(dict(tested=len(report['records']), classes=len(known), exhausted=exhausted, best=best)))


if __name__ == '__main__':
    main()
