#!/usr/bin/env python3
"""Explore every physical-letter gauge on a selected strict crossing seed."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from relabelled_gaussian_sums import (
    RelabelingPressure, affine_relabeling, charge_projection, gaussian_tensor,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', type=int, default=4)
    parser.add_argument('--family', choices=('gaussian', 'indicator', 'parity', 'swap'), default='gaussian')
    parser.add_argument('--labels', type=int, nargs='*', default=[])
    parser.add_argument('--seconds', type=float, default=600)
    parser.add_argument('--maxfun', type=int, default=250)
    parser.add_argument('--restarts', type=int, default=6)
    parser.add_argument('--height-bound', type=float, default=24)
    parser.add_argument('--seed', type=int, default=786214)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    n, q = args.modes, 1 << args.modes
    base = gaussian_tensor(tuple(range(n)), n)
    if args.family == 'gaussian':
        values = tuple((value, 0) for value in base)
    elif args.family == 'indicator':
        if len(args.labels) != 1:
            raise ValueError('An indicator projection requires one label.')
        values = charge_projection(base, [int(label == args.labels[0]) for label in range(q)])
    elif args.family == 'parity':
        values = charge_projection(base, [label.bit_count() % 2 for label in range(q)])
    else:
        if len(args.labels) != 2:
            raise ValueError('A transposition mixture requires two labels.')
        permutation = list(range(q))
        a, b = args.labels
        permutation[a], permutation[b] = permutation[b], permutation[a]
        values = affine_relabeling(base, permutation)
    model = RelabelingPressure(values, tuple(range(n)), n)
    dimension = 8 * (q - 1) + 1
    started = time.monotonic()
    best = None
    history = []
    failures = 0
    rng = np.random.default_rng(args.seed)

    def save():
        result = {'status': 'exploratory floating gauge search; no exact rate claimed',
                  'modes': n, 'family': args.family, 'labels': args.labels,
                  'queries': model.queries, 'failures': failures,
                  'seconds': time.monotonic() - started, 'best': best,
                  'restarts': history, 'random_seed': args.seed,
                  'height_bound': args.height_bound,
                  'parameter_layout': 'two places, four legs, nonzero physical labels; final log(s)'}
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(result, indent=2) + '\n')
        temporary.replace(args.output)

    def objective(parameters):
        nonlocal best, failures
        if time.monotonic() - started > args.seconds:
            raise TimeoutError
        offsets = np.zeros((2, 4, q))
        offsets[:, :, 1:] = parameters[:-1].reshape((2, 4, q - 1))
        try:
            value, (derivative, scalar_derivative), detail = model.evaluate(parameters[-1], offsets, True)
        except ArithmeticError:
            failures += 1
            return 100.0, np.zeros(dimension)
        if best is None or value < best['rho'] - 1e-10:
            best = {'rho': value, 'parameters': parameters.tolist(), **detail}
            save()
            print(json.dumps({'query': model.queries, 'rho': value,
                              'seconds': time.monotonic() - started}), flush=True)
        return value, np.r_[derivative[:, :, 1:].reshape(-1), scalar_derivative]

    try:
        for theta in (0.3, 0.658022, 1.0):
            objective(np.r_[np.zeros(dimension - 1), theta])
        for restart in range(args.restarts):
            initial = np.asarray(best['parameters'])
            if restart:
                initial[:-1] += rng.normal(0, 0.15 if restart % 3 else 0.7, dimension - 1)
                initial[-1] *= np.exp(rng.normal(0, 0.2))
            run = minimize(objective, initial, jac=True, method='L-BFGS-B',
                           bounds=[(-args.height_bound, args.height_bound)] * (dimension - 1) + [(0.005, 3.0)],
                           options={'maxfun': args.maxfun, 'maxiter': args.maxfun,
                                    'maxls': 20, 'ftol': 1e-13, 'gtol': 1e-7})
            history.append({'restart': restart, 'success': bool(run.success),
                            'message': str(run.message), 'evaluations': int(run.nfev),
                            'rho': float(run.fun)})
            save()
    except TimeoutError:
        history.append({'status': 'time budget reached'})
    finally:
        save()
    print(json.dumps({'best_rho': None if best is None else best['rho'],
                      'queries': model.queries, 'seconds': time.monotonic() - started}), flush=True)


if __name__ == '__main__':
    main()
