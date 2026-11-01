#!/usr/bin/env python3
"""Search complete symbol potentials on actual restricted Pascal tensors."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from pascal_letter_gauges import PascalLetterPressure
from permuted_pascal_kernels import permuted_pascal_kernel
from weighted_pascal_flags import change_flag


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--images', type=int, nargs='+', required=True)
    parser.add_argument('--alphabet', type=int, nargs='+')
    parser.add_argument('--backend', choices=('literal', 'sectors'), default='literal')
    parser.add_argument('--precision', type=int, default=512)
    parser.add_argument('--mode-start', type=Path,
                        help='Lift a retained best mode-gauge record, including its constant flag.')
    parser.add_argument('--restarts', type=int, default=4)
    parser.add_argument('--iterations', type=int, default=300)
    parser.add_argument('--seconds', type=float, default=900)
    parser.add_argument('--seed', type=int, default=20261101)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if min(args.restarts, args.iterations, args.seconds) <= 0:
        parser.error('Search budgets must be positive.')
    seed = permuted_pascal_kernel(args.images)
    mode = None
    if args.mode_start:
        mode = json.loads(args.mode_start.read_text())
        if mode['exponent_images'] != args.images:
            parser.error('The mode-start record belongs to a different exponent permutation.')
        if mode['constant_flag'] is not None:
            seed = change_flag(seed, mode['constant_flag'])
    if args.backend == 'sectors':
        from pascal_coordinate_sectors import SectorPascalLetterPressure
        alphabet = list(range(1 << seed['n'])) if args.alphabet is None else args.alphabet
        pressure = SectorPascalLetterPressure(seed, alphabet, precision=args.precision)
    else:
        pressure = PascalLetterPressure(seed, alphabet=args.alphabet, precision=args.precision)
    if mode is not None:
        natural = pressure.lift_mode_parameters(mode['best']['parameters'])
    else:
        natural = pressure.fitted_natural_start()
    pressure.evaluate(natural)
    baseline = pressure.last.copy()
    best = baseline.copy()
    rng = np.random.default_rng(args.seed)
    control = natural.copy()
    control[:-1] += rng.normal(size=len(control) - 1) * .0001
    _, derivative = pressure.evaluate(control, True)
    error = 0.
    for coordinate in (0, 1, pressure.q, 4 * pressure.q, 8 * pressure.q):
        delta = np.zeros_like(control)
        delta[coordinate] = 1e-7
        finite = (pressure.evaluate(control + delta) - pressure.evaluate(control - delta)) / 2e-7
        error = max(error, abs(float(finite - derivative[coordinate])))
    if error > 1e-4:
        raise ArithmeticError('The independent letter gradient control failed: ' + str(error))
    root = Path(__file__).resolve().parents[1]
    sources = ['tools/search_pascal_letter_gauges.py', 'tools/pascal_letter_gauges.py',
               'tools/permuted_pascal_kernels.py', 'tools/weighted_pascal_flags.py',
               'tools/tensor_pascal_box_crossings.py', 'tools/native_matrix_smith_f4.py',
               'tools/native_matrix_smith_f4.cpp', 'tools/native_rectangular_smith_f4.cpp',
               'tools/native_weighted_smith_f4.py']
    if args.backend == 'sectors':
        sources += ['tools/pascal_coordinate_sectors.py']
    report = {'status': 'exploratory floating search', 'exponent_images': args.images,
              'alphabet': pressure.alphabet, 'constant_flag': seed.get('constant_flag'),
              'backend': args.backend, 'precision': args.precision,
              'random_seed': args.seed, 'baseline': baseline, 'best': best, 'runs': [],
              'gradient_max_error': error, 'rejected_queries': 0, 'budget_exhausted': False,
              'scope': 'Finite numerical optimization of complete original-polynomial profiles '
              'with arbitrary symbol potentials and both physical pair costs. No exact pressure '
              'certificate, universal exclusion, or stronger factoring exponent is claimed.',
              'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest()
                                for name in sources}}
    if hasattr(pressure, 'construction'):
        report['construction'] = pressure.construction
    if args.mode_start:
        report['mode_start_sha256'] = sha256(args.mode_start.read_bytes()).hexdigest()
    start = time.monotonic()

    def save():
        report.update(best=best, queries=pressure.queries, seconds=time.monotonic() - start)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(report, indent=2) + '\n')
        temporary.replace(args.output)

    save()
    print(json.dumps({'baseline_rho': baseline['rho'], 'letters': pressure.q,
                      'gradient_max_error': error}), flush=True)
    for trial in range(args.restarts):
        if time.monotonic() - start >= args.seconds:
            report['budget_exhausted'] = True
            break
        scale = natural[-1]
        initial = np.r_[natural[:-1] / (scale * seed['outer']), math.log(scale)]
        if trial:
            initial[:-1] += rng.normal(size=len(initial) - 1) * (.03 if trial % 2 else .3)

        def objective(parameters):
            nonlocal best
            if time.monotonic() - start >= args.seconds:
                raise TimeoutError('The registered search budget is exhausted.')
            current_scale = math.exp(parameters[-1])
            current = np.r_[current_scale * seed['outer'] * parameters[:-1], current_scale]
            try:
                value, derivative = pressure.evaluate(current, True)
            except ArithmeticError:
                report['rejected_queries'] += 1
                return 2., np.zeros_like(parameters)
            if value < best['rho'] - 1e-12:
                best = pressure.last.copy()
                save()
            transformed = np.r_[derivative[:-1] * current_scale * seed['outer'],
                                derivative @ current]
            return value, transformed

        try:
            fitted = minimize(objective, initial, jac=True, method='SLSQP',
                              bounds=[(-12., 12.)] * (len(initial) - 1)
                              + [(math.log(.0001), math.log(2.))],
                              options={'maxiter': args.iterations, 'ftol': 1e-11})
            report['runs'].append({'trial': trial, 'rho': float(fitted.fun),
                                  'iterations': int(fitted.nit), 'success': bool(fitted.success),
                                  'message': str(fitted.message)})
        except TimeoutError:
            report['budget_exhausted'] = True
            break
        save()
        print(json.dumps({'trial': trial, 'best_rho': best['rho'],
                          'queries': pressure.queries, 'seconds': report['seconds']}), flush=True)
    save()


if __name__ == '__main__':
    main()
