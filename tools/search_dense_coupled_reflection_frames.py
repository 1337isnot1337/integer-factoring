#!/usr/bin/env python3
"""Complete natural-profile screen, then bounded dense-frame gauge searches.

The finite frame catalogue and polynomial coefficients are exact.  Native
valuations, Gibbs roots, and continuous optimization are search evidence.
No numerical result from this script is an exponent certificate.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path
import signal
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from asymmetric_crossing_pressure import AsymmetricPressure
from dense_coupled_reflection_frames import rank_one_terms, transformed_from_terms, digest, save
from native_rectangular_smith_f4 import RectangularOracle, _logsumexp
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize
from search_physical_code_gauges import start_params


class Deadline(Exception):
    pass


def alarm(signum, frame):
    raise Deadline()


def profile_rho(orders, energy, scale):
    n = len(orders) // 2
    constant = (_logsumexp([-scale * v for v in orders[:n]])
                + _logsumexp([-scale * v for v in orders[n:]]))
    low, high = 1e-12, 8.0

    def equation(root):
        return constant - 2 * root * _logsumexp([-scale * value / root for value in energy])

    if equation(low) <= 0:
        raise ArithmeticError('Nonpositive Gibbs root in the natural-profile screen.')
    while equation(high) > 0:
        high *= 2
    for _ in range(65):
        middle = (low + high) / 2
        if equation(middle) > 0:
            low = middle
        else:
            high = middle
    return (low + high) / 4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--refine', type=int, default=48)
    parser.add_argument('--trial-seconds', type=float, default=12)
    parser.add_argument('--iterations', type=int, default=220)
    parser.add_argument('--seed', type=int, default=312132)
    args = parser.parse_args()
    if args.seconds <= 0 or args.refine < 0 or args.trial_seconds <= 0:
        parser.error('Use positive time budgets and a nonnegative refinement count.')
    signal.signal(signal.SIGALRM, alarm)
    catalogue = json.loads(args.catalogue.read_text())
    base = tensor(4)
    _, factors = bell_normalize(base)
    assert catalogue['source_seed_sha256'] == digest(base['values'])
    assert catalogue['full_q4_rank_one_catalogue']
    energy = np.asarray(catalogue['energy'], float)
    natural = lambda scale: start_params(energy * scale, scale)
    unit = natural(1.0)
    baseline_oracle = RectangularOracle(base['values'], 4, 4)
    baseline_orders = baseline_oracle.types(unit)[3]
    fitted = minimize_scalar(lambda scale: profile_rho(baseline_orders, energy, scale),
                             bounds=(.005, 1), method='bounded', options={'xatol': 1e-9})
    baseline_scale, baseline_rho = float(fitted.x), float(fitted.fun)
    assert abs(baseline_oracle.types(natural(baseline_scale))[0] - baseline_rho) < 1e-10
    source_names = (Path(__file__).name, 'dense_coupled_reflection_frames.py',
                    'reflection_laurent_physical_frames.py', 'asymmetric_crossing_pressure.py',
                    'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    report = {'scope': 'complete exact finite catalogue; numerical natural-profile screen and bounded full gauge/type refinements',
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'requested_seconds': args.seconds, 'seed': args.seed,
              'trial_seconds': args.trial_seconds, 'iterations': args.iterations,
              'catalogue_sha256': sha256(args.catalogue.read_bytes()).hexdigest(),
              'baseline_rho': baseline_rho, 'baseline_scale': baseline_scale,
              'baseline_orders': baseline_orders, 'alphabet': 4, 'degree': 24,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in source_names},
              'profiles': [], 'screen': [], 'refinements': [],
              'counters': {'native_profiles': 0, 'resolved_objectives': 0, 'unresolved_objectives': 0}}
    started = time.monotonic()
    profiles = {}
    stopped = False
    for generator in catalogue['records']:
        terms = rank_one_terms(base, generator['u'], generator['v'], factors)
        for candidate in generator['candidates']:
            if time.monotonic() - started >= args.seconds:
                stopped = True
                break
            changed = transformed_from_terms(base, terms, candidate['pole'])
            assert digest(changed['values']) == candidate['seed_sha256']
            oracle = RectangularOracle(changed['values'], 4, 4)
            raw = oracle.evaluate(unit, [.25] * 4, [.25] * 4, want_spectra=True)[2]
            # Integer polynomial valuations and integral natural heights.
            key = tuple(round(value) for value in raw)
            assert max(abs(x - y) for x, y in zip(raw, key)) < 1e-8
            report['counters']['native_profiles'] += 1
            if key not in profiles:
                fitted = minimize_scalar(lambda scale: profile_rho(key, energy, scale),
                                         bounds=(.005, 1), method='bounded',
                                         options={'xatol': 1e-7, 'maxiter': 35})
                profiles[key] = len(report['profiles'])
                report['profiles'].append({'orders': list(key), 'scale': float(fitted.x),
                                           'rho': float(fitted.fun)})
            profile_index = profiles[key]
            report['screen'].append({'generator_index': generator['index'], 'pole': candidate['pole'],
                                     'profile_index': profile_index,
                                     'seed_sha256': candidate['seed_sha256']})
        if len(report['screen']) and generator['index'] % 250 == 0:
            report.update(screen_seconds=time.monotonic() - started,
                          unique_natural_profiles=len(report['profiles']),
                          best_natural_rho=min(x['rho'] for x in report['profiles']))
            save(args.output, report)
            print(json.dumps({'stage': 'screen', 'generators': generator['index'] + 1,
                              'frames': len(report['screen']), 'profiles': len(profiles),
                              'best_rho': report['best_natural_rho']}), flush=True)
        if stopped:
            break
    report.update(screen_seconds=time.monotonic() - started,
                  complete_screen=len(report['screen']) == catalogue['polynomial_frames'],
                  unique_natural_profiles=len(report['profiles']))

    ranked = sorted(report['screen'], key=lambda case: report['profiles'][case['profile_index']]['rho'])
    selected, used_profiles = [], set()
    for case in ranked:
        generator = catalogue['records'][case['generator_index']]
        if sum(bool(value) for value in generator['u'] + generator['v']) == 2:
            continue  # The coordinate family is handled in a separate search.
        if case['profile_index'] in used_profiles:
            continue
        used_profiles.add(case['profile_index'])
        selected.append(case)
        if len(selected) >= args.refine:
            break
    if args.refine == 0:
        selected = []
    rng = np.random.default_rng(args.seed)
    warm, best_rho, best_support = natural(baseline_scale), math.inf, math.inf
    for case in selected:
        remaining = args.seconds - (time.monotonic() - started)
        if remaining <= 0:
            break
        generator = catalogue['records'][case['generator_index']]
        terms = rank_one_terms(base, generator['u'], generator['v'], factors)
        changed = transformed_from_terms(base, terms, case['pole'])
        oracle = AsymmetricPressure(changed, 1)
        initial = natural(report['profiles'][case['profile_index']]['scale'])
        valid = {'balanced': [math.inf, None, None], 'asymmetric': [math.inf, None, None]}
        phase = 'balanced'
        record = dict(case, u=generator['u'], v=generator['v'], optimizations=[])

        def objective(params):
            try:
                if phase == 'balanced':
                    value, gradient, details = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4])
                    value, gradient = value / math.sqrt(2), gradient[:-1] / 2
                else:
                    value, gradient, details = oracle.evaluate(params)
                report['counters']['resolved_objectives'] += 1
                if value < valid[phase][0]:
                    valid[phase][:] = [value, np.asarray(params).copy(), details]
                return value, gradient
            except (ArithmeticError, ValueError):
                report['counters']['unresolved_objectives'] += 1
                return 1e6, np.zeros(len(params))

        begin = time.monotonic()
        signal.setitimer(signal.ITIMER_REAL, min(args.trial_seconds, remaining))
        try:
            for first in (initial, warm, initial * .3,
                          initial + np.r_[rng.normal(0, .1, 24), 0]):
                fit = minimize(objective, first, jac=True, method='SLSQP',
                               bounds=[(-24, 24)] * 24 + [(.0005, 2)],
                               options={'maxiter': args.iterations, 'ftol': 2e-11})
                record['optimizations'].append({'phase': phase, 'iterations': int(fit.nit),
                                                'evaluations': int(fit.nfev), 'success': bool(fit.success)})
            phase = 'asymmetric'
            first = valid['balanced'][1] if valid['balanced'][1] is not None else initial
            for angle in (.3, 1.27):
                fit = minimize(objective, np.r_[first / math.sqrt(2), angle], jac=True,
                               method='SLSQP', bounds=[(-24, 24)] * 24 + [(.0005, 2), (.03, math.pi / 2 - .03)],
                               options={'maxiter': args.iterations, 'ftol': 2e-11})
                record['optimizations'].append({'phase': phase, 'iterations': int(fit.nit),
                                                'evaluations': int(fit.nfev), 'success': bool(fit.success)})
            record['status'] = 'all six optimization loops completed'
        except Deadline:
            record['status'] = 'trial deadline in ' + phase + '; best resolved calls retained'
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        for mode, (value, params, details) in valid.items():
            if params is not None:
                record[mode] = {'value': value, 'params': params.tolist(), 'details': details}
        if valid['balanced'][0] < best_rho:
            best_rho, warm = valid['balanced'][0], valid['balanced'][1].copy()
        best_support = min(best_support, valid['asymmetric'][0])
        record['elapsed_seconds'] = time.monotonic() - begin
        report['refinements'].append(record)
        report.update(elapsed_seconds=time.monotonic() - started, best_refined_rho=best_rho,
                      best_support_value=None if math.isinf(best_support) else best_support)
        save(args.output, report)
        print(json.dumps({'stage': 'refine', 'completed': len(report['refinements']),
                          'best_rho': best_rho, 'best_support': report['best_support_value']}), flush=True)
        if best_rho < 1 / math.sqrt(2) - 1e-6 or best_support < 1 - 1e-6:
            report['stopping_reason'] = 'numerical threshold candidate needs independent exact replay'
            break
    report.update(completed=True, elapsed_seconds=time.monotonic() - started,
                  finished_utc=datetime.now(timezone.utc).isoformat(),
                  stopping_reason=report.get('stopping_reason', 'requested screen and refinement limits reached'))
    save(args.output, report)
    print(json.dumps({key: report[key] for key in ('complete_screen', 'unique_natural_profiles',
                                                   'elapsed_seconds', 'stopping_reason')}))


if __name__ == '__main__':
    main()
