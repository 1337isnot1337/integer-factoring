#!/usr/bin/env python3
"""Complete q4 constant companions f=w^-k+c of dense rank-one frames."""
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
from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import add, bell_normalize, scale
from search_dense_coupled_reflection_frames import profile_rho, Deadline, alarm
from search_physical_code_gauges import start_params


def translate_terms(terms, constant):
    """T(f+c) by exact binomial coefficients in characteristic two."""
    powers = [1]
    for _ in range(4):
        powers.append(F4_PRODUCT[powers[-1]][constant])
    output = [{} for _ in range(5)]
    for old in range(5):
        for new in range(old + 1):
            if old & new != new:  # Lucas parity for binomial(old,new).
                continue
            coefficient = powers[old - new]
            if not coefficient:
                continue
            for index, value in terms[old].items():
                changed = add(output[new].get(index, (0, 0)), scale(value, coefficient))
                if changed != (0, 0):
                    output[new][index] = changed
                elif index in output[new]:
                    del output[new][index]
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path, default=Path('evidence/dense-rank-one-coupled-frames.json'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--refine', type=int, default=80)
    parser.add_argument('--trial-seconds', type=float, default=12)
    parser.add_argument('--iterations', type=int, default=200)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, alarm)
    catalogue = json.loads(args.catalogue.read_text())
    base = tensor(4)
    _, factors = bell_normalize(base)
    assert catalogue['source_seed_sha256'] == digest(base['values'])
    energy = np.asarray([0, 12, 4, 16], float)
    natural = lambda s: start_params(energy * s, s)
    unit = natural(1.0)
    base_oracle = RectangularOracle(base['values'], 4, 4)
    base_orders = base_oracle.evaluate(unit, [.25] * 4, [.25] * 4, want_spectra=True)[2]
    fitted = minimize_scalar(lambda s: profile_rho(base_orders, energy, s), bounds=(.005, 1),
                             method='bounded', options={'xatol': 1e-9})
    baseline_scale, baseline_rho = float(fitted.x), float(fitted.fun)
    names = (Path(__file__).name, 'dense_coupled_reflection_frames.py',
             'reflection_laurent_physical_frames.py', 'search_dense_coupled_reflection_frames.py',
             'asymmetric_crossing_pressure.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp')
    report = {'scope': 'complete exact constant-companion catalogue; numerical complete-profile search',
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'requested_seconds': args.seconds, 'requested_refinements': args.refine,
              'trial_seconds': args.trial_seconds, 'iterations': args.iterations,
              'catalogue_sha256': sha256(args.catalogue.read_bytes()).hexdigest(),
              'baseline_rho': baseline_rho, 'baseline_scale': baseline_scale,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'profiles': [], 'screen': [], 'refinements': [],
              'counters': {'proposed': 0, 'duplicate_tensors': 0, 'resolved_objectives': 0,
                           'unresolved_objectives': 0}}
    started = time.monotonic()
    profiles, seen = {}, set()
    complete = True
    for generator in catalogue['records']:
        if time.monotonic() - started > args.seconds * .45:
            complete = False
            break
        terms = rank_one_terms(base, generator['u'], generator['v'], factors)
        for constant in (1, 2, 3):
            translated = translate_terms(terms, constant)
            for pole in generator['pole_orders']:
                report['counters']['proposed'] += 1
                changed = transformed_from_terms(base, translated, pole)
                seed_digest = digest(changed['values'])
                if seed_digest in seen:
                    report['counters']['duplicate_tensors'] += 1
                    continue
                seen.add(seed_digest)
                oracle = RectangularOracle(changed['values'], 4, 4)
                raw = oracle.evaluate(unit, [.25] * 4, [.25] * 4, want_spectra=True)[2]
                key = tuple(round(value) for value in raw)
                assert max(abs(a - b) for a, b in zip(raw, key)) < 1e-8
                if key not in profiles:
                    fitted = minimize_scalar(lambda s: profile_rho(key, energy, s), bounds=(.005, 1),
                                             method='bounded', options={'xatol': 1e-7, 'maxiter': 35})
                    profiles[key] = len(report['profiles'])
                    report['profiles'].append({'orders': list(key), 'rho': float(fitted.fun),
                                               'scale': float(fitted.x)})
                report['screen'].append({'generator_index': generator['index'], 'pole': pole,
                                         'constant': constant, 'seed_sha256': seed_digest,
                                         'profile_index': profiles[key]})
        if generator['index'] % 250 == 0:
            report.update(elapsed_seconds=time.monotonic() - started,
                          best_natural_rho=min(record['rho'] for record in report['profiles']))
            save(args.output, report)
            print(json.dumps({'generators': generator['index'] + 1, 'screen': len(report['screen']),
                              'profiles': len(profiles), 'best_rho': report['best_natural_rho']}), flush=True)
    report.update(complete_companion_screen=complete, screen_seconds=time.monotonic() - started)
    ranked = sorted(report['screen'], key=lambda row: report['profiles'][row['profile_index']]['rho'])
    selected, used = [], set()
    for case in ranked:
        if case['profile_index'] not in used:
            used.add(case['profile_index'])
            selected.append(case)
        if len(selected) >= args.refine:
            break
    if args.refine == 0:
        selected = []
    best_rho, best_support = math.inf, math.inf
    warm = natural(baseline_scale)
    for candidate in selected:
        remaining = args.seconds - (time.monotonic() - started)
        if remaining <= 0:
            break
        generator = catalogue['records'][candidate['generator_index']]
        terms = rank_one_terms(base, generator['u'], generator['v'], factors)
        changed = transformed_from_terms(base, translate_terms(terms, candidate['constant']), candidate['pole'])
        oracle = AsymmetricPressure(changed, 1)
        initial = natural(report['profiles'][candidate['profile_index']]['scale'])
        phase = 'balanced'
        valid = {'balanced': [math.inf, None, None], 'asymmetric': [math.inf, None, None]}
        record = dict(candidate, u=generator['u'], v=generator['v'], optimizations=[])

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
            for first in (initial, warm, initial * .3):
                fit = minimize(objective, first, jac=True, method='SLSQP',
                               bounds=[(-24, 24)] * 24 + [(.0005, 2)],
                               options={'maxiter': args.iterations, 'ftol': 2e-11})
                record['optimizations'].append({'phase': phase, 'iterations': int(fit.nit),
                                                'evaluations': int(fit.nfev), 'success': bool(fit.success)})
            phase = 'asymmetric'
            first = valid['balanced'][1] if valid['balanced'][1] is not None else initial
            for angle in (.3, 1.27):
                fit = minimize(objective, np.r_[first / math.sqrt(2), angle], jac=True, method='SLSQP',
                               bounds=[(-24, 24)] * 24 + [(.0005, 2), (.03, math.pi / 2 - .03)],
                               options={'maxiter': args.iterations, 'ftol': 2e-11})
                record['optimizations'].append({'phase': phase, 'iterations': int(fit.nit),
                                                'evaluations': int(fit.nfev), 'success': bool(fit.success)})
            record['status'] = 'all five optimization loops completed'
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
        report.update(elapsed_seconds=time.monotonic() - started, best_rho=best_rho,
                      best_support_value=None if math.isinf(best_support) else best_support)
        save(args.output, report)
        print(json.dumps({'refinements': len(report['refinements']), 'best_rho': best_rho,
                          'best_support': report['best_support_value']}), flush=True)
        if best_rho < 1 / math.sqrt(2) - 1e-6 or best_support < 1 - 1e-6:
            report['stopping_reason'] = 'numerical threshold candidate requires exact independent replay'
            break
    report.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat(),
                  elapsed_seconds=time.monotonic() - started,
                  stopping_reason=report.get('stopping_reason', 'declared finite screen/refinement limits reached'))
    save(args.output, report)
    print(json.dumps({key: report[key] for key in ('completed', 'elapsed_seconds', 'stopping_reason')}))


if __name__ == '__main__':
    main()
