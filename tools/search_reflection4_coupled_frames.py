#!/usr/bin/env python3
"""Complete finite coupled-pole catalogue and bounded full-gauge pressure search."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import signal
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from asymmetric_crossing_pressure import AsymmetricPressure
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import (
    bell_normalize, bell_restore, coupled_terms, coupled_admissible_poles,
    coupled_transformed,
)
from search_physical_code_gauges import start_params


class Deadline(Exception):
    pass


def alarm(signum, frame):
    raise Deadline()


def save(path, report):
    head = json.dumps({k: v for k, v in report.items() if k != 'cases'}, indent=2)
    value = head[:-2] + ',\n  "cases": [\n'
    value += ',\n'.join(json.dumps(r, separators=(',', ':')) for r in report['cases'])
    value += '\n  ]\n}\n'
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(value)
    temp.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--trial-seconds', type=float, default=12)
    parser.add_argument('--iterations', type=int, default=220)
    parser.add_argument('--seed', type=int, default=312121)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, alarm)
    base = tensor(4)
    bell, factors = bell_normalize(base)
    catalogue = []
    candidates = [({'kind': 'identity'}, bell)]
    for target in range(4):
        for source in range(4):
            if target == source:
                continue
            for scalar in (1, 2, 3):
                terms = coupled_terms(bell, target, source, scalar)
                allowed = coupled_admissible_poles(terms, base['degree'])
                catalogue.append({'target': target, 'source': source, 'scalar': scalar,
                                  'tested_k': [1, base['degree']], 'allowed_k': allowed})
                for k in allowed:
                    changed = coupled_transformed(bell, target, source, -k, scalar=scalar, terms=terms)
                    assert changed['degree'] <= base['degree']
                    candidates.append(({'kind': 'coupled-negative-frame', 'target': target,
                                        'source': source, 'scalar': scalar, 'pole_order': k}, changed))
    rng = np.random.default_rng(args.seed)
    natural = lambda s: start_params(np.asarray([0, 12, 4, 16]) * s, s)
    baseline = AsymmetricPressure(base, 1)

    def balanced(oracle, params):
        value, gradient, detail = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4])
        return value / math.sqrt(2), gradient[:-1] / 2, detail

    scale = minimize_scalar(lambda s: balanced(baseline, natural(s))[0], bounds=(.005, 1),
                            method='bounded', options={'xatol': 1e-10}).x
    initial = natural(scale)
    names = (Path(__file__).name, 'reflection_laurent_physical_frames.py',
             'asymmetric_crossing_pressure.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp')
    report = {'schema_version': 1, 'scope': 'exact finite admissibility; numerical bounded pressure search',
              'seed': args.seed, 'requested_seconds': args.seconds, 'trial_seconds': args.trial_seconds,
              'iteration_limit': args.iterations, 'alphabet': 4, 'degree': base['degree'],
              'admissibility_catalogue': catalogue, 'admissible_frames': len(candidates) - 1,
              'planned_coordinate_cases': 2 * len(candidates),
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'counters': {'calls': 0, 'unresolved': 0}, 'cases': []}
    started = time.monotonic()
    best_rho = best_support = math.inf
    warm = initial.copy()
    finished = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for gate, seed_bell in candidates:
        for coordinates in ('original', 'bell'):
            if time.monotonic() - started >= args.seconds:
                finished = False
                break
            seed = bell_restore(seed_bell, factors) if coordinates == 'original' else seed_bell
            oracle = AsymmetricPressure(seed, 1)
            record = dict(gate, coordinates=coordinates, degree=seed['degree'],
                          monomial_content=seed_bell.get('monomial_content', 0),
                          seed_sha256=sha256(json.dumps(seed['values'], separators=(',', ':')).encode()).hexdigest(),
                          endpoint_entries=256, optimizations=[])
            valid = {'balanced': [math.inf, None, None], 'asymmetric': [math.inf, None, None]}
            phase = 'balanced'
            begin = time.monotonic()

            def objective(params):
                report['counters']['calls'] += 1
                try:
                    value, gradient, details = (balanced(oracle, params) if phase == 'balanced'
                                                else oracle.evaluate(params))
                    if value < valid[phase][0]:
                        valid[phase][:] = [value, np.asarray(params, float).copy(), details]
                    return value, gradient
                except (ArithmeticError, ValueError):
                    report['counters']['unresolved'] += 1
                    return 1e6, np.zeros(len(params))

            budget = min(args.trial_seconds, max(.01, args.seconds - (begin - started)))
            signal.setitimer(signal.ITIMER_REAL, budget)
            try:
                starts = [initial.copy(), warm.copy(), initial * .25,
                          initial + np.r_[rng.normal(0, .12, 24), 0.0]]
                for start in starts:
                    fit = minimize(objective, start, jac=True, method='SLSQP',
                                   bounds=[(-24, 24)] * 24 + [(.0005, 2)],
                                   options={'maxiter': args.iterations, 'ftol': 2e-11})
                    record['optimizations'].append({'phase': phase, 'iterations': int(fit.nit),
                                                    'evaluations': int(fit.nfev), 'success': bool(fit.success),
                                                    'message': str(fit.message)})
                phase = 'asymmetric'
                selected = valid['balanced'][1]
                if selected is None:
                    selected = initial
                for angle in (.3, 1.27):
                    fit = minimize(objective, np.r_[selected / math.sqrt(2), angle], jac=True,
                                   method='SLSQP', bounds=[(-24, 24)] * 24 + [(.0005, 2), (.03, math.pi / 2 - .03)],
                                   options={'maxiter': args.iterations, 'ftol': 2e-11})
                    record['optimizations'].append({'phase': phase, 'iterations': int(fit.nit),
                                                    'evaluations': int(fit.nfev), 'success': bool(fit.success),
                                                    'message': str(fit.message)})
                record['status'] = 'all six optimization loops completed'
            except Deadline:
                record['status'] = 'trial deadline in ' + phase + '; best resolved calls retained'
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            for mode in ('balanced', 'asymmetric'):
                value, params, details = valid[mode]
                if params is not None:
                    record[mode] = {'value': value, 'params': params.tolist(), 'details': details}
            if valid['balanced'][1] is not None and valid['balanced'][0] < best_rho:
                best_rho, warm = valid['balanced'][0], valid['balanced'][1].copy()
            best_support = min(best_support, valid['asymmetric'][0])
            record['elapsed_seconds'] = time.monotonic() - begin
            report['cases'].append(record)
            report.update(completed=False, elapsed_seconds=time.monotonic() - started,
                          completed_cases=len(report['cases']), best_rho=best_rho,
                          best_support_value=None if math.isinf(best_support) else best_support)
            save(args.output, report)
            print(json.dumps({k: report[k] for k in ('completed_cases', 'elapsed_seconds', 'best_rho', 'best_support_value')}), flush=True)
            if best_rho < 1 / math.sqrt(2) - 1e-6 or best_support < 1 - 1e-6:
                finished = False
                report['stopping_reason'] = 'numerical threshold candidate requires exact independent replay'
                break
        if not finished:
            break
    report.update(completed=True, complete_catalogue_optimization=finished,
                  elapsed_seconds=time.monotonic() - started,
                  stopping_reason=report.get('stopping_reason', 'all coordinate cases completed' if finished else 'time budget completed'))
    save(args.output, report)


if __name__ == '__main__':
    main()
