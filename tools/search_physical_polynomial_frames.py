#!/usr/bin/env python3
"""Numerical complete-profile search through exact polynomial physical frames."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import signal
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from asymmetric_crossing_pressure import AsymmetricPressure
from physical_polynomial_frames import apply_word
from reflection_dual_rail_gaussian import tensor
from search_physical_code_gauges import start_params


class Deadline(Exception):
    pass


def alarm(signum, frame):
    raise Deadline()


def save(path, report):
    head = json.dumps({k: v for k, v in report.items() if k != 'cases'}, indent=2)
    data = head[:-2] + ',\n  "cases": [\n'
    data += ',\n'.join(json.dumps(r, separators=(',', ':')) for r in report['cases'])
    data += '\n  ]\n}\n'
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(data)
    temporary.replace(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seconds', type=float, default=7200)
    p.add_argument('--trial-seconds', type=float, default=15)
    p.add_argument('--iterations', type=int, default=160)
    p.add_argument('--seed', type=int, default=312114)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if min(args.seconds, args.trial_seconds, args.iterations) <= 0:
        p.error('positive budgets are required')
    signal.signal(signal.SIGALRM, alarm)
    rng = random.Random(args.seed)
    noise = np.random.default_rng(args.seed)
    base = tensor(4)
    energy = np.asarray([0, 12, 4, 16])
    natural = lambda s: start_params(energy * s, s)
    initial_oracle = AsymmetricPressure(base, 1)

    def balanced(oracle, params):
        value, gradient, detail = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4])
        return value / math.sqrt(2), gradient[:-1] / 2, detail

    scale = minimize_scalar(lambda s: balanced(initial_oracle, natural(s))[0],
                            bounds=(.005, 1), method='bounded', options={'xatol': 1e-10}).x
    natural_start = natural(scale)
    alphabet_pairs = [(i, j) for i in range(4) for j in range(4) if i != j]
    names = (Path(__file__).name, 'physical_polynomial_frames.py', 'asymmetric_crossing_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'reflection_dual_rail_gaussian.py')
    report = {'schema_version': 1, 'scope': 'numerical search; every polynomial frame and endpoint is exact',
              'seed': args.seed, 'requested_seconds': args.seconds, 'trial_seconds': args.trial_seconds,
              'iteration_limit': args.iterations, 'alphabet': 4, 'base_degree': base['degree'],
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'counters': {'calls': 0, 'unresolved': 0}, 'cases': []}
    started, index = time.monotonic(), 0
    best_balanced = best_asymmetric = None
    elites, seen = [], set()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    while time.monotonic() - started < args.seconds:
        if index == 0:
            word, label = [], 'identity'
        elif index <= 144:
            offset = index - 1
            target, source = alphabet_pairs[offset % 12]
            pair = (offset // 12) % 2
            degree = 1 + (offset // 24) % 6
            scalar = 1 + (offset // 48) % 3
            terms = [[0, scalar], [degree, scalar]]
            word = [{'pair': pair, 'target': target, 'source': source, 'terms': terms}]
            label = 'single-endpoint-trivial-shear'
        else:
            word = list(rng.choice(elites)['word']) if elites and index % 3 else []
            for _ in range(rng.randint(1, 2 if word else 4)):
                target, source = rng.choice(alphabet_pairs)
                degree = rng.randint(1, 8)
                scalar = rng.randint(1, 3)
                terms = [[degree, scalar]]
                if rng.randrange(3):
                    terms.insert(0, [0, scalar if rng.randrange(2) else rng.randint(1, 3)])
                word.append({'pair': rng.randrange(2), 'target': target, 'source': source, 'terms': terms})
            label = 'composed-polynomial-shears'
        index += 1
        seed = apply_word(base, word)
        if seed['degree'] >= 160:
            continue
        digest = sha256(json.dumps(seed['values'], separators=(',', ':')).encode()).hexdigest()
        if digest in seen:
            continue
        seen.add(digest)
        cut = 0 if len(report['cases']) % 9 == 8 else 1
        oracle = AsymmetricPressure(seed, cut)
        calls_before = dict(report['counters'])
        case_start = time.monotonic()
        candidates = [natural_start.copy(), natural_start * .3]
        if best_balanced is not None:
            candidates.append(np.asarray(best_balanced['balanced']['params'], float))
        if len(report['cases']) % 4 == 3:
            candidates.append(natural_start + np.r_[noise.normal(0, .15, 24), 0.0])
        record = {'generation': index - 1, 'kind': label, 'word': word, 'degree': seed['degree'],
                  'seed_sha256': digest, 'cut_at_infinity': cut, 'endpoint_entries': 256, 'optimizations': []}
        valid_balanced, valid_asym = [math.inf, None, None], [math.inf, None, None]
        phase = 'balanced'

        def objective(params):
            report['counters']['calls'] += 1
            try:
                if phase == 'balanced':
                    value, gradient, details = balanced(oracle, params)
                    target = valid_balanced
                else:
                    value, gradient, details = oracle.evaluate(params)
                    target = valid_asym
                if value < target[0]:
                    target[:] = [value, np.asarray(params, float).copy(), details]
                return value, gradient
            except (ArithmeticError, ValueError):
                report['counters']['unresolved'] += 1
                return 1e6, np.zeros(len(params))

        budget = min(args.trial_seconds, max(.01, args.seconds - (case_start - started)))
        signal.setitimer(signal.ITIMER_REAL, budget)
        try:
            scored = [(objective(x)[0], j) for j, x in enumerate(candidates)]
            start = candidates[min(scored)[1]]
            fit = minimize(objective, start, jac=True, method='SLSQP',
                           bounds=[(-24, 24)] * 24 + [(.0005, 2)],
                           options={'maxiter': args.iterations, 'ftol': 2e-11})
            record['optimizations'].append({'phase': phase, 'iterations': int(fit.nit),
                                            'evaluations': int(fit.nfev), 'success': bool(fit.success),
                                            'message': str(fit.message)})
            phase = 'asymmetric'
            start = (valid_balanced[1] if valid_balanced[1] is not None else start) / math.sqrt(2)
            start = np.r_[start, .35 if len(report['cases']) % 2 else 1.22]
            fit = minimize(objective, start, jac=True, method='SLSQP',
                           bounds=[(-24, 24)] * 24 + [(.0005, 2), (.03, math.pi / 2 - .03)],
                           options={'maxiter': args.iterations, 'ftol': 2e-11})
            record['optimizations'].append({'phase': phase, 'iterations': int(fit.nit),
                                            'evaluations': int(fit.nfev), 'success': bool(fit.success),
                                            'message': str(fit.message)})
            record['status'] = 'optimization loops completed'
        except Deadline:
            record['status'] = 'trial deadline in ' + phase + '; best resolved calls retained'
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        if valid_balanced[1] is not None:
            record['balanced'] = {'rho': valid_balanced[0], 'params': valid_balanced[1].tolist(),
                                  'details': valid_balanced[2]}
            if best_balanced is None or valid_balanced[0] < best_balanced['balanced']['rho']:
                best_balanced = record
            elites.append(record)
            elites.sort(key=lambda r: r['balanced']['rho'])
            del elites[12:]
        if valid_asym[1] is not None:
            record['asymmetric'] = {'support_value': valid_asym[0], 'params': valid_asym[1].tolist(),
                                    'details': valid_asym[2]}
            if best_asymmetric is None or valid_asym[0] < best_asymmetric['asymmetric']['support_value']:
                best_asymmetric = record
        record.update(elapsed_seconds=time.monotonic() - case_start,
                      calls={k: report['counters'][k] - calls_before[k] for k in calls_before})
        report['cases'].append(record)
        report.update(completed=False, elapsed_seconds=time.monotonic() - started,
                      completed_cases=len(report['cases']),
                      best_rho=None if best_balanced is None else best_balanced['balanced']['rho'],
                      best_support_value=None if best_asymmetric is None else best_asymmetric['asymmetric']['support_value'])
        save(args.output, report)
        print(json.dumps({k: report[k] for k in ('completed_cases', 'elapsed_seconds', 'best_rho', 'best_support_value')}), flush=True)
        if ((best_balanced is not None and report['best_rho'] < 1 / math.sqrt(2) - 1e-6)
                or (best_asymmetric is not None and report['best_support_value'] < 1 - 1e-6)):
            report['stopping_reason'] = 'numerical threshold candidate needs exact independent replay'
            break
    report.update(completed=True, elapsed_seconds=time.monotonic() - started,
                  stopping_reason=report.get('stopping_reason', 'requested time budget completed'))
    save(args.output, report)


if __name__ == '__main__':
    main()
