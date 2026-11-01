#!/usr/bin/env python3
"""Bounded sixteen-letter physical-basis search outside the Gaussian frames.

Exact F4 gates act on two of the four physical modes. Complete local
profiles and endpoints are used; floating pressure is exploratory only.
"""
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
from physical_split_crossings import identity_split
from reflection_dual_rail_gaussian import tensor, physical_bit_energies, profile_factors
from search_physical_code_gauges import start_params
from search_reflection_physical_bases import change_basis, generate_basis_pair, row_shear, canonical_basis


class TrialDeadline(Exception):
    pass


def alarm_handler(signum, frame):
    raise TrialDeadline('per-trial wall-time limit')


def lift(gate, modes):
    """Lift one exact four-letter gate; bit0 of the gate is modes[0]."""
    other = 15 ^ sum(1 << b for b in modes)
    out = [[0] * 16 for _ in range(16)]
    for i in range(16):
        for j in range(16):
            if i & other == j & other:
                a = ((i >> modes[0]) & 1) + 2 * ((i >> modes[1]) & 1)
                b = ((j >> modes[0]) & 1) + 2 * ((j >> modes[1]) & 1)
                out[i][j] = gate[a][b]
    return out


def save(path, data):
    head = {k: v for k, v in data.items() if k != 'cases'}
    result = json.dumps(head, indent=2)[:-2] + ',\n  "cases": [\n'
    result += ',\n'.join(json.dumps(x, separators=(',', ':')) for x in data['cases'])
    result += '\n  ]\n}\n'
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(result)
    tmp.replace(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seconds', type=float, default=7200)
    p.add_argument('--trial-seconds', type=float, default=120)
    p.add_argument('--iterations', type=int, default=100)
    p.add_argument('--seed', type=int, default=312104)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if min(args.seconds, args.trial_seconds, args.iterations) <= 0:
        p.error('all budgets must be positive')
    signal.signal(signal.SIGALRM, alarm_handler)
    rng = random.Random(args.seed)
    noise = np.random.default_rng(args.seed)
    base = tensor(8)
    energies = np.asarray([sum(e for b, e in enumerate(physical_bit_energies(8)) if s >> b & 1)
                           for s in range(16)], float)
    natural = lambda scale: start_params(energies * scale, scale)
    baseline = AsymmetricPressure(base, 1)

    def pressure(oracle, params):
        rate, gradient, details = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4])
        return rate / math.sqrt(2), gradient[:-1] / 2, details

    exact = baseline.oracle.prefixes(natural(1.0), 1)['orders']
    subset = sorted(sum(v for i, v in enumerate(profile_factors(8)) if mask >> i & 1)
                    for mask in range(256))
    assert exact == [x - 32 for x in subset] + [x + 32 for x in subset]
    scalar = minimize_scalar(lambda x: pressure(baseline, natural(x))[0],
                             bounds=(.005, .3), method='bounded',
                             options={'xatol': 1e-10})
    seed_start = natural(float(scalar.x))
    baseline_value = pressure(baseline, seed_start)[0]
    names = (Path(__file__).name, 'asymmetric_crossing_pressure.py',
             'search_reflection_physical_bases.py', 'physical_split_crossings.py',
             'reflection_dual_rail_gaussian.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp')
    report = {'schema_version': 1, 'scope': 'bounded numerical search; no certified rate or optimality claim',
              'alphabet': 16, 'interval': 8, 'degree': base['degree'], 'seed': args.seed,
              'requested_seconds': args.seconds, 'trial_seconds': args.trial_seconds,
              'iteration_limit': args.iterations, 'baseline_rho': baseline_value,
              'baseline_scale': float(scalar.x),
              'controls': {'exact_natural_profile_orders': 512, 'endpoint_entries_per_seed': 65536},
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'counters': {'objective_calls': 0, 'unresolved_calls': 0}, 'cases': []}
    started = time.monotonic()
    pairs = [(0, 1), (0, 3), (1, 2), (2, 3), (0, 2), (1, 3)]
    seen, elites, best = set(), [], None
    index = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    while time.monotonic() - started < args.seconds:
        if index == 0:
            label, left, right, pair = 'identity', identity_split(16), identity_split(16), None
        else:
            # Independent four-letter gates include odd/even mixtures and are
            # not restricted to exterior powers or one-particle frames.
            gate_index = (1 + ((index // 4) % 6) if index % 4 == 0
                          else 7 + ((index - 1) // len(pairs)))
            pair = pairs[(index - 1) % len(pairs)]
            kind, gx, gy = generate_basis_pair(gate_index, 4, rng, [])
            right_pair = pairs[(index + 2) % len(pairs)] if index % 3 == 0 else pair
            left, right, label = lift(gx, pair), lift(gy, right_pair), 'two-mode-' + kind
            if index > 36 and elites and index % 3:
                ref = rng.choice(elites)
                left, right = ref['basis_x'], ref['basis_y']
                for _ in range(rng.randint(1, 2)):
                    i, j = rng.sample(range(16), 2)
                    if rng.randrange(2):
                        left = row_shear(left, i, j, rng.choice((1, 2, 3)))
                    else:
                        right = row_shear(right, i, j, rng.choice((1, 2, 3)))
                label = 'full-letter-neighbor'
        index += 1
        key = canonical_basis(left), canonical_basis(right)
        if key in seen:
            continue
        seen.add(key)
        seed = change_basis(base, left, right)
        oracle = AsymmetricPressure(seed, 1)
        valid = [math.inf, None, None]
        calls_before = dict(report['counters'])
        begin = time.monotonic()

        def objective(params):
            report['counters']['objective_calls'] += 1
            try:
                value, gradient, details = pressure(oracle, params)
                if value < valid[0]:
                    valid[:] = [value, np.asarray(params, float).copy(), details]
                return value, gradient
            except (ArithmeticError, ValueError):
                report['counters']['unresolved_calls'] += 1
                return 1e6, np.zeros(len(params))

        initial = seed_start.copy()
        if best is not None and len(report['cases']) % 3 == 1:
            initial = np.asarray(best['best']['params'], float)
        if len(report['cases']) % 3 == 2:
            initial[:-1] += noise.normal(0, .06, len(initial) - 1)
        status, result_data = 'completed optimization', {}
        signal.setitimer(signal.ITIMER_REAL, min(args.trial_seconds, max(.01, args.seconds - (begin - started))))
        try:
            candidates = [initial, seed_start * .2, seed_start * .5]
            chosen = min((objective(start)[0], k) for k, start in enumerate(candidates))
            initial_value, chosen_index = chosen
            initial = candidates[chosen_index]
            result = minimize(objective, initial, jac=True, method='SLSQP',
                              bounds=[(-16, 16)] * (len(initial) - 1) + [(.0005, .6)],
                              options={'maxiter': args.iterations, 'ftol': 1e-10})
            result_data = {'initial_rho': initial_value, 'success': bool(result.success),
                           'message': str(result.message), 'iterations': int(result.nit),
                           'evaluations': int(result.nfev)}
        except TrialDeadline:
            status = 'trial time limit; best resolved call retained'
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        record = {'generation': index - 1, 'kind': label, 'mode_pair': pair,
                  'basis_x': left, 'basis_y': right, 'status': status,
                  'optimizer': result_data, 'elapsed_seconds': time.monotonic() - begin,
                  'counters': {k: report['counters'][k] - calls_before[k] for k in calls_before}}
        if valid[1] is not None:
            record['best'] = {'rho': valid[0], 'params': valid[1].tolist(), 'details': valid[2]}
            elites.append(record)
            elites.sort(key=lambda r: r['best']['rho'])
            del elites[10:]
            if best is None or valid[0] < best['best']['rho']:
                best = record
        report['cases'].append(record)
        report.update(completed=False, elapsed_seconds=time.monotonic() - started,
                      completed_cases=len(report['cases']), best_rho=None if best is None else best['best']['rho'])
        save(args.output, report)
        print(json.dumps({'cases': len(report['cases']), 'best_rho': report['best_rho'],
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)
        if best is not None and best['best']['rho'] < 1 / math.sqrt(2) - 1e-6:
            report['stopping_reason'] = 'numerical threshold candidate requires exact independent replay'
            break
    report.update(completed=True, elapsed_seconds=time.monotonic() - started,
                  stopping_reason=report.get('stopping_reason', 'requested time budget completed'))
    save(args.output, report)


if __name__ == '__main__':
    main()
