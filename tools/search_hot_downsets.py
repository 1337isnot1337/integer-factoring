#!/usr/bin/env python3
"""Bounded natural-pressure search over arbitrary hot physical downsets.

Every proposal remains closed under deleting bits. Complete charge-block
profiles determine the recorded costs. Optimized costs are numerical;
the finite search does not certify global optimality or a factoring gain.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from math import exp, log
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

from hot_code_charge_profiles import hot_energy, natural_profiles
from search_fused_mixed_weighted import atomic_write


def pressure(profile):
    orders = np.asarray(profile['orders_w'], dtype=float)
    first, second = (np.asarray(profile[key], dtype=float) for key in ('energies_x', 'energies_y'))

    def objective(log_scale):
        scale = exp(log_scale)
        target = 2 * float(logsumexp(-scale * orders))
        equation = lambda root: target - root * (float(logsumexp(-scale * first / root))
                                                   + float(logsumexp(-scale * second / root)))
        return brentq(equation, 1e-10, 4, xtol=1e-13) / 2

    fit = minimize_scalar(objective, bounds=(log(.005), log(1.5)), method='bounded',
                          options={'xatol': 1e-9, 'maxiter': 70})
    return {'rho': float(fit.fun), 'scale': exp(float(fit.x)),
            'completed': bool(fit.success), 'evaluations': int(fit.nfev)}


def neighbors(code, modes):
    code = set(code)
    removable = [value for value in code if value and not any(
        (value | (1 << bit)) in code for bit in range(modes) if not value >> bit & 1)]
    candidates = {value | (1 << bit) for value in code for bit in range(modes)
                  if not value >> bit & 1}
    addable = [value for value in candidates - code if all(
        (value ^ (1 << bit)) in code for bit in range(modes) if value >> bit & 1)]
    result = []
    for removed in removable:
        for added in addable:
            if removed & added == removed:
                continue
            result.append(tuple(sorted((code - {removed}) | {added})))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=24)
    parser.add_argument('--modes', type=int, default=10)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--seed', type=int, default=938241)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 2 <= args.size <= min(128, 1 << args.modes) or not 1 <= args.modes <= 16:
        raise ValueError('Use a supported finite physical code and at most sixteen modes.')
    rng = np.random.default_rng(args.seed)
    source_names = (Path(__file__).name, 'hot_code_charge_profiles.py', 'hot_code_charge_profiles.cpp',
                    'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'arguments': {key: str(value) if isinstance(value, Path) else value
                                            for key, value in vars(args).items()},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in source_names},
              'started_utc': datetime.now(timezone.utc).isoformat(), 'records': [],
              'duplicate_proposals': 0, 'completed': False, 'unresolved': 0}
    labels = sorted(range(1 << args.modes), key=lambda value: (hot_energy(value), value))
    current = tuple(sorted(labels[:args.size]))
    visited, best, current_rho = {}, None, float('inf')
    began, deadline = time.monotonic(), time.monotonic() + args.seconds
    queue = [current]
    while time.monotonic() < deadline:
        if queue:
            proposal = queue.pop(0)
        else:
            available = neighbors(current, args.modes)
            rng.shuffle(available)
            proposal = next((value for value in available if value not in visited), None)
            if proposal is None:
                report['duplicate_proposals'] += len(available)
                if best is None:
                    break
                # A fresh positive additive energy still produces a downset.
                weights = np.asarray([8 * a + 4 for a in range(args.modes)]) * np.exp(rng.normal(0, .35, args.modes))
                ordered = sorted(range(1 << args.modes), key=lambda value: (
                    sum(weights[a] for a in range(args.modes) if value >> a & 1), value))
                current = tuple(sorted(ordered[:args.size]))
                proposal = current
                current_rho = float('inf')
        if proposal in visited:
            report['duplicate_proposals'] += 1
            if best is not None:
                current = tuple(best['letters'])
                current_rho = best['pressure']['rho']
            continue
        start = time.monotonic()
        record = {'index': len(report['records']), 'letters': list(proposal), 'completed': False}
        try:
            profile = natural_profiles(proposal)
            result = pressure(profile)
            record.update(pressure=result, orders_w=profile['orders_w'], energies=profile['energies_x'],
                          block_sizes=[block['size'] for block in profile['blocks']],
                          precision_t=[block['precision_t'] for block in profile['blocks']],
                          completed=True)
            value = result['rho']
            if best is None or value < best['pressure']['rho']:
                best = record
                report['best_record_index'] = record['index']
                print(json.dumps({'new_best': value, 'letters': proposal,
                                  'elapsed_seconds': time.monotonic() - began}), flush=True)
            temperature = .00015 * (0.2 + 0.8 * max(0, (deadline - time.monotonic()) / args.seconds))
            if value < current_rho or rng.random() < exp(min(0, (current_rho - value) / temperature)):
                current, current_rho = proposal, value
            if record['index'] == 0:
                queue = neighbors(current, args.modes)
                rng.shuffle(queue)
        except (ArithmeticError, ValueError) as exception:
            record['error'] = f'{type(exception).__name__}: {exception}'
            report['unresolved'] += 1
        record['seconds'] = time.monotonic() - start
        visited[proposal] = record['index']
        report['records'].append(record)
        report.update(elapsed_seconds=time.monotonic() - began, attempted_trials=len(report['records']),
                      completed_trials=sum(row['completed'] for row in report['records']))
        atomic_write(args.output, report)
    report.update(completed=True, elapsed_seconds=time.monotonic() - began,
                  finished_utc=datetime.now(timezone.utc).isoformat(),
                  budget_completed=time.monotonic() >= deadline)
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
