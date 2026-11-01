#!/usr/bin/env python3
"""Bounded full-letter gauge search inside the exact rational compound cone.

Profiles and Gibbs-optimized costs are numerical. Every queried gauge
is rounded and repaired in the exact integer cone; a proposed pressure
improvement still needs independent exact order and cost certificates.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from math import exp
from pathlib import Path
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import hstack

from hot_code_charge_profiles import hot_energy
from hot_compound_pressure import evaluate
from search_fused_mixed_weighted import SearchDeadline, atomic_write
from search_hot_downsets import pressure
from weighted_hot_code_oracle import HotWeightedOracle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=24)
    parser.add_argument('--modes', type=int, default=10)
    parser.add_argument('--seconds', type=float, default=1200)
    parser.add_argument('--seed', type=int, default=861732)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 2 <= args.size <= 1 << args.modes or not 1 <= args.modes <= 16:
        raise ValueError('Use a nontrivial supported hot physical code.')
    rng = np.random.default_rng(args.seed)
    labels = sorted(range(1 << args.modes), key=lambda value: (hot_energy(value), value))[:args.size]
    oracle = HotWeightedOracle(labels)
    active = oracle.active_heights
    a, b = oracle.cone(1)
    constraint = hstack((a[:, active], -b[:, None]), format='csr').astype(np.int64)
    denominator = 10 ** 10
    occupation = np.r_[[labels[index - oracle.offsets[leg]].bit_count()
                         for leg in range(4) for index in active
                         if oracle.offsets[leg] <= index < oracle.offsets[leg + 1]], 0]
    assert len(occupation) == len(active) + 1
    zero = np.zeros(oracle.height_count)
    _, _, _, baseline_orders = oracle.log_profile_with_gradient(zero, 1)
    initial = pressure({'orders_w': baseline_orders.tolist(), 'energies_x': oracle.energies_x.tolist(),
                        'energies_y': oracle.energies_y.tolist()})
    current = np.r_[np.zeros(len(active)), initial['scale']]
    names = (Path(__file__).name, 'weighted_hot_code_oracle.py', 'hot_compound_pressure.py',
             'hot_code_charge_profiles.py', 'hot_code_charge_profiles.cpp',
             'selected_gaussian_coefficients.cpp', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'arguments': {key: str(value) if isinstance(value, Path) else value
                                            for key, value in vars(args).items()},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'started_utc': datetime.now(timezone.utc).isoformat(), 'letters': labels,
              'constraint_count': constraint.shape[0], 'rational_denominator': denominator,
              'baseline': initial, 'queries': [], 'iterations': [], 'completed': False,
              'unresolved': 0, 'linear_program_failures': 0}
    began, deadline = time.monotonic(), time.monotonic() + args.seconds

    def query(values):
        if time.monotonic() >= deadline:
            raise SearchDeadline()
        integer = np.rint(np.asarray(values) * denominator).astype(np.int64)
        residue = constraint @ integer
        minimum = int(np.min(residue)) if len(residue) else 0
        if minimum < 0:
            integer += ((-minimum + 1) // 2) * occupation
            residue = constraint @ integer
        if len(residue) and int(np.min(residue)) < 0:
            raise ArithmeticError('The rational cone repair failed.')
        values = integer.astype(float) / denominator
        heights = np.zeros(oracle.height_count)
        heights[active] = values[:-1]
        started = time.monotonic()
        try:
            rho, dh, ds, details = evaluate(oracle, heights, float(values[-1]))
        except (ArithmeticError, ValueError):
            report['unresolved'] += 1
            raise
        gradient = np.r_[dh[active], ds]
        record = {'rho': rho, 'seconds': time.monotonic() - started,
                  'minimum_cone_numerator': int(np.min(residue)) if len(residue) else 0,
                  'parameters_sha256': sha256(integer.tobytes()).hexdigest()}
        report['queries'].append(record)
        if 'best' not in report or rho < report['best']['rho']:
            report['best'] = details | {'rational_parameters': integer.tolist(),
                                        'rational_denominator': denominator,
                                        'minimum_cone_numerator': record['minimum_cone_numerator']}
            print(json.dumps({'best_rho': rho, 'queries': len(report['queries']),
                              'elapsed_seconds': time.monotonic() - began}), flush=True)
        return values, rho, gradient

    radius, restart, restart_needed = .15, 0, False
    try:
        current, current_rho, gradient = query(current)
        while time.monotonic() < deadline:
            if restart_needed:
                restart += 1
                heights = []
                amplitude = .02 * 2 ** (restart % 4)
                for _ in range(4):
                    potentials = rng.uniform(0, amplitude, args.size)
                    heights.extend(sum(potentials[j] for j, subset in enumerate(labels)
                                       if subset and subset & value == subset) for value in labels)
                candidate = np.r_[np.asarray(heights)[active], initial['scale'] * exp(rng.normal(0, .1))]
                current, current_rho, gradient = query(candidate)
                radius, restart_needed = .15, False
            scale_radius = radius / (1 + max(oracle.energies_x))
            bounds = [(max(-radius, -32 - value), min(radius, 32 - value)) for value in current[:-1]]
            bounds.append((max(-scale_radius, .005 - current[-1]), min(scale_radius, 1 - current[-1])))
            result = linprog(gradient, A_ub=-constraint, b_ub=constraint @ current,
                             bounds=bounds, method='highs')
            if not result.success:
                report['linear_program_failures'] += 1
                restart_needed = True
                continue
            predicted = float(gradient @ result.x)
            accepted, trials = False, 0
            if predicted < -1e-10:
                for exponent in range(9):
                    trials += 1
                    values, rho, next_gradient = query(current + result.x * 2.0 ** -exponent)
                    if rho < current_rho - 1e-11:
                        current, current_rho, gradient = values, rho, next_gradient
                        accepted = True
                        break
            report['iterations'].append({'restart': restart, 'rho': current_rho, 'predicted_change': predicted,
                                         'radius': radius, 'accepted': accepted, 'line_search_trials': trials})
            if accepted:
                radius = min(.6, radius * 1.15)
            else:
                radius *= .35
                if radius < .0003 or predicted >= -1e-10:
                    restart_needed = True
            report['elapsed_seconds'] = time.monotonic() - began
            atomic_write(args.output, report)
    except SearchDeadline:
        pass
    report.update(completed=True, elapsed_seconds=time.monotonic() - began,
                  budget_completed=time.monotonic() >= deadline,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
