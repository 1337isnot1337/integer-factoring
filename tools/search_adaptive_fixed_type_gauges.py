#!/usr/bin/env python3
"""Explore a local-gauge catalogue on one fixed pair of physical types.

The numerical local rate at coefficient fraction d is optimized separately
at each end. The resulting gauges then bound every fraction, so their lower
envelope is a valid catalogue. Numerical maxima and optimizers do not certify
an all-width crossing; a proposed gain needs exact finite-block rank bounds.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp, softmax

from native_rectangular_smith_f4 import RectangularOracle
from reflection_dual_rail_gaussian import physical_bit_energies, tensor
from search_fused_mixed_weighted import SearchDeadline, atomic_write
from search_physical_code_gauges import start_params


class LocalFixedTypeRates:
    """Complete local Smith orders and chamber derivatives at both ends."""

    def __init__(self, seed, type_x, type_y):
        self.seed = seed
        self.qx, self.qy = seed['q_x'], seed['q_u']
        self.n, self.degree = self.qx * self.qy, seed['degree']
        self.oracle = RectangularOracle(seed['values'], self.qx, self.qy)
        self.alphabets = (self.qx, self.qy, self.qy, self.qx)
        self.offsets = np.cumsum([0] + [q - 1 for q in self.alphabets[:-1]])
        self.stride = sum(q - 1 for q in self.alphabets)
        self.types = tuple(np.asarray(p, float) for p in (type_x, type_y))
        if any(len(p) != q or min(p) <= 0 or abs(sum(p) - 1) > 1e-10
               for p, q in zip(self.types, (self.qx, self.qy))):
            raise ValueError('Use positive, normalized physical types.')
        self.marginals = np.concatenate([self.types[0 if leg in (0, 3) else 1][1:]
                                         for leg in range(4)])
        self.entropy = -sum(float(p @ np.log(p)) for p in self.types)

    def profile(self, side, variables):
        """Variables are log-scaled heights and the positive local scale."""
        variables = np.asarray(variables, float)
        if side not in (0, 1) or len(variables) != self.stride + 1:
            raise ValueError('Use side 0=infinity or 1=zero and all local heights.')
        parameters = np.zeros(2 * self.stride + 1)
        parameters[side * self.stride:(side + 1) * self.stride] = variables[:-1]
        parameters[-1] = variables[-1]
        data = self.oracle.prefixes(parameters, cut=1)
        orders = np.asarray(data['orders'][side * self.n:(side + 1) * self.n], float)
        features = np.zeros((self.n, self.stride + 1))
        previous = 0
        for index in range(self.n):
            at = side * self.n + index
            row, column = data['rows'][at], data['cols'][at]
            x, second = divmod(row, self.qy)
            third, v = divmod(column, self.qx)
            letters = (x, third, second, v) if side == 0 else (x, second, third, v)
            for leg, letter in enumerate(letters):
                if letter:
                    features[index, self.offsets[leg] + letter - 1] += 1
            cumulative = data['base_prefix'][at]
            features[index, -1] = cumulative - previous
            previous = cumulative
        if side == 0:
            # The native end uses z^-1 without the homogeneous z^D scalar.
            orders += self.degree * variables[-1]
            features[:, -1] += self.degree
        if max(abs(features @ variables - orders)) > 2e-7:
            raise ArithmeticError('Complete local prefix features do not replay the orders.')
        return orders, features

    def evaluate(self, side, fraction, variables):
        orders, features = self.profile(side, variables)
        probabilities = softmax(-orders)
        mean = float(self.marginals @ variables[:-1])
        value = float(logsumexp(-orders)) + mean + float(variables[-1]) * fraction
        gradient = -(probabilities @ features)
        gradient[:-1] += self.marginals
        gradient[-1] += fraction
        return value, gradient

    def catalogue_record(self, side, variables):
        orders, _ = self.profile(side, variables)
        scale = float(variables[-1])
        mean = float(self.marginals @ variables[:-1])
        return {'side': side, 'heights': (np.asarray(variables[:-1]) / scale).tolist(),
                'orders': (orders / scale).tolist(), 'mean_shift': mean / scale,
                'effective_orders': ((orders - mean) / scale).tolist(),
                'parameters': list(map(float, variables))}


def local_curve(profile, fractions, entropy, temperatures):
    """Direct Chernoff envelope; retains exact-boundary positive jet costs."""
    orders = np.asarray(profile['effective_orders'], float)
    intercept = logsumexp(-temperatures[:, None] * orders[None, :], axis=1)
    rates = np.min(intercept[:, None] + temperatures[:, None] * fractions[None, :], axis=0)
    rates = np.minimum(entropy, rates)
    # A negative exponential rate can only bound an eventually zero rank.
    # Keep it as -infinity, while the cold boundary itself retains log of
    # the minimum-channel multiplicity. This is search guidance only.
    rates[fractions < min(orders) - 1e-9] = -np.inf
    tied = np.isclose(fractions, min(orders), atol=1e-9, rtol=0)
    rates[tied] = min(entropy, math.log(sum(abs(orders - min(orders)) < 1e-9)))
    return rates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, choices=(4, 8), default=4)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--fractions', type=int, default=17)
    parser.add_argument('--restarts', type=int, default=2)
    parser.add_argument('--seconds', type=float, default=600)
    parser.add_argument('--trial-seconds', type=float, default=30)
    parser.add_argument('--iterations', type=int, default=160)
    parser.add_argument('--seed', type=int, default=312731)
    args = parser.parse_args()
    base = tensor(args.interval)
    q = base['q_x']
    energies = np.asarray([sum(v for bit, v in enumerate(physical_bit_energies(args.interval))
                               if letter >> bit & 1) for letter in range(q)], float)
    native = RectangularOracle(base['values'], q, q)
    start = lambda scale: start_params(energies * scale, scale)
    optimum = minimize_scalar(lambda scale: native.types(start(scale))[0],
                              bounds=(.001, 1), method='bounded', options={'xatol': 1e-11})
    baseline = start(float(optimum.x))
    baseline_rho, px, py, _, _, _ = native.types(baseline)
    oracle = LocalFixedTypeRates(base, px, py)
    local_baselines = [np.r_[baseline[side * oracle.stride:(side + 1) * oracle.stride], baseline[-1]]
                       for side in range(2)]
    initial_profiles = [oracle.catalogue_record(side, local_baselines[side]) for side in range(2)]
    fractions = np.linspace(0, base['degree'], 2049)
    temperatures = np.r_[0., np.geomspace(1e-6, 80, 2049)]
    curves = [local_curve(p, fractions, oracle.entropy, temperatures) for p in initial_profiles]
    finite = np.isfinite(curves[0] + curves[1][::-1])
    baseline_grid_rho = float(np.max((curves[0] + curves[1][::-1])[finite]) / (2 * oracle.entropy))
    if abs(baseline_grid_rho - baseline_rho) > 1e-4:
        raise ArithmeticError('The initial adaptive envelope does not match the two-end control.')
    # Verify representative derivatives away from Smith chamber walls.
    noise = np.random.default_rng(args.seed)
    controls = []
    for side in range(2):
        point = local_baselines[side].copy()
        point[:-1] += noise.normal(0, .002, len(point) - 1)
        fraction = base['degree'] * .43
        value, derivative = oracle.evaluate(side, fraction, point)
        errors = []
        for coordinate in range(len(point)):
            delta = np.zeros_like(point)
            delta[coordinate] = 2e-7
            difference = (oracle.evaluate(side, fraction, point + delta)[0]
                          - oracle.evaluate(side, fraction, point - delta)[0]) / 4e-7
            errors.append(abs(float(derivative[coordinate]) - difference))
        if max(errors) > 3e-6:
            raise ArithmeticError('A local fixed-type derivative control failed.')
        controls.append({'side': side, 'coordinates': len(point), 'maximum_error': max(errors)})
    names = (Path(__file__).name, 'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'reflection_dual_rail_gaussian.py', 'search_physical_code_gauges.py',
             'search_fused_mixed_weighted.py')
    report = {'scope': __doc__, 'arguments': vars(args) | {'output': str(args.output)},
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'degree': base['degree'], 'alphabet': q, 'type_x': px, 'type_y': py,
              'entropy_sum': oracle.entropy, 'baseline_rho': baseline_rho,
              'baseline_grid_rho': baseline_grid_rho, 'gradient_controls': controls,
              'profiles': [[initial_profiles[0]], [initial_profiles[1]]], 'records': [],
              'completed_trials': 0, 'completed': False}
    began, deadline = time.monotonic(), time.monotonic() + args.seconds
    local_targets = np.linspace(.02, .98, args.fractions) * base['degree']
    targets = [(side, float(fraction), restart) for restart in range(args.restarts)
               for fraction in sorted(local_targets, key=lambda x: abs(x - base['degree'] / 2))
               for side in range(2)]
    best_rho = baseline_grid_rho
    for side, fraction, restart in targets:
        if time.monotonic() >= deadline:
            break
        point = local_baselines[side].copy()
        if restart:
            point[:-1] += noise.normal(0, .15, len(point) - 1)
            point[-1] *= noise.uniform(.5, 1.5)
        best = [math.inf, point.copy()]
        unresolved = 0
        trial_deadline = min(deadline, time.monotonic() + args.trial_seconds)

        def objective(values):
            nonlocal unresolved
            if time.monotonic() >= trial_deadline:
                raise SearchDeadline()
            try:
                value, gradient = oracle.evaluate(side, fraction, values)
            except (ArithmeticError, ValueError):
                unresolved += 1
                return 1e6, np.zeros(len(values))
            if value < best[0]:
                best[:] = [value, np.asarray(values).copy()]
            return value, gradient

        completed, fit = True, None
        try:
            objective(point)
            fit = minimize(objective, point, jac=True, method='SLSQP',
                           bounds=[(-48, 48)] * oracle.stride + [(1e-5, 4)],
                           options={'maxiter': args.iterations, 'ftol': 2e-10})
        except SearchDeadline:
            completed = False
        profile = oracle.catalogue_record(side, best[1])
        report['profiles'][side].append(profile)
        curves[side] = np.minimum(curves[side], local_curve(profile, fractions, oracle.entropy, temperatures))
        total = curves[0] + curves[1][::-1]
        index = int(np.argmax(total))
        rho = float(total[index] / (2 * oracle.entropy))
        best_rho = min(best_rho, rho)
        record = {'side': side, 'target_fraction': fraction, 'restart': restart,
                  'local_bound': best[0], 'profile_index': len(report['profiles'][side]) - 1,
                  'optimization_completed': completed,
                  'optimizer_success': bool(fit.success) if fit is not None else False,
                  'optimizer_message': str(fit.message) if fit is not None else 'Declared local or global budget reached.',
                  'evaluations': int(fit.nfev) if fit is not None else None,
                  'iterations': int(fit.nit) if fit is not None else None, 'unresolved_calls': unresolved,
                  'catalogue_grid_rho': rho, 'critical_fraction_infinity': float(fractions[index])}
        report['records'].append(record)
        report.update(elapsed_seconds=time.monotonic() - began, attempted_trials=len(report['records']),
                      completed_trials=sum(record['optimization_completed'] for record in report['records']),
                      best_grid_rho=best_rho)
        atomic_write(args.output, report)
        print(json.dumps(record | {'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(completed=True, completed_plan=len(report['records']) == len(targets),
                  elapsed_seconds=time.monotonic() - began,
                  stopping_reason='All planned local targets completed.' if len(report['records']) == len(targets)
                  else 'Declared time budget reached between local targets.')
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
