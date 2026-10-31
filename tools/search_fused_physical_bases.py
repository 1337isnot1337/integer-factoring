#!/usr/bin/env python3
"""Budgeted physical-basis and full-symbol mixed-place searches for fused seeds."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from f16_identity_operator_profiles import MUL, inverse_matrix
from fused_six_vertex_crossings import power
from fused_physical_basis_oracle import FusedPhysicalBasisOracle, identity, quantum_pascal, vandermonde
from search_fused_mixed_weighted import MixedPressure, SearchDeadline, atomic_write


class ScaledHeightPressure(MixedPressure):
    """Use h=t*H so line searches retain complete native series profiles.

    The optimized bounded H coordinates describe arbitrary physical-letter
    heights within the declared search window. Actual h and their complete
    type shifts are still what the mathematical criterion evaluates.
    """
    def initial(self, baseline, rng, start):
        result = super().initial(baseline, rng, start)
        scale = result[self.scale_index]
        for cut in self.height_slices:
            for current in cut:
                result[current] /= scale
        return result

    def evaluate(self, values, smoothing=0.01, record=True):
        actual = np.asarray(values, dtype=float).copy()
        scale = float(values[self.scale_index])
        for cut in self.height_slices:
            for current in cut:
                actual[current] *= scale
        previous = self.best
        value, derivative = super().evaluate(actual, smoothing, record)
        scale_derivative = derivative[self.scale_index]
        for cut in self.height_slices:
            for current in cut:
                scale_derivative += float(np.asarray(values[current]) @ derivative[current])
                derivative[current] *= scale
        derivative[self.scale_index] = scale_derivative
        if self.best is not previous:
            self.best['optimizer_parameters'] = np.asarray(values).tolist()
            self.best['height_parameterization'] = 'Actual physical heights h equal outer log scale t times bounded optimized coordinates H.'
        return value, derivative


def structured_pairs(q, theta):
    unit = identity(q)
    named = [('identity', unit), ('ordinary_pascal', quantum_pascal(q, 1)),
             ('quantum_pascal', quantum_pascal(q, power(theta, 2))),
             ('inverse_quantum_pascal', inverse_matrix(quantum_pascal(q, power(theta, 2)))),
             ('evaluation', vandermonde(range(q))),
             ('inverse_evaluation', inverse_matrix(vandermonde(range(q)))),
             ('geometric_evaluation', vandermonde([power(theta, i) for i in range(q)])),
             ('inverse_geometric_evaluation', inverse_matrix(vandermonde([power(theta, i) for i in range(q)])))]
    result, seen = [], set()
    for name, value in named:
        for suffix, first, second in [('common', value, value), ('x_only', value, unit), ('y_only', unit, value)]:
            key = tuple(map(tuple, first)), tuple(map(tuple, second))
            if key not in seen:
                result.append((name + '_' + suffix, key[0], key[1]))
                seen.add(key)
    return result


def mutate(first, second, q, field_size, rng, steps):
    pair = [[list(row) for row in first], [list(row) for row in second]]
    for _ in range(steps):
        side = int(rng.integers(2))
        i, j = map(int, rng.choice(q, size=2, replace=False))
        scalar = int(rng.integers(1, field_size))
        pair[side][i] = [a ^ MUL[scalar][b] for a, b in zip(pair[side][i], pair[side][j])]
    return tuple(tuple(row) for row in pair[0]), tuple(tuple(row) for row in pair[1])


def compact_best(best):
    return {key: value for key, value in best.items() if key != 'place_records'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('baseline', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--seed', type=int, default=312014)
    parser.add_argument('--maxiter', type=int, default=140)
    parser.add_argument('--maxfun', type=int, default=1400)
    parser.add_argument('--alphabets', type=int, nargs='+', default=[3, 4])
    parser.add_argument('--trials', type=int, default=1000000)
    args = parser.parse_args()
    started = time.monotonic()
    deadline = started + args.seconds
    baseline_bytes = args.baseline.read_bytes()
    source_seeds = json.loads(baseline_bytes)['seeds']
    baselines = []
    for seed in source_seeds:
        if seed['seed']['alphabet'] in args.alphabets:
            baselines.append({'m': seed['seed']['m'], 'theta': seed['seed']['theta'],
                              'alphabet': seed['seed']['alphabet'],
                              'partitions': sorted(seed['completed_partitions'], key=lambda row: row['rho'])[:3]})
    if not baselines:
        raise ValueError('No requested alphabet is present in the baseline.')
    rng = np.random.default_rng(args.seed)
    families = {seed['alphabet']: structured_pairs(seed['alphabet'], seed['theta']) for seed in baselines}
    winners, records = {}, []
    names = (Path(__file__).name, 'fused_physical_basis_oracle.py', 'search_fused_mixed_weighted.py',
             'compressed_fused_six_vertex.py', 'f16_weighted_local_oracle.py',
             'fused_six_vertex_crossings.py', 'f16_identity_operator_profiles.py', 'f4_lattice_linear.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp',
             'native_weighted_smith_f4.py')
    directory = Path(__file__).parent
    result = {'scope': 'Finite floating full-symbol basis/gauge search, with exact field arithmetic and complete type/degree costs; no global optimum or factoring improvement certified.',
              'arguments': vars(args) | {'baseline': str(args.baseline), 'output': str(args.output)},
              'baseline_sha256': sha256(baseline_bytes).hexdigest(),
              'baseline_partitions': baselines,
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names},
              'records': records, 'winners': winners, 'budget_completed': False}
    for trial in range(args.trials):
        if time.monotonic() >= deadline:
            break
        seed = baselines[trial % len(baselines)]
        q, m, theta = seed['alphabet'], seed['m'], seed['theta']
        round_index = trial // len(baselines)
        family = families[q]
        if round_index < len(family):
            basis_name, first, second = family[round_index]
        else:
            if str(q) in winners and round_index % 5:
                parent = winners[str(q)]
                first, second = parent['first_basis'], parent['second_basis']
                basis_name = 'winner_elementary_row_mutation'
            else:
                basis_name, first, second = family[int(rng.integers(len(family)))]
                basis_name += '_restart_mutation'
            first, second = mutate(first, second, q, 4 if q == 3 else 16, rng,
                                   1 + (round_index % (2 * q)))
        partition_index = round_index % len(seed['partitions'])
        baseline = seed['partitions'][partition_index]
        trial_started = time.monotonic()
        oracle = FusedPhysicalBasisOracle(m, theta, first, second)
        pressure = ScaledHeightPressure(oracle, baseline, 'full', False, True, deadline)
        initial = pressure.initial(baseline, rng, round_index % 4)
        if str(q) in winners and round_index >= len(family) and round_index % 3:
            previous = winners[str(q)]['best']
            if len(previous['optimizer_parameters']) == pressure.dimension:
                initial = np.asarray(previous['optimizer_parameters'], dtype=float)
        outcomes, complete, interruption = [], True, None
        try:
            pressure.evaluate(initial)
            for smoothing in (0.035, 0.006, 0.001):
                optimized = minimize(lambda values: pressure.evaluate(values, smoothing), initial,
                                     jac=True, method='L-BFGS-B', bounds=pressure.bounds(),
                                     options={'maxiter': args.maxiter, 'maxfun': args.maxfun,
                                              'ftol': 1e-11, 'gtol': 3e-7, 'maxls': 30})
                outcomes.append({'smoothing': smoothing, 'optimizer_success': bool(optimized.success),
                                 'optimizer_message': str(optimized.message), 'iterations': int(optimized.nit),
                                 'function_evaluations': int(optimized.nfev)})
                initial = np.asarray(pressure.best['optimizer_parameters'])
        except SearchDeadline:
            complete, interruption = False, 'Declared wall-clock budget reached during this trial.'
        except (ArithmeticError, ValueError) as error:
            complete, interruption = False, f'{type(error).__name__}: {error}'
        record = {'trial': trial, 'm': m, 'theta': theta, 'alphabet': q,
                  'basis_name': basis_name, 'first_basis': first, 'second_basis': second,
                  'baseline_partition_index': partition_index, 'completed_trial': complete,
                  'queries': pressure.queries, 'elapsed_seconds': time.monotonic() - trial_started,
                  'outcomes': outcomes, 'best': compact_best(pressure.best) if pressure.best else None,
                  'interruption': interruption}
        records.append(record)
        if pressure.best and (str(q) not in winners or pressure.best['rho'] < winners[str(q)]['best']['rho']):
            winners[str(q)] = record | {'best': pressure.best}
        result.update({'elapsed_seconds': time.monotonic() - started,
                       'completed_trials': sum(item['completed_trial'] for item in records),
                       'attempted_trials': len(records), 'budget_completed': time.monotonic() >= deadline})
        atomic_write(args.output, result)
        print(json.dumps({'trial': trial, 'q': q, 'basis': basis_name, 'complete': complete,
                          'rho': pressure.best['rho'] if pressure.best else None,
                          'best_rho_for_seed': winners.get(str(q), {}).get('best', {}).get('rho'),
                          'queries': pressure.queries, 'elapsed': result['elapsed_seconds']}), flush=True)
        if time.monotonic() >= deadline:
            break
    result.update({'elapsed_seconds': time.monotonic() - started,
                   'completed_trials': sum(item['completed_trial'] for item in records),
                   'attempted_trials': len(records), 'budget_completed': time.monotonic() >= deadline,
                   'trial_limit_completed': len(records) >= args.trials})
    atomic_write(args.output, result)


if __name__ == '__main__':
    main()
