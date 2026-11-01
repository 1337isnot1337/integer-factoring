#!/usr/bin/env python3
"""Budgeted exploratory mixed-place fused searches with complete type costs.

Each completed trial is saved atomically. Both place sets remain disjoint,
conjugate roots share a gauge within their declared group, and one pair of
physical types is used at every place. Floating results are not certificates.
"""
from __future__ import annotations

import argparse
import ast
from hashlib import sha256
import json
from math import log
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import softmax, logsumexp

from fused_weighted_place_profiles import FusedPlaceOracle


class SearchDeadline(Exception):
    pass


def point_key(label):
    if label is None:
        return None
    if isinstance(label, int):
        return (label, 1)
    if not label.startswith('root:'):
        raise ValueError('Unrecognized exact closed-place label.')
    return tuple(ast.literal_eval(label.split(':')[1]))


def grouped_places(places, allocations, include_inactive=False):
    groups = {}
    for label, allocation in zip(places, allocations):
        if allocation <= 1e-7 and not include_inactive:
            continue
        key = point_key(label)
        if key not in groups:
            groups[key] = {'point': key, 'labels': [], 'allocation': 0.0}
        groups[key]['labels'].append(label)
        groups[key]['allocation'] += max(allocation, 0.0)
    if include_inactive:
        for group in groups.values():
            group['allocation'] = 0.92 * group['allocation'] + 0.08 / len(groups)
    total = sum(group['allocation'] for group in groups.values())
    if total <= 0:
        raise ValueError('At least one positive allocation is required.')
    for group in groups.values():
        group['allocation'] /= total
    return list(groups.values())


class MixedPressure:
    def __init__(self, oracle, baseline, mode, include_inactive, train_allocations, deadline):
        self.oracle, self.mode, self.deadline = oracle, mode, deadline
        self.q = oracle.alphabet
        self.groups = [grouped_places(baseline['row_places'], baseline['row_allocations'], include_inactive),
                       grouped_places(baseline['column_places'], baseline['column_allocations'], include_inactive)]
        row_labels = {str(label) for group in self.groups[0] for label in group['labels']}
        column_labels = {str(label) for group in self.groups[1] for label in group['labels']}
        assert row_labels.isdisjoint(column_labels)
        self.train_allocations = train_allocations
        self.queries = 0
        self.best = None
        self.height_size = 3 if mode == 'linear' else 4 * (self.q - 1)
        self.height_slices = []
        offset = 0
        for cut_groups in self.groups:
            slices = []
            for _ in cut_groups:
                slices.append(slice(offset, offset + self.height_size))
                offset += self.height_size
            self.height_slices.append(slices)
        self.type_slice = slice(offset, offset + 2 * (self.q - 1))
        offset += 2 * (self.q - 1)
        self.allocation_slices = []
        for groups in self.groups:
            size = len(groups) - 1 if train_allocations else 0
            self.allocation_slices.append(slice(offset, offset + size))
            offset += size
        self.scale_index = offset
        self.dimension = offset + 1
        self.active_heights = np.asarray([leg * self.q + a for leg in range(4)
                                          for a in range(1, self.q)], dtype=int)

    def height_vector(self, parameters):
        if self.mode == 'linear':
            return np.asarray(self.oracle.linear_heights((*parameters, 0.0)), dtype=float)
        result = np.zeros(4 * self.q)
        result[self.active_heights] = parameters
        return result

    def height_gradient(self, gradient):
        if self.mode == 'linear':
            return np.asarray([gradient[leg * self.q:(leg + 1) * self.q] @ np.arange(self.q)
                               for leg in range(3)])
        return gradient[self.active_heights]

    def initial(self, baseline, rng, start):
        values = np.zeros(self.dimension)
        scale = baseline['outer_log_scale']
        if start:
            scale *= (0.6, 1.0, 1.7, 2.5)[start % 4]
        values[self.scale_index] = np.clip(scale, 0.004, 3.5)
        for cut, groups in enumerate(self.groups):
            if self.train_allocations and len(groups) > 1:
                allocations = np.asarray([group['allocation'] for group in groups])
                logits = np.log(allocations[:-1] / allocations[-1])
                values[self.allocation_slices[cut]] = np.clip(logits, -6, 6)
            if start:
                for index, _ in enumerate(groups):
                    if self.mode == 'linear':
                        initial = rng.normal(0, 0.6 * scale, self.height_size)
                    else:
                        initial = rng.normal(0, 0.12 * scale, self.height_size)
                    values[self.height_slices[cut][index]] = initial
        if start:
            slopes = ((0.5, 0.5), (1.0, -1.0), (-0.7, -0.3), (1.5, 0.3), (0.0, 0.0))[start % 5]
            type_logits = np.concatenate([slope * np.arange(1, self.q) for slope in slopes])
            values[self.type_slice] = np.clip(type_logits + rng.normal(0, 0.05, len(type_logits)), -5, 5)
        return values

    def bounds(self):
        bounds = [(-12.0, 12.0)] * self.dimension
        for index in range(self.type_slice.start, self.type_slice.stop):
            bounds[index] = (-5.0, 5.0)
        for current in self.allocation_slices:
            for index in range(current.start, current.stop):
                bounds[index] = (-6.0, 6.0)
        bounds[self.scale_index] = (0.002, 4.0)
        return bounds

    def evaluate(self, values, smoothing=0.01, record=True):
        if time.monotonic() >= self.deadline:
            raise SearchDeadline()
        q = self.q
        scale = float(values[self.scale_index])
        type_values = values[self.type_slice].reshape(2, q - 1)
        types = [softmax(np.r_[0.0, row]) for row in type_values]
        entropies = [-float(probabilities @ np.log(probabilities)) for probabilities in types]
        entropy = sum(entropies)
        if entropy < 1e-8:
            raise ArithmeticError('Physical entropy is unresolved.')
        mean_vector = np.concatenate((types[0], types[1], types[1], types[0]))
        total = self.oracle.degree * scale
        exact_total = total
        derivative = np.zeros(self.dimension)
        derivative[self.scale_index] = self.oracle.degree
        type_derivative = [np.zeros(q), np.zeros(q)]
        place_records = []
        for cut, groups in enumerate(self.groups):
            if self.train_allocations and len(groups) > 1:
                allocations = softmax(np.r_[values[self.allocation_slices[cut]], 0.0])
            else:
                allocations = np.asarray([group['allocation'] for group in groups])
            costs, local_data = [], []
            for index, (group, allocation) in enumerate(zip(groups, allocations)):
                local_scale = scale * len(group['labels']) / allocation
                heights = self.height_vector(values[self.height_slices[cut][index]])
                log_sum, gradient_h, gradient_scale, orders = self.oracle.log_profile_with_gradient(
                    cut, group['point'], heights, local_scale)
                shift = float(mean_vector @ heights)
                costs.append(log_sum + shift)
                local_data.append((heights, gradient_h + mean_vector, gradient_scale,
                                   local_scale, orders, shift, log_sum))
            costs = np.asarray(costs)
            weights = softmax(costs / smoothing)
            total += smoothing * float(logsumexp(costs / smoothing))
            exact_total += float(np.max(costs))
            allocation_gradient = np.zeros(len(groups))
            cut_records = []
            for index, (group, allocation, weight, data) in enumerate(zip(groups, allocations, weights, local_data)):
                heights, grad_h, grad_scale, local_scale, orders, shift, log_sum = data
                derivative[self.height_slices[cut][index]] += weight * self.height_gradient(grad_h)
                derivative[self.scale_index] += weight * grad_scale * local_scale / scale
                allocation_gradient[index] = -weight * grad_scale * local_scale
                type_derivative[0] += weight * (heights[:q] + heights[3 * q:])
                type_derivative[1] += weight * (heights[q:2 * q] + heights[2 * q:3 * q])
                cut_records.append({'point': None if group['point'] is None else list(group['point']),
                                    'conjugate_labels': group['labels'], 'group_allocation': float(allocation),
                                    'allocation_per_root': float(allocation / len(group['labels'])),
                                    'local_log_scale': float(local_scale), 'heights': heights.tolist(),
                                    'restored_orders': orders.tolist(), 'type_shift': shift,
                                    'log_profile_sum': float(log_sum), 'complete_cost': float(costs[index])})
            if self.train_allocations and len(groups) > 1:
                transformed = allocation_gradient - allocations * np.sum(allocation_gradient)
                derivative[self.allocation_slices[cut]] += transformed[:-1]
            place_records.append(cut_records)
        rho = total / (2 * entropy)
        exact_rho = exact_total / (2 * entropy)
        derivative /= 2 * entropy
        type_result = []
        for probabilities, h, gradient in zip(types, entropies, type_derivative):
            cost_gradient = probabilities * (gradient - float(probabilities @ gradient)) / (2 * entropy)
            entropy_gradient = -probabilities * (np.log(probabilities) + h)
            type_result.extend((cost_gradient - rho * entropy_gradient / entropy)[1:])
        derivative[self.type_slice] = type_result
        self.queries += 1
        if record and (self.best is None or exact_rho < self.best['rho']):
            self.best = {'rho': float(exact_rho), 'smoothed_rho': float(rho),
                         'critical_pressure': float(exact_total - np.sqrt(2) * entropy),
                         'entropy_sum': entropy, 'paired_types': [row.tolist() for row in types],
                         'outer_log_scale': scale, 'mode': self.mode,
                         'parameters': values.tolist(), 'place_records': place_records,
                         'all_physical_costs_included': True}
        return float(rho), derivative


def atomic_write(path, result):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(result, indent=2) + '\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--seed', type=int, default=311602)
    parser.add_argument('--maxiter', type=int, default=160)
    parser.add_argument('--maxfun', type=int, default=1400)
    parser.add_argument('--trials', type=int, default=1000000)
    args = parser.parse_args()
    started = time.monotonic()
    deadline = started + args.seconds
    input_bytes = args.baseline.read_bytes()
    baselines = json.loads(input_bytes)['seeds']
    rng = np.random.default_rng(args.seed)
    oracles = {}
    records = []
    source_names = ('search_fused_mixed_weighted.py', 'fused_weighted_place_profiles.py',
                    'check_fused_closed_profiles.py', 'check_fused_multipoint_profiles.py',
                    'fused_six_vertex_crossings.py', 'f16_identity_operator_profiles.py', 'f4_lattice_linear.py')
    directory = Path(__file__).parent
    result = {'scope': 'Budgeted floating mixed-place full-cost search; no global optimum or factoring improvement certified.',
              'arguments': vars(args) | {'baseline': str(args.baseline), 'output': str(args.output)},
              'baseline_sha256': sha256(input_bytes).hexdigest(),
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in source_names},
              'records': records, 'budget_completed': False, 'trial_limit_completed': False}
    interrupted = None
    for trial in range(args.trials):
        if time.monotonic() >= deadline:
            break
        seed_entry = baselines[trial % len(baselines)]
        m, theta = seed_entry['seed']['m'], seed_entry['seed']['theta']
        if (m, theta) not in oracles:
            oracles[m, theta] = FusedPlaceOracle(m, theta)
        oracle = oracles[m, theta]
        round_index = trial // len(baselines)
        ranked = sorted(seed_entry['completed_partitions'], key=lambda row: row['rho'])
        baseline = ranked[(round_index // 8) % min(3, len(ranked))]
        mode = 'linear' if round_index % 2 == 0 else 'full'
        include_inactive = round_index % 6 == 5
        train_allocations = round_index % 4 != 0
        pressure = MixedPressure(oracle, baseline, mode, include_inactive, train_allocations, deadline)
        initial = pressure.initial(baseline, rng, round_index)
        trial_started = time.monotonic()
        outcomes = []
        complete = True
        try:
            pressure.evaluate(initial)
            for smoothing in (0.04, 0.008, 0.0015):
                optimized = minimize(lambda values: pressure.evaluate(values, smoothing), initial,
                                     jac=True, method='L-BFGS-B', bounds=pressure.bounds(),
                                     options={'maxiter': args.maxiter, 'maxfun': args.maxfun,
                                              'ftol': 1e-11, 'gtol': 2e-7, 'maxls': 35})
                outcomes.append({'smoothing': smoothing, 'optimizer_success': bool(optimized.success),
                                 'optimizer_message': str(optimized.message),
                                 'iterations': int(optimized.nit), 'function_evaluations': int(optimized.nfev)})
                initial = np.asarray(pressure.best['parameters'])
        except SearchDeadline:
            complete = False
            interrupted = 'The requested wall-clock budget expired during this trial.'
        except (ArithmeticError, ValueError) as error:
            complete = False
            interrupted = f'Numerical trial stopped: {error}'
        record = {'trial': trial, 'm': m, 'theta': theta, 'alphabet': oracle.alphabet,
                  'mode': mode, 'include_inactive': include_inactive,
                  'train_allocations': train_allocations, 'parameter_count': pressure.dimension,
                  'completed_trial': complete, 'outcomes': outcomes, 'queries': pressure.queries,
                  'elapsed_seconds': time.monotonic() - trial_started,
                  'best': pressure.best, 'interruption': None if complete else interrupted}
        records.append(record)
        best_records = [row for row in records if row['best'] is not None]
        result.update(elapsed_seconds=time.monotonic() - started,
                      completed_trials=sum(row['completed_trial'] for row in records),
                      best=min(best_records, key=lambda row: row['best']['rho']) if best_records else None)
        atomic_write(args.output, result)
        print(json.dumps({'trial': trial, 'q': oracle.alphabet, 'mode': mode,
                          'rho': None if pressure.best is None else pressure.best['rho'],
                          'queries': pressure.queries, 'complete': complete,
                          'seconds': time.monotonic() - started}), flush=True)
        if time.monotonic() >= deadline:
            break
    result.update(elapsed_seconds=time.monotonic() - started,
                  budget_completed=time.monotonic() >= deadline,
                  trial_limit_completed=len(records) == args.trials and all(row['completed_trial'] for row in records))
    atomic_write(args.output, result)


if __name__ == '__main__':
    main()
