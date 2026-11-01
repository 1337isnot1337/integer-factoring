#!/usr/bin/env python3
"""Exact finite Bell jets while retaining regular reflection frames.

The default lane fixes every original column-infinity row target and
relaxes the row-zero targets inside one explicitly bounded regular frame.
An explicit infinity-window parameter permits a separately charged trade.
Candidates are complete polynomial tensors, with no paired-trace or
derivative constraint imposed. Numerical pressure is exploratory.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path
import random
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from affine_place_profile_oracle import AffinePlaceOracle, primitive_seed
from entrywise_reflection_constraints import entries_for_rectangle
from finite_jet_smith_fibers import FiniteJetFiber
from finite_precision_smith_frames import build_lattice
from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import digest
from search_entrywise_reflection_jets import basis_catalogue, bell_values
from search_full_alphabet_two_cut_lattices import root_for_profiles
from search_fused_mixed_weighted import atomic_write


def mixed_fixed_gauge_score(row_zero, finite_row, column_infinity, degree=24):
    row_zero, finite_row, column_infinity = map(lambda x: np.asarray(x, dtype=float),
                                              (row_zero, finite_row, column_infinity))
    energy = np.asarray([0., 6., 2., 8.])
    best = None

    def branch(orders, scale, mean=0.):
        exponents = -scale * orders
        normalizer = float(logsumexp(exponents))
        probabilities = np.exp(exponents - normalizer)
        return normalizer + scale * mean, mean - float(probabilities @ orders)

    def evaluate(values, smoothing=0.):
        nonlocal best
        t, beta = math.exp(float(values[0])), float(values[1])
        alpha = 1 / (1 + math.exp(-float(values[2])))
        logits = -beta * energy
        logz = float(logsumexp(logits))
        probabilities = np.exp(logits - logz)
        mean = float(probabilities @ energy)
        variance = float(probabilities @ ((energy - mean) ** 2))
        typed_mean, entropy = 2 * mean, 2 * (logz + beta * mean)
        if entropy <= 1e-10:
            return 1e10, np.zeros(3)
        first, dfirst = branch(row_zero, t / alpha, typed_mean)
        second, dsecond = branch(finite_row, t / (1 - alpha))
        column, dcolumn = branch(column_infinity, t, typed_mean)
        actual = (degree * t + max(first, second) + column) / (2 * entropy)
        if best is None or actual < best['rho']:
            best = {'rho': actual, 'outer_log_scale': t, 'row_allocations': [alpha, 1 - alpha],
                    'column_allocations': [1.], 'paired_type': probabilities.tolist(),
                    'beta': beta, 'entropy_sum': entropy, 'parameters': list(map(float, values))}
        if smoothing:
            logweights = np.asarray([first, second]) / smoothing
            denominator = float(logsumexp(logweights))
            weights = np.exp(logweights - denominator)
            row = smoothing * denominator
        else:
            weights = np.asarray([float(first >= second), float(first < second)])
            row = max(first, second)
        numerator = degree * t + row + column
        dmean, dentropy = -2 * variance, -2 * beta * variance
        derivative = np.asarray([
            degree * t + weights[0] * (t / alpha) * dfirst
            + weights[1] * (t / (1 - alpha)) * dsecond + t * dcolumn,
            (weights[0] * t / alpha + t) * dmean - numerator * dentropy / entropy,
            -weights[0] * (t / alpha) * (1 - alpha) * dfirst
            + weights[1] * (t / (1 - alpha)) * alpha * dsecond,
        ]) / (2 * entropy)
        return numerator / (2 * entropy), derivative

    for allocation in (.95, .65, .3):
        values = np.asarray([math.log(.329011123), .450009, math.log(allocation / (1 - allocation))])
        for smoothing in (.02, .001):
            fit = minimize(lambda v: evaluate(v, smoothing), values, jac=True, method='L-BFGS-B',
                           bounds=[(-8, 2), (-3, 3), (-9, 9)],
                           options={'maxiter': 80, 'ftol': 3e-11, 'gtol': 2e-7, 'maxls': 30})
            values = np.asarray(best['parameters'])
    return best


def templates(matrices):
    for order in (1, 2, 4, 6, 8, 12):
        for matrix_name, matrix in matrices.items():
            target = bell_values(matrix, 0)
            for scalar in (3, 1, 2):
                requests = [(entry, 2, jet, F4_PRODUCT[scalar][target[entry]] if jet == 0 else 0)
                            for jet in range(order) for entry in range(256)]
                yield {'kind': 'row Bell jets', 'point': 2, 'order': order,
                       'matrix': matrix_name, 'scalar': scalar}, requests
    halves = list(combinations(range(4), 2))
    for order in (2, 4, 6):
        for first in halves:
            for second in halves:
                entries = entries_for_rectangle((2, 0, first, second, order))
                requests = [(entry, 2, jet, 0) for jet in range(order) for entry in entries]
                yield {'kind': 'row rectangle jets', 'point': 2, 'order': order,
                       'selected_rows': first, 'selected_columns': second}, requests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1200)
    parser.add_argument('--lower-window', type=int, default=6)
    parser.add_argument('--infinity-window', type=int, default=0)
    parser.add_argument('--upper-window', type=int, default=8)
    parser.add_argument('--samples', type=int, default=2)
    parser.add_argument('--seed', type=int, default=312523)
    args = parser.parse_args()
    started, rng = time.monotonic(), random.Random(args.seed)
    deadline = started + args.seconds
    base = tensor(4)
    h, zero = [0, 6, 2, 8], [0] * 4
    lattice = build_lattice(base['values'], 4, 24, [zero, zero, h, h, h, h, zero, zero], args.upper_window)
    floors = tuple(max(0, value - args.lower_window) for value in lattice.orders)
    fiber = FiniteJetFiber(lattice, floors)
    fixed_targets = tuple(max(floors[i], lattice.orders[i] - args.infinity_window) for i in range(16)) + tuple(floors[16:])
    initial_pivots = fiber.solve([], fixed_targets)
    assert initial_pivots is not None
    matrices = basis_catalogue(base)
    matrices['alternating'] = [[int(j == (i ^ 1)) for j in range(4)] for i in range(4)]
    matrices['upper'] = [[int(i == j) ^ (2 if j == i + 1 else 0) for j in range(4)] for i in range(4)]
    names = (Path(__file__).name, 'finite_jet_smith_fibers.py', 'finite_precision_smith_frames.py',
             'full_alphabet_two_cut_lattices.py', 'entrywise_reflection_constraints.py',
             'search_full_alphabet_two_cut_lattices.py', 'affine_place_profile_oracle.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')
    report = {'scope': __doc__, 'arguments': vars(args) | {'output': str(args.output)},
              'source_tensor_sha256': digest(base['values']), 'physical_heights': h,
              'frame_reference_orders': lattice.orders, 'common_floors': floors,
              'fixed_initial_targets': fixed_targets, 'fiber_variables': fiber.variables,
              'fixed_frame_solution_dimension': fiber.variables - len(initial_pivots),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'records': [], 'retained_seeds': [], 'completed': False, 'complete_template_list': False}
    seen = set()
    for index, (metadata, requests) in enumerate(templates(matrices)):
        if time.monotonic() >= deadline:
            break
        began = time.monotonic()
        equations = [fiber.jet_equation(*request) for request in requests]
        pivots = fiber.solve(equations, fixed_targets)
        record = metadata | {'target_index': index, 'feasible': pivots is not None}
        if pivots is not None:
            record['initial_solution_dimension'] = fiber.variables - len(pivots)
            targets = list(fixed_targets)
            order = list(range(0 if args.infinity_window else 16, 32))
            rng.shuffle(order)
            for row in order:
                while targets[row] < lattice.caps[row]:
                    candidate = fiber.increment(pivots, row, targets[row])
                    if candidate is None:
                        break
                    pivots = candidate
                    targets[row] += 1
            record.update(final_solution_dimension=fiber.variables - len(pivots), row_targets=targets)
            samples = []
            for sample in range(args.samples):
                values = fiber.sample(pivots, None if sample == 0 else rng)
                fiber.verify_jets(values, requests)
                key = digest(values)
                if key in seen:
                    continue
                seen.add(key)
                profiles = lattice.check_witness(values, targets)
                old_cost = root_for_profiles(profiles, [2 * value for value in h], [2 * value for value in h])
                seed = {'q_x': 4, 'q_u': 4, 'dimensions': [4] * 4, 'values': values, 'degree': 24,
                        'strict_endpoint': True, 'denominator_degree': 0,
                        'construction': 'Finite Hasse equations inside a complete-entry regular reflection fiber'}
                oracle = AffinePlaceOracle(seed)
                finite = list(map(round, oracle.log_profile_with_gradient(0, (2, 1), [0] * 16)[3]))
                score = mixed_fixed_gauge_score(profiles[1], finite, [value + 24 for value in profiles[0]])
                primitive, common = primitive_seed(seed)
                item = {'target_index': index, 'sample': sample, 'tensor_sha256': key,
                        'degree_cap': 24, 'actual_degree': primitive['degree'], 'common_factor': common,
                        'two_place_profiles': profiles, 'finite_row_orders': finite,
                        'old_natural_cost': old_cost, 'mixed_fixed_gauge': score, 'seed': primitive,
                        'metadata': metadata, 'row_targets': targets}
                samples.append({key: value for key, value in item.items() if key != 'seed'})
                report['retained_seeds'].append(item)
                retained = report['retained_seeds']
                best_cost = sorted(retained, key=lambda row: row['mixed_fixed_gauge']['rho'])[:80]
                low_degree = sorted(retained, key=lambda row: (row['actual_degree'], row['mixed_fixed_gauge']['rho']))[:20]
                report['retained_seeds'] = list({row['tensor_sha256']: row for row in best_cost + low_degree}.values())
            record['samples'] = samples
        record['seconds'] = time.monotonic() - began
        report['records'].append(record)
        report.update(elapsed_seconds=time.monotonic() - started, completed_targets=len(report['records']),
                      feasible_targets=sum(row['feasible'] for row in report['records']), distinct_seeds=len(seen))
        atomic_write(args.output, report)
        print(json.dumps({'target': index, 'kind': metadata['kind'], 'order': metadata['order'],
                          'feasible': record['feasible'], 'best_mixed_rho': min((row['mixed_fixed_gauge']['rho'] for row in report['retained_seeds']), default=None),
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    else:
        report['complete_template_list'] = True
    report.update(completed=True, elapsed_seconds=time.monotonic() - started,
                  stopping_reason='complete finite template list' if report['complete_template_list'] else 'declared wall-clock cap')
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
