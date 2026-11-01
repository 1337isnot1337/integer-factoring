#!/usr/bin/env python3
"""Build exact local gauge catalogues from bounded typed jet optimization.

The physical basis and uniform paired types are fixed throughout. Rational
gauges and complete profiles are exact; jet-density pressure and all search
decisions remain numerical until a separate all-length certificate exists.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp, softmax

from search_ternary_physical_bases import change_basis, oracles_for


BASIS = [[1, 0, 0], [0, 2, 1], [2, 2, 2]]
SLOPES = ([0, 0, 0, 0, .5, 1, .5, 1],
          [-.5, -1, -.5, -1, 0, 0, 0, 0],
          [-1, -2, -1, -2, 1, 2, 1, 2])
UNIFORM = np.full(3, 1 / 3)


def local_objective(parameters, oracle, density, gradient=True):
    rate, heights = float(parameters[0]), np.asarray(parameters[1:])
    value, derivative, _, _ = oracle.cost(rate, heights, UNIFORM, UNIFORM, True)
    value += rate * density
    derivative[0] += density
    return (float(value), derivative) if gradient else float(value)


def exact_integral_profile(oracle, ramification, heights):
    """Complete weighted Schur profile using only integral minor valuations."""
    rows, columns = list(range(9)), list(range(9))
    row_heights = [sum(int(a) * b for a, b in zip(row, heights)) for row in oracle.row_counts]
    column_heights = [sum(int(a) * b for a, b in zip(row, heights)) for row in oracle.column_counts]
    rm = cm = total_height = previous = 0
    orders = []
    for _ in range(9):
        candidates = []
        for r in rows:
            for c in columns:
                value = int(oracle.valuations[rm | (1 << r), cm | (1 << c)])
                if value < 10000:
                    candidates.append((ramification * value + row_heights[r] + column_heights[c], r, c))
        _, r, c = min(candidates)
        rows.remove(r)
        columns.remove(c)
        rm |= 1 << r
        cm |= 1 << c
        total_height += row_heights[r] + column_heights[c]
        current = ramification * int(oracle.valuations[rm, cm]) + total_height
        orders.append(current - previous)
        previous = current
    assert orders == sorted(orders)
    assert orders == list(map(int, oracle.profiles(ramification, np.array(heights))))
    return orders


def make_gauge(oracle, ramification, original_heights):
    heights = [round(float(value) * ramification) for value in original_heights]
    orders = exact_integral_profile(oracle, ramification, heights)
    mean = Fraction(sum(heights), 3)
    minimum = min(orders)
    centered = Fraction(mean - minimum)
    key = (tuple(order - minimum for order in orders), centered.numerator, centered.denominator)
    return key, {'height_denominator': ramification, 'integral_heights': heights,
                 'orders': orders, 'mean_shift': [mean.numerator, mean.denominator]}


def effective_orders(gauge):
    mean = gauge['mean_shift'][0] / gauge['mean_shift'][1]
    return (np.array(gauge['orders']) - mean) / gauge['height_denominator']


def curve(gauge, density):
    orders = effective_orders(gauge)
    if density >= float(np.mean(orders)) - 1e-13:
        return math.log(9), 0.0
    minimum = float(min(orders))
    if density < minimum - 1e-13:
        return -math.inf, math.inf
    if abs(density - minimum) < 1e-13:
        return math.log(int(np.count_nonzero(np.abs(orders - minimum) < 1e-13))), math.inf
    def function(scale):
        return float(logsumexp(-scale * orders) + scale * density)
    high = 1.0
    while float(softmax(-high * orders) @ orders) > density and high < 1e6:
        high *= 2
    result = minimize_scalar(function, bounds=(0, high), method='bounded', options={'xatol': 1e-11})
    return min(math.log(9), float(result.fun)), float(result.x)


def save(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=900)
    parser.add_argument('--grid', type=int, default=20)
    parser.add_argument('--restarts', type=int, default=8)
    parser.add_argument('--iterations', type=int, default=400)
    parser.add_argument('--denominator', type=int, default=16)
    parser.add_argument('--seed', type=int, default=81131)
    parser.add_argument('--output', type=Path, default=Path('/tmp/ternary-gauge-catalogues.json'))
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    oracles, minor_count = oracles_for(change_basis(np.array(BASIS), np.array(BASIS)))
    catalogues, keys = [[], [], []], [set(), set(), set()]
    for place, oracle in enumerate(oracles):
        for slopes in (SLOPES[place], [0] * 8):
            key, gauge = make_gauge(oracle, args.denominator, slopes)
            keys[place].add(key)
            catalogues[place].append(gauge)
    # Derivatives are checked away from selected-minor and row-cost ties.
    largest_error = 0.0
    for oracle in oracles:
        x = np.r_[.73, rng.normal(0, .2, 8)]
        value, gradient = local_objective(x, oracle, .37)
        numerical = np.array([(local_objective(x + np.eye(9)[i] * 1e-6, oracle, .37, False)
                               - local_objective(x - np.eye(9)[i] * 1e-6, oracle, .37, False)) / 2e-6
                              for i in range(9)])
        largest_error = max(largest_error, float(np.max(np.abs(gradient - numerical))))
    assert largest_error < 2e-5
    files = (Path(__file__).name, 'search_ternary_physical_bases.py', 'ternary_mixed_pressure.py')
    report = {'scope': 'fixed-basis uniform-type rational gauge catalogues; numerical pressure search only',
              'basis': BASIS, 'paired_type_x': [[1, 3]] * 3, 'paired_type_y': [[1, 3]] * 3,
              'places': [2, 'infinity', 0], 'cuts': [0, 0, 1], 'seed': args.seed,
              'requested_seconds': args.seconds, 'height_denominator': args.denominator,
              'initial_minor_valuations': minor_count, 'gradient_controls': 3,
              'maximum_gradient_error': largest_error,
              'source_sha256': {name: sha256((Path(__file__).parent / name).read_bytes()).hexdigest() for name in files},
              'catalogues': catalogues, 'records': []}
    begun, batch = time.monotonic(), 0
    starts = [[np.r_[1.0, np.array(slopes)]] for slopes in SLOPES]
    maximum_improvement = 0.0
    while time.monotonic() - begun < args.seconds:
        if batch < 3 * (args.grid - 1):
            place = batch % 3
            fraction = (batch // 3 + 1) / args.grid
        else:
            place, fraction = int(rng.integers(0, 3)), float(rng.uniform(.005, .995))
        density = fraction * (2 if place == 2 else 1)
        oracle = oracles[place]
        reference, scale = curve(catalogues[place][0], density)
        clean_start = np.r_[scale, scale * np.array(SLOPES[place])]
        retained = [reference, clean_start.copy()]
        record = {'place_index': place, 'density': density, 'clean_pressure': reference, 'runs': []}
        for restart in range(args.restarts):
            if time.monotonic() - begun >= args.seconds:
                break
            if restart == 0:
                start = clean_start
            elif restart % 3 == 1:
                start = retained[1] + rng.normal(0, .2, 9)
            elif restart % 3 == 2:
                start = np.r_[float(rng.uniform(.1, 5)), rng.normal(0, 1.5, 8)]
            else:
                start = starts[place][int(rng.integers(0, len(starts[place])))] + rng.normal(0, .1, 9)
            start[0] = min(64, max(1e-4, start[0]))
            def evaluate(parameters):
                value, derivative = local_objective(parameters, oracle, density)
                if value < retained[0]:
                    retained[:] = [value, np.asarray(parameters).copy()]
                return value, derivative
            result = minimize(evaluate, start, jac=True, method='SLSQP',
                              bounds=[(1e-5, 64)] + [(-128, 128)] * 8,
                              options={'maxiter': args.iterations, 'ftol': 2e-12})
            record['runs'].append({'restart': restart, 'value': float(result.fun),
                                   'iterations': int(result.nit), 'evaluations': int(result.nfev),
                                   'success': bool(result.success), 'message': str(result.message)})
            for parameters in (result.x, retained[1]):
                key, gauge = make_gauge(oracle, args.denominator, parameters[1:] / parameters[0])
                if key not in keys[place]:
                    keys[place].add(key)
                    catalogues[place].append(gauge)
            record.update(best_pressure=float(retained[0]), best_parameters=retained[1].tolist())
            maximum_improvement = max(maximum_improvement, reference - retained[0])
            report.update(maximum_continuous_improvement=maximum_improvement,
                          elapsed_seconds=time.monotonic() - begun)
            report['records'] = report['records'][:batch] + [record]
            save(args.output, report)
        starts[place].append(retained[1].copy())
        print(json.dumps({'batches': batch + 1, 'catalogue_sizes': list(map(len, catalogues)),
                          'maximum_continuous_improvement': maximum_improvement,
                          'elapsed_seconds': time.monotonic() - begun}), flush=True)
        batch += 1
    report.update(complete=True, elapsed_seconds=time.monotonic() - begun)
    save(args.output, report)


if __name__ == '__main__':
    main()
