#!/usr/bin/env python3
"""Search over endpoint-preserving physical bases for ternary seeds.

Numerical optimization only; no asymptotic certificate is asserted. Exact
polynomial minors supply every valuation. Progress is written atomically.
"""
from __future__ import annotations

import argparse
from itertools import combinations, product
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, softmax

import ternary_mixed_pressure as original
from check_multipoint_hermite_crossings import mixing_seed, local_profiles
from check_weighted_multipoint_hermite import local_gauged_profile


def inverse(matrix):
    size = len(matrix)
    work = np.concatenate((np.asarray(matrix, dtype=np.int64) % 3, np.eye(size, dtype=np.int64)), axis=1)
    for j in range(size):
        found = next((i for i in range(j, size) if work[i, j]), None)
        if found is None:
            raise ValueError('singular basis')
        work[[j, found]] = work[[found, j]]
        work[j] = work[j] * pow(int(work[j, j]), -1, 3) % 3
        for i in range(size):
            if i != j:
                work[i] = (work[i] - work[i, j] * work[j]) % 3
    return work[:, size:]


def change_basis(a, b):
    seed = mixing_seed()
    dense = np.zeros((3, 3, 3, 3, 3), dtype=np.int64)
    for key, value in seed.items():
        dense[key][:len(value)] = value
    c, d = inverse(b).T, inverse(a).T
    transformed = np.einsum('ax,by,cu,dv,xyuvk->abcdk', a, b, c, d, dense, optimize=True) % 3
    result = {key: transformed[key].tolist() for key in product(range(3), repeat=4)
              if np.any(transformed[key])}
    for x, y, u, v in product(range(3), repeat=4):
        assert int(sum(result.get((x, y, u, v), [])) % 3) == int(x == v and y == u)
    return result


def all_minor_valuations(seed, cut, points):
    """All minors by vectorized exact Laplace recursion over F3[z]."""
    n = 9
    labels = list(product(range(3), repeat=2))
    matrix = np.zeros((n, n, 3), dtype=np.int64)
    for i, (a, b) in enumerate(labels):
        for j, (c, d) in enumerate(labels):
            key = (a, b, c, d) if cut == 0 else (a, c, b, d)
            value = seed.get(key, [])
            matrix[i, j, :len(value)] = value
    result = {point: np.full((1 << n, 1 << n), 10000, dtype=np.int16) for point in points}
    for values in result.values():
        values[0, 0] = 0
    previous = {0: np.ones((1, 1), dtype=np.int64)}
    previous_columns = {0: 0}
    minor_count = 0
    for k in range(1, n + 1):
        subsets = list(combinations(range(n), k))
        subset_array = np.array(subsets, dtype=np.int64)
        masks = np.array([sum(1 << i for i in subset) for subset in subsets], dtype=np.int64)
        removed = np.array([[previous_columns[int(mask) ^ (1 << int(j))] for j in subset]
                            for mask, subset in zip(masks, subsets)], dtype=np.int64)
        translations = {point: np.array([[math.comb(i, j) * pow(point, i - j, 3) % 3 if i >= j else 0
                                          for j in range(2 * k + 1)] for i in range(2 * k + 1)], dtype=np.int64)
                        for point in points if point is not None and point != 0}
        current = {}
        for subset, mask in zip(subsets, masks):
            last = subset[-1]
            prior = previous[int(mask) ^ (1 << last)]
            polynomials = np.zeros((len(subsets), 2 * k + 1), dtype=np.int64)
            for position in range(k):
                subminors = prior[removed[:, position]]
                entries = matrix[last, subset_array[:, position]]
                sign = -1 if (k - 1 + position) % 2 else 1
                for exponent in range(3):
                    polynomials[:, exponent:exponent + 2 * k - 1] += (
                        sign * entries[:, exponent, None] * subminors)
            polynomials %= 3
            current[int(mask)] = polynomials
            nonzero = np.any(polynomials, axis=1)
            minor_count += int(np.count_nonzero(nonzero))
            for point in points:
                local = polynomials[:, ::-1] if point is None else (
                    polynomials if point == 0 else polynomials @ translations[point] % 3)
                orders = np.argmax(local != 0, axis=1)
                result[point][int(mask), masks] = np.where(nonzero, orders, 10000)
        previous = current
        previous_columns = {int(mask): i for i, mask in enumerate(masks)}
    return result, minor_count


class DenseOracle:
    def __init__(self, valuations, point, cut):
        self.valuations, self.point, self.cut = valuations, point, cut
        labels = list(product(range(3), repeat=2))
        row_legs, column_legs = ((0, 1), (2, 3)) if cut == 0 else ((0, 2), (1, 3))
        self.row_counts = np.zeros((9, 8), dtype=np.int64)
        self.column_counts = np.zeros((9, 8), dtype=np.int64)
        for index, letters in enumerate(labels):
            for leg, symbol in zip(row_legs, letters):
                if symbol:
                    self.row_counts[index, 2 * leg + symbol - 1] += 1
            for leg, symbol in zip(column_legs, letters):
                if symbol:
                    self.column_counts[index, 2 * leg + symbol - 1] += 1

    def profiles(self, rate, weights, derivatives=False):
        row_heights, column_heights = self.row_counts @ weights, self.column_counts @ weights
        rows, columns = list(range(9)), list(range(9))
        row_mask = column_mask = 0
        previous, previous_gradient = 0.0, np.zeros(9)
        counts = np.zeros(8)
        orders, gradients = [], []
        for _ in range(9):
            rm = np.array([row_mask | (1 << i) for i in rows])
            cm = np.array([column_mask | (1 << j) for j in columns])
            valuation = self.valuations[rm[:, None], cm[None, :]].astype(float)
            candidates = rate * valuation + row_heights[rows, None] + column_heights[None, columns]
            candidates[valuation >= 10000] = np.inf
            flat = int(np.argmin(candidates))
            i, j = np.unravel_index(flat, candidates.shape)
            r, c = rows.pop(i), columns.pop(j)
            row_mask |= 1 << r
            column_mask |= 1 << c
            counts += self.row_counts[r] + self.column_counts[c]
            base = int(self.valuations[row_mask, column_mask])
            assert base < 10000
            current_gradient = np.concatenate(([base], counts.copy()))
            current = base * rate + counts @ weights
            orders.append(current - previous)
            gradients.append(current_gradient - previous_gradient)
            previous, previous_gradient = current, current_gradient
        values = np.array(orders)
        if derivatives:
            return values, np.array(gradients)
        return values

    def cost(self, rate, weights, px, py, derivatives=False):
        orders, gradients = self.profiles(rate, weights, True)
        means = np.concatenate((px[1:], py[1:], py[1:], px[1:]))
        value = logsumexp(-orders) + means @ weights
        if not derivatives:
            return value
        gradient = -softmax(-orders) @ gradients
        gradient[1:] += means
        physical = weights.reshape((4, 2))
        return (value, gradient, np.concatenate(([0], physical[0] + physical[3])),
                np.concatenate(([0], physical[1] + physical[2])))


def oracles_for(seed):
    row, row_count = all_minor_valuations(seed, 0, (2, None))
    column, column_count = all_minor_valuations(seed, 1, (0,))
    return ([DenseOracle(row[2], 2, 0), DenseOracle(row[None], None, 0),
             DenseOracle(column[0], 0, 1)], row_count + column_count)


def controls(seed, oracles, rng, count=3):
    checked = 0
    for oracle in oracles:
        assert sorted(oracle.profiles(1, np.zeros(8))) == local_profiles(seed, 3, 2, oracle.point, 3)[oracle.cut]
        for _ in range(count):
            weights = rng.integers(-2, 3, size=8)
            expected = local_gauged_profile(seed, 3, 2, oracle.point, 3, oracle.cut,
                                           [(0, int(weights[2 * i]), int(weights[2 * i + 1])) for i in range(4)])
            actual = sorted(oracle.profiles(1, weights))
            assert actual == expected, (oracle.point, oracle.cut, weights.tolist(), actual, expected)
            checked += 1
    for _ in range(count):
        vector = rng.normal(0, .3, size=30)
        vector[0] -= .7
        value, gradient = original.objective(vector, True)
        numerical = np.array([
            (original.objective(vector + np.eye(30)[i] * 1e-6)
             - original.objective(vector - np.eye(30)[i] * 1e-6)) / 2e-6
            for i in range(30)])
        assert np.max(np.abs(gradient - numerical)) < 2e-5
    return checked


def write_atomic(path, record):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(record, indent=2) + '\n')
    temporary.replace(path)


def random_basis(rng, index):
    eye = np.eye(3, dtype=np.int64)
    pascal = np.array([[1, 0, 0], [1, 1, 0], [1, 2, 1]], dtype=np.int64)
    evaluations = np.array([[1, 0, 0], [1, 1, 1], [1, 2, 1]], dtype=np.int64)
    special = [(eye, eye), (pascal, pascal), (inverse(pascal), inverse(pascal)),
               (evaluations, evaluations), (inverse(evaluations), inverse(evaluations)),
               (eye, pascal), (pascal, eye), (eye, inverse(pascal)),
               (inverse(pascal), eye), (pascal, inverse(pascal))]
    if index < len(special):
        return special[index]
    output = []
    for _ in range(2):
        while True:
            matrix = rng.integers(0, 3, size=(3, 3), dtype=np.int64)
            try:
                inverse(matrix)
            except ValueError:
                continue
            output.append(matrix)
            break
    if index % 3 == 0:
        output[1] = output[0].copy()
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--seed', type=int, default=21831)
    parser.add_argument('--restarts', type=int, default=6)
    parser.add_argument('--controls-only', action='store_true')
    parser.add_argument('--warm', action='store_true',
                        help='start from the exact uniform-type polynomial-basis family')
    parser.add_argument('--output', type=Path, default=Path('/tmp/ternary-basis-progress.json'))
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    base = np.zeros(30)
    base[0] = math.log(.51)
    clean = np.zeros(30)
    temperature = math.log((1 + math.sqrt(13)) / 2)
    clean[0] = math.log(temperature)
    clean[6:14] = temperature * np.array([0, 0, 0, 0, 1, 2, 1, 2])
    clean[14:22] = temperature * np.array([-1, -2, -1, -2, 0, 0, 0, 0])
    clean[22:30] = temperature * np.array([-1, -2, -1, -2, 1, 2, 1, 2])
    clean_basis = np.array([[1, 0, 0], [0, 2, 1], [2, 2, 2]], dtype=np.int64)
    started, record = time.monotonic(), {
        'schema_version': 1, 'seed': args.seed, 'requested_seconds': args.seconds,
        'scope': 'bounded floating basis/gauge/type search; no certificate or global optimality claim',
        'baseline_rho': original.objective(base),
        'best_rho': original.objective(base), 'best_vector': base.tolist(), 'records': []}
    if args.warm:
        original.ORACLES, _ = oracles_for(change_basis(clean_basis, clean_basis))
        record.update(best_rho=original.objective(clean), best_vector=clean.tolist(),
                      clean_family_rho=original.objective(clean),
                      best_basis_a=clean_basis.tolist(), best_basis_b=clean_basis.tolist(),
                      best_basis_index='exact uniform-type family')
    index = 0
    while time.monotonic() - started < args.seconds:
        a, b = random_basis(rng, index)
        if args.warm and index % 5 == 0:
            a, b = clean_basis.copy(), clean_basis.copy()
        seed = change_basis(a, b)
        preprocessing = time.monotonic()
        oracles, minor_count = oracles_for(seed)
        original.ORACLES = oracles
        profile_controls = controls(seed, oracles, rng, 3 if index < 5 else 1)
        basis_record = {'index': index, 'basis_a': a.tolist(), 'basis_b': b.tolist(),
                        'exact_nonzero_minor_valuations': minor_count,
                        'exact_weighted_profile_controls': profile_controls,
                        'finite_difference_gradient_controls': 3 if index < 5 else 1,
                        'preprocessing_seconds': time.monotonic() - preprocessing, 'runs': []}
        assert abs(original.objective(base) - record['baseline_rho']) < 1e-9
        if args.controls_only:
            record['records'].append(basis_record)
            record['elapsed_seconds'] = time.monotonic() - started
            write_atomic(args.output, record)
            print(json.dumps(basis_record), flush=True)
            if index >= 4:
                break
            index += 1
            continue
        local_best, local_vector = original.objective(base), base.copy()
        if args.warm and original.objective(clean) < local_best:
            local_best, local_vector = original.objective(clean), clean.copy()
        for run in range(args.restarts):
            if time.monotonic() - started >= args.seconds:
                break
            if run == 0:
                start = local_vector + rng.normal(0, .015, 30)
            elif run % 3 == 1:
                start = local_vector + rng.normal(0, .2, 30)
            elif run % 3 == 2:
                start = base + rng.normal(0, .7, 30)
            else:
                start = np.array(record['best_vector']) + rng.normal(0, .1, 30)
            result = minimize(lambda vector: original.objective(vector, True), start, jac=True,
                              method='L-BFGS-B',
                              bounds=[(-4, 1.5), (-3.5, 3.5)] + [(-5, 5)] * 4 + [(-8, 8)] * 24,
                              options={'maxiter': 700, 'maxls': 40, 'ftol': 1e-11, 'gtol': 1e-7})
            basis_record['runs'].append({'run': run, 'rho': float(result.fun),
                                         'iterations': int(result.nit), 'success': bool(result.success)})
            if result.fun < local_best:
                local_best, local_vector = float(result.fun), result.x.copy()
            if result.fun < record['best_rho']:
                record.update(best_rho=float(result.fun), best_vector=result.x.tolist(),
                              best_basis_a=a.tolist(), best_basis_b=b.tolist(), best_basis_index=index)
                print(json.dumps({'basis': index, 'run': run, 'best_rho': record['best_rho']}), flush=True)
            basis_record['best_rho'] = local_best
            basis_record['best_vector'] = local_vector.tolist()
            record['records'] = record['records'][:index] + [basis_record]
            record['elapsed_seconds'] = time.monotonic() - started
            write_atomic(args.output, record)
        print(json.dumps({'basis': index, 'best_rho': local_best,
                          'preprocessing_seconds': basis_record['preprocessing_seconds'],
                          'elapsed_seconds': time.monotonic() - started}), flush=True)
        index += 1
    record['elapsed_seconds'] = time.monotonic() - started
    record['complete'] = True
    write_atomic(args.output, record)
    print(json.dumps({key: record[key] for key in ('best_rho', 'elapsed_seconds', 'complete')}), flush=True)


if __name__ == '__main__':
    main()
