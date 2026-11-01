#!/usr/bin/env python3
"""Endpoint-only affine two-letter frames with a complete Gaussian hot tail.

L(t)=I+(1+t)A and P(t)=I+(1+t)B for binary two-by-two A,B.
The full low-mode entropy is retained. The 8-mode active matrix has
order256 and full generic rank256. All coefficients and both Laurent
profiles are computed completely; numerical pressure fits are separate.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

from native_matrix_smith_f4 import MatrixOracle
from tensor_literal_full_tail_core import images, left_action, right_action, number_projector, pmul


def affine(mask):
    if not 0 <= mask < 16:
        raise ValueError('Expected a four-bit matrix mask.')
    return [[int(i == j) ^ (3 if mask >> (2*i+j) & 1 else 0) for j in range(2)] for i in range(2)]


def endomorphism(matrix, row_heights=None, column_heights=None, reverse_parameter=False, transpose=False):
    if transpose:
        matrix = [list(row) for row in zip(*matrix)]
    row_heights = [0, 0] if row_heights is None else row_heights
    column_heights = [0, 0] if column_heights is None else column_heights
    terms = [[[] for _ in range(2)] for _ in range(2)]
    for i in range(2):
        for j in range(2):
            value = matrix[i][j]
            for exponent in range(value.bit_length()):
                if value >> exponent & 1:
                    terms[i][j].append((-exponent if reverse_parameter else exponent)
                                       +row_heights[i]+column_heights[j])
    clearing = max(0, -min(e for row in terms for cell in row for e in cell))
    polynomial = [[sum(1 << (e+clearing) for e in cell) for cell in row] for row in terms]
    return polynomial, clearing


def local_action(matrix, operators, polynomial, right):
    a, b = polynomial[0]
    c, d = polynomial[1]
    coefficients = [('identity', a), ('number', a ^ d), ('annihilator', b), ('creator', c)]
    size = len(matrix)
    out = [[0]*size for _ in range(size)]
    for kind, coefficient in coefficients:
        if not coefficient:
            continue
        if kind == 'identity':
            changed = matrix
        elif kind == 'number':
            changed = number_projector(matrix, operators, 1, right)
        else:
            operator = operators[0 if kind == 'creator' else 1]
            changed = (right_action if right else left_action)(matrix, operator)
        for i, row in enumerate(changed):
            for j, value in enumerate(row):
                if value:
                    out[i][j] ^= pmul(coefficient, value)
    return out


def core_matrix(M, a_mask, b_mask, cut, low_height=1):
    L, P = affine(a_mask), affine(b_mask)
    operators, diagonal, active = images(M, 1)
    virtual, negative_ambient = [0, low_height], [0, -1]
    if cut == 'zero':
        left = endomorphism(L, virtual, negative_ambient)
        right = endomorphism(P)
        maps = [left, left, right, right]
    elif cut == 'infinity':
        maps = [endomorphism(L, reverse_parameter=True),
                endomorphism(P, virtual, negative_ambient, True, True),
                endomorphism(P, negative_ambient, virtual, True),
                endomorphism(L, reverse_parameter=True, transpose=True)]
    else:
        raise ValueError('Unknown place.')
    size = 256
    matrix = [[0]*size for _ in range(size)]
    for state in range(size):
        matrix[state][state] = 1 << sum(value for i, value in enumerate(diagonal) if state >> i & 1)
    for side, (polynomial, _) in enumerate(maps):
        matrix = local_action(matrix, operators[side][0], polynomial, side >= 2)
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            assert value.bit_count() % 2 == int(i == j)
    clearing = sum(value for _, value in maps)
    local_degree_sum = sum(max(value.bit_length()-1 for row in polynomial for value in row)
                           for polynomial, _ in maps)
    degree = max(value.bit_length()-1 for row in matrix for value in row)
    assert degree <= 32*M+local_degree_sum
    return matrix, clearing, degree, local_degree_sum, active


def profile(M, a_mask, b_mask, cut, low_height=1):
    matrix, clearing, degree, local_degree, active = core_matrix(M, a_mask, b_mask, cut, low_height)
    converted = [[(value, 0) for value in row] for row in matrix]
    precision = 128
    while precision <= degree:
        precision *= 2
    attempts = []
    while precision <= 8192:
        try:
            values = MatrixOracle(converted, 256, precision=precision).profile()['orders']
            attempts.append([precision, True])
            break
        except ArithmeticError:
            attempts.append([precision, False])
            precision *= 2
    else:
        raise ArithmeticError('The full-rank affine active profile exceeded its precision budget.')
    assert len(values) == 256 and all(x == round(x) for x in values)
    cleared = sorted(map(round, values))
    result = {'cut': cut, 'orders': [x-clearing for x in cleared], 'cleared_orders': cleared,
              'clearing_power': clearing, 'active_ambient_order': 256, 'active_generic_rank': 256,
              'coefficient_degree': degree, 'coefficient_degree_bound': 32*M+local_degree,
              'local_degree_sum': local_degree, 'precision_attempts': attempts,
              'active_full_mode_indices': active}
    if a_mask == b_mask == 0 and low_height == 1:
        _, diagonal, _ = images(M, 1)
        expected = sorted(sum(value for i, value in enumerate(diagonal) if state >> i & 1)
                          for state in range(256))
        assert result['orders'] == expected
        result['complete_gaussian_subset_sum_control'] = True
    return result


def pressure(M, zero, infinity, low_height):
    zero, infinity = np.asarray(zero, dtype=float), np.asarray(infinity, dtype=float)
    spectator = np.asarray([2*j+1 for j in range(2, M-2)]+[2*M+2*j+1 for j in range(2, M-2)])
    physical = np.asarray([low_height]+[2*j+1 for j in range(1, M)])
    queries = 0

    def objective(log_scale):
        nonlocal queries
        queries += 1
        scale = math.exp(log_scale)
        target = float(logsumexp(-scale*zero)+logsumexp(-scale*infinity)
                       +2*np.logaddexp(0, -scale*spectator).sum())
        equation = lambda rho: target-4*rho*float(np.logaddexp(0, -scale*physical/rho).sum())
        if equation(1e-8) <= 0:
            raise ArithmeticError('The declared positive pressure root is unresolved.')
        high = 2
        while equation(high) > 0:
            high *= 2
            if high > 2**20:
                raise ArithmeticError('The pressure root exceeded its search budget.')
        return float(brentq(equation, 1e-8, high, xtol=1e-13))

    grid = np.linspace(math.log(1e-6), math.log(2), 101)
    values = [objective(x) for x in grid]
    candidates = [{'rho': values[i], 'scale_t': math.exp(float(grid[i])), 'source': 'grid'}
                  for i in (0, len(grid)-1)]
    for i in range(1, len(grid)-1):
        if values[i] <= min(values[i-1], values[i+1]):
            fit = minimize_scalar(objective, bounds=(grid[i-1], grid[i+1]), method='bounded',
                                  options={'xatol': 1e-10, 'maxiter': 80})
            candidates.append({'rho': float(fit.fun), 'scale_t': math.exp(float(fit.x)),
                               'source': 'refined local minimum', 'completed': bool(fit.success)})
    return {**min(candidates, key=lambda value: value['rho']), 'queries': queries,
            'grid_points': len(grid), 'scale_range': [1e-6, 2], 'local_candidates': candidates,
            'global_minimum_certified': False}


def run_case(parameters):
    M, a_mask, b_mask, low_height = parameters
    started = time.monotonic()
    cuts = [profile(M, a_mask, b_mask, cut, low_height) for cut in ('zero', 'infinity')]
    fit = pressure(M, cuts[0]['orders'], cuts[1]['orders'], low_height)
    return {'M': M, 'a_mask': a_mask, 'b_mask': b_mask, 'low_heights': [0, low_height],
            'L_coefficients': affine(a_mask), 'P_coefficients': affine(b_mask),
            'profiles': cuts, 'pressure_fit': fit,
            'physical_alphabet_dimension': 1 << M,
            'spectator_orders': [2*j+1 for j in range(2, M-2)]+[2*M+2*j+1 for j in range(2, M-2)],
            'physical_tail_heights': [2*j+1 for j in range(1, M)],
            'seconds': time.monotonic()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', type=int, nargs='+', default=[16])
    parser.add_argument('--a-masks', type=int, nargs='+', default=list(range(16)))
    parser.add_argument('--b-masks', type=int, nargs='+', default=list(range(16)))
    parser.add_argument('--low-height', type=int, default=1)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.workers <= 4 or args.low_height < 0:
        raise ValueError('Use one through four workers and a nonnegative low height.')
    tasks = [(M, a, b, args.low_height) for M in args.modes for a in args.a_masks for b in args.b_masks]
    assert len(set(tasks)) == len(tasks)
    names = [Path(__file__).name, 'tensor_literal_full_tail_core.py',
             'check_tensor_hot_compounds_f16.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_weighted_smith_f4.py', 'native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'parameters': {k: v for k, v in vars(args).items() if k != 'output'},
              'declared_cases': len(tasks), 'completed': False, 'records': [],
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def retain(record):
        report['records'].append(record)
        report['elapsed_seconds'] = time.monotonic()-started
        temporary = args.output.with_suffix(args.output.suffix+'.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':'))+'\n')
        temporary.replace(args.output)
        print(json.dumps({'completed_cases': len(report['records']), 'declared_cases': len(tasks),
                          'M': record['M'], 'a_mask': record['a_mask'], 'b_mask': record['b_mask'],
                          'rho': record['pressure_fit']['rho'], 'seconds': record['seconds']}), flush=True)

    if args.workers == 1:
        for task in tasks:
            retain(run_case(task))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run_case, task) for task in tasks]
            for future in as_completed(futures):
                retain(future.result())
    assert len(report['records']) == len(tasks)
    report['records'].sort(key=lambda value: (value['M'], value['a_mask'], value['b_mask']))
    report['completed'] = True
    report['elapsed_seconds'] = time.monotonic()-started
    args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
