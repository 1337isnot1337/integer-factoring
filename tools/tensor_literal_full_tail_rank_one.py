#!/usr/bin/env python3
"""Complete F4 rank-one low-mode frames with a full Gaussian hot tail.

The physical maps may mix the vacuum and occupied lowest mode. A fixed
coordinate inclusion embeds the actual weighted rectangular maps as low
endomorphisms, retaining all Laurent clearing shifts. The literal 8-mode
active matrix has order256 and generic rank64, independently of M.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

import tensor_hot_compounds_f16 as field
from native_matrix_smith_f4 import MatrixOracle
from tensor_literal_full_tail_core import images, left_action, right_action, number_projector
from tensor_literal_rank_one_cell_search import all_cases, orbit


def maps(case):
    infinity, a, b = case
    if infinity:
        p, left = [0, 1], [b, 1]
    else:
        p, left = [1, a], [1 ^ field.mul(a, b), b]
    assert field.mul(p[0], left[0]) ^ field.mul(p[1], left[1]) == 1
    assert all(0 <= c < 4 for c in p+left)
    return p, left


def apply_endomorphism(matrix, operators, terms, right):
    size = len(matrix[0])
    out = [[[0]*size for _ in range(size)] for _ in range(2)]
    for kind, coefficient, power in terms:
        if not coefficient:
            continue
        assert power >= 0
        if kind in ('00', '11'):
            transformed = [number_projector(plane, operators, int(kind == '11'), right) for plane in matrix]
        else:
            operator = operators[0 if kind == '10' else 1]
            action = right_action if right else left_action
            transformed = [action(plane, operator) for plane in matrix]
        for i in range(size):
            for j in range(size):
                a, b = transformed[0][i][j], transformed[1][i][j]
                if coefficient == 1:
                    low, high = a, b
                elif coefficient == 2:
                    low, high = b, a ^ b
                elif coefficient == 3:
                    low, high = a ^ b, a
                else:
                    raise ValueError('Expected an F4 coefficient.')
                out[0][i][j] ^= low << power
                out[1][i][j] ^= high << power
    return out


def core_matrix(M, case, cut):
    p, left = maps(case)
    (rx, ry, cv, cu), diagonal, active = images(M, 1)
    size = 1 << len(diagonal)
    matrix = [[[0]*size for _ in range(size)] for _ in range(2)]
    for state in range(size):
        matrix[0][state][state] = 1 << sum(order for i, order in enumerate(diagonal) if state >> i & 1)
    lpole, ppole = int(left[1] != 0), int(p[1] != 0)
    left_plain = [('00', left[0], 0), ('01', left[1], 0)]
    right_plain = [('00', p[0], 0), ('10', p[1], 0)]
    left_weighted = [('00', left[0], lpole), ('01', left[1], max(0, lpole-1))]
    p_row_weighted = [('00', p[0], ppole), ('01', p[1], max(0, ppole-1))]
    p_column_weighted = [('00', p[0], ppole), ('10', p[1], max(0, ppole-1))]
    l_column_plain = [('00', left[0], 0), ('10', left[1], 0)]
    if cut == 'zero':
        specifications = [(rx[0], left_weighted, False), (ry[0], left_weighted, False),
                          (cv[0], right_plain, True), (cu[0], right_plain, True)]
        clearing = 2*lpole
    elif cut == 'infinity':
        specifications = [(rx[0], left_plain, False), (ry[0], p_row_weighted, False),
                          (cv[0], p_column_weighted, True), (cu[0], l_column_plain, True)]
        clearing = 2*ppole
    else:
        raise ValueError('Unknown cut.')
    for operators, terms, right in specifications:
        matrix = apply_endomorphism(matrix, operators, terms, right)
    degree = max(value.bit_length()-1 for plane in matrix for row in plane for value in row)
    assert degree <= 32*M+2
    endpoint_rank = 0
    for i, j in ((i, j) for i in range(size) for j in range(size)):
        expected = int(i == j and not (i & 1) and not (i & (1 << 4)))
        assert matrix[0][i][j].bit_count() % 2 == expected
        assert matrix[1][i][j].bit_count() % 2 == 0
        if i == j:
            endpoint_rank += expected
    assert endpoint_rank == 64
    return matrix, clearing, degree, active


def core_profile(M, case, cut):
    matrix, clearing, degree, active = core_matrix(M, case, cut)
    size = len(matrix[0])
    converted = [[(matrix[0][i][j], matrix[1][i][j]) for j in range(size)] for i in range(size)]
    precision = 128
    while precision <= degree:
        precision *= 2
    attempts = []
    while precision <= 8192:
        try:
            values = MatrixOracle(converted, 64, precision=precision).profile()['orders']
            attempts.append([precision, True])
            break
        except ArithmeticError:
            attempts.append([precision, False])
            precision *= 2
    else:
        raise ArithmeticError('The complete active profile exceeded the series budget.')
    assert len(values) == 64 and all(value == round(value) for value in values)
    cleared = sorted(map(round, values))
    return {'cut': cut, 'orders': [value-clearing for value in cleared],
            'cleared_orders': cleared, 'clearing_power': clearing,
            'active_ambient_order': 256, 'active_generic_rank': 64,
            'coefficient_degree': degree, 'coefficient_degree_bound': 32*M+2,
            'precision_attempts': attempts, 'active_full_mode_indices': active}


def pressure(M, zero, infinity):
    zero, infinity = np.asarray(zero, dtype=float), np.asarray(infinity, dtype=float)
    spectators = np.asarray([2*j+1 for j in range(2, M-2)]
                           + [2*M+2*j+1 for j in range(2, M-2)], dtype=float)
    physical = np.asarray([2*j+1 for j in range(1, M)], dtype=float)
    queries = 0

    def objective(log_scale):
        nonlocal queries
        queries += 1
        scale = math.exp(log_scale)
        target = float(logsumexp(-scale*zero)+logsumexp(-scale*infinity)
                       +2*np.logaddexp(0, -scale*spectators).sum())
        equation = lambda rho: target-4*rho*float(np.logaddexp(0, -scale*physical/rho).sum())
        if equation(1e-8) <= 0:
            raise ArithmeticError('The positive pressure root is not resolved at this scale.')
        high = 2
        while equation(high) > 0:
            high *= 2
            if high > 2**20:
                raise ArithmeticError('Pressure root exceeded the numerical search budget.')
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
    return {**min(candidates, key=lambda c: c['rho']), 'queries': queries,
            'grid_points': len(grid), 'scale_range': [1e-6, 2], 'local_candidates': candidates,
            'global_minimum_certified': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', type=int, nargs='+', default=[4, 8, 16, 32, 64, 128])
    parser.add_argument('--case', type=int, nargs=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cases = [tuple(args.case)] if args.case is not None else all_cases(4)
    names = [Path(__file__).name, 'tensor_literal_full_tail_core.py',
             'tensor_literal_rank_one_cell_search.py', 'check_tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16.py', 'tensor_orthogonal_scaling.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'native_weighted_smith_f4.py', 'native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'field_size': 4, 'modes': args.modes,
              'cases': cases, 'complete_projective_catalogue': args.case is None,
              'assigned_low_virtual_height': 0, 'completed': False,
              'in_progress': None, 'records': [],
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        report['elapsed_seconds'] = time.monotonic()-started
        temporary = args.output.with_suffix(args.output.suffix+'.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':'))+'\n')
        temporary.replace(args.output)

    for M in args.modes:
        for case in cases:
            before = time.monotonic()
            p, left = maps(case)
            cuts = []
            for cut in ('zero', 'infinity'):
                report['in_progress'] = {'M': M, 'case': case, 'cut': cut}
                save()
                cuts.append(core_profile(M, case, cut))
            fit = pressure(M, cuts[0]['orders'], cuts[1]['orders'])
            report['records'].append({'M': M, 'case': case, 'frobenius_orbit': orbit(case),
                                      'embedding': p, 'retraction': left, 'profiles': cuts,
                                      'spectator_orders': [2*j+1 for j in range(2, M-2)]
                                                         + [2*M+2*j+1 for j in range(2, M-2)],
                                      'physical_tail_heights': [2*j+1 for j in range(1, M)],
                                      'physical_alphabet_dimension': 1 << (M-1),
                                      'pressure_fit': fit, 'seconds': time.monotonic()-before})
            report['in_progress'] = None
            save()
            print(json.dumps({'M': M, 'case': case, 'rho': fit['rho'],
                              'scale_t': fit['scale_t'], 'seconds': time.monotonic()-before}), flush=True)
    report['completed'] = True
    save()


if __name__ == '__main__':
    main()
