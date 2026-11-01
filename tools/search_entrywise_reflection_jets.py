#!/usr/bin/env python3
"""Fast exact higher-jet searches using entrywise polynomial interpolation."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path
import time

import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure
from entrywise_reflection_constraints import field_dot, hasse_row, solve_target
from f16_identity_operator_profiles import inverse_matrix
from physical_split_crossings import F4_PRODUCT, identity_split
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize, digest
from screen_affine_finite_places import bitplane_rank, rows
from search_fused_mixed_weighted import atomic_write
from search_physical_code_gauges import start_params


def bell_values(matrix, cut):
    inverse = inverse_matrix(matrix)
    result = []
    for entry in range(256):
        x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
        first, second = (matrix[x][y], inverse[u][v]) if cut == 0 else (matrix[x][u], inverse[y][v])
        result.append(F4_PRODUCT[first][second])
    return result


def bell_jets(matrix, cut, point, order, trace_power=16):
    values = bell_values(matrix, cut)
    scalar = (1 | (1 << trace_power), 0)
    coefficients = []
    for jet in range(order):
        multiplier = field_dot(hasse_row(point, jet, trace_power), scalar)
        coefficients.append([F4_PRODUCT[multiplier][value] for value in values])
    return {'point': point, 'values': coefficients}


def basis_catalogue(base):
    _, factors = bell_normalize(base)
    evaluation = []
    for letter in range(4):
        row, value = [], 1
        for _ in range(4):
            row.append(value)
            value = F4_PRODUCT[value][letter]
        evaluation.append(row)
    return {'identity': identity_split(4), 'source_Bell': factors['first_factor'],
            'source_Bell_inverse': factors['inverse_factor'],
            'Pascal': [[int(j & i == j) for j in range(4)] for i in range(4)],
            'evaluation': evaluation}


def tasks(matrices):
    protections = [(8, 16), (12, 20), (0, 0), (12, 0), (16, 0), (0, 16),
                   (4, 12), (8, 12), (12, 16), (16, 24)]
    pairs = [('source_Bell', 'source_Bell'), ('source_Bell', 'identity'),
             ('identity', 'source_Bell'), ('identity', 'identity'),
             ('source_Bell', 'source_Bell_inverse'), ('Pascal', 'Pascal'),
             ('evaluation', 'evaluation')]
    for order in (4, 6, 2, 8, 10, 12):
        for cap in (16, 24, 20, 22):
            for cutoffs in protections:
                for mode in ('exact', 'zero-pattern'):
                    for row_basis, column_basis in pairs:
                        yield {'kind': 'Bell jets', 'cap': cap, 'cutoffs': cutoffs, 'mode': mode,
                               'order': order, 'row_basis': row_basis, 'column_basis': column_basis,
                               'conditions': [], 'prescribed_jets': [
                                   bell_jets(matrices[row_basis], 0, 2, order),
                                   bell_jets(matrices[column_basis], 1, 3, order)]}
    supports = [(first, tuple(i for i in range(4) if i not in first)) for first in combinations(range(4), 2)]
    for order in (4, 8, 2, 6, 12):
        for cap in (16, 24, 20, 22):
            for cutoffs in protections:
                for mode in ('exact', 'zero-pattern'):
                    for first, second in supports:
                        yield {'kind': 'complementary rectangles', 'cap': cap, 'cutoffs': cutoffs,
                               'mode': mode, 'order': order,
                               'conditions': [(2, 0, first, second, order), (3, 1, first, second, order)],
                               'prescribed_jets': []}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--targets', type=int, default=100000)
    args = parser.parse_args()
    base = tensor(4)
    matrices = basis_catalogue(base)
    names = (Path(__file__).name, 'entrywise_reflection_constraints.py',
             'solve_companion_affine_degree.py', 'reflection_dual_rail_gaussian.py',
             'asymmetric_crossing_pressure.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Exact entrywise higher-jet constraints and full paired traces; finite ranks and numerical pressure are diagnostics, not a factoring certificate.',
              'arguments': vars(args) | {'output': str(args.output)}, 'matrices': matrices,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'cases': [], 'retained_seeds': [], 'feasible_targets': 0, 'unique_seeds': 0,
              'unresolved_pressure_calls': 0, 'completed': False}
    began, seen = time.monotonic(), set()
    natural = start_params(np.asarray([0, 12, 4, 16]) * .2326459924883792, .2326459924883792)
    exhausted = True
    for index, target in enumerate(tasks(matrices)):
        if index >= args.targets or time.monotonic() - began >= args.seconds:
            exhausted = False
            break
        trial_start = time.monotonic()
        result = solve_target(base, target['cap'], target['conditions'], *target['cutoffs'],
                              trace_power=16, protection_mode=target['mode'],
                              prescribed_jets=target['prescribed_jets'])
        seed = result.pop('seed', None)
        record = {key: value for key, value in target.items() if key != 'prescribed_jets'}
        record.update(result, target_index=index, exact_solve_seconds=time.monotonic() - trial_start)
        if seed is not None:
            report['feasible_targets'] += 1
            key = digest(seed['values'])
            record['seed_sha256'] = key
            if key not in seen:
                seen.add(key)
                finite_ranks = [[bitplane_rank(rows(seed['values'], cut, point)) for cut in (0, 1)]
                                for point in (2, 3)]
                if target['kind'] == 'Bell jets':
                    assert finite_ranks[0][0] == finite_ranks[1][1] == 1
                pressure = AsymmetricPressure(seed, 1)
                scores = []
                for multiple in (1, .5, .2):
                    try:
                        scores.append(pressure.evaluate(np.r_[natural * multiple / math.sqrt(2), math.pi / 4],
                                                        gradient=False)[0] / math.sqrt(2))
                    except (ArithmeticError, ValueError):
                        report['unresolved_pressure_calls'] += 1
                best = min(scores) if scores else None
                record.update(actual_degree=seed['degree'], finite_ranks=finite_ranks, natural_rho=best)
                report['retained_seeds'].append(record | {'seed': seed})
                retained = report['retained_seeds']
                low_cost = sorted(retained, key=lambda row: row.get('natural_rho') or 1e6)[:40]
                high_order = sorted(retained, key=lambda row: (-row['order'], row.get('natural_rho') or 1e6))[:40]
                protected = sorted(retained, key=lambda row: (-sum(row['cutoffs']), row.get('natural_rho') or 1e6))[:40]
                report['retained_seeds'] = list({row['seed_sha256']: row for row in low_cost + high_order + protected}.values())
        report['cases'].append(record)
        if index % 8 == 7 or seed is not None:
            report.update(elapsed_seconds=time.monotonic() - began,
                          completed_targets=len(report['cases']), unique_seeds=len(seen))
            atomic_write(args.output, report)
            print(json.dumps({'targets': len(report['cases']), 'feasible': report['feasible_targets'],
                              'unique': len(seen), 'seconds': report['elapsed_seconds'],
                              'last_kind': target['kind'], 'last_order': target['order']}), flush=True)
    report.update(completed=True, elapsed_seconds=time.monotonic() - began,
                  completed_targets=len(report['cases']), unique_seeds=len(seen),
                  target_queue_exhausted=exhausted)
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
