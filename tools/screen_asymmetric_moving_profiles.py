#!/usr/bin/env python3
"""Finite asymmetric supporting-line screen of frozen moving-frame profiles.

Balanced pressure exclusions do not exclude these supporting weights.
All fit values here are numerical; a strict improvement needs rational
replay and the complete asymmetric transfer.
"""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize

ROOT = Path(__file__).resolve().parents[1]
ENERGIES = np.array([0, 4, 12, 16, 20, 24, 28, 32], float)


def log_count(gaps, scale):
    units = int(np.count_nonzero(gaps == 0))
    return math.log(units) + math.log1p(float(np.exp(-scale * gaps[gaps > 0]).sum()) / units)


def support(sides, tilt, angle):
    a, b = math.cos(angle), math.sin(angle)
    minima = [min(side) for side in sides]
    gaps = [np.array(side, float) - low for side, low in zip(sides, minima)]
    raw = a * log_count(gaps[0], tilt / a) + b * log_count(gaps[1], tilt / b) - tilt * sum(minima)
    def equation(rate):
        return raw - 2 * rate * math.log1p(float(np.exp(-tilt * ENERGIES[1:] / rate).sum()))
    upper = 4.
    while equation(upper) > 0:
        upper *= 2
    return float(brentq(equation, 1e-10, upper, xtol=1e-12))


def native_controls():
    from asymmetric_crossing_pressure import AsymmetricPressure
    from check_coherent_hot_boundary import logical, natural_profiles
    from moving_physical_frames import compress_frame
    from restrict_implicit_reflection import restricted_tensor
    from search_moving_interior_crossings import pressure
    from search_physical_code_gauges import start_params
    raw = restricted_tensor(16, list(map(logical, [0, 1, 2, 3, 4, 5, 8, 6, 9])))
    shifts = [(32 - int(value)) // 4 for value in ENERGIES]
    errors = []
    for am, bm in ((128, 128), (170, 85)):
        a = [(int(am >> i & 1), 0) for i in range(8)]
        b = [(int(bm >> i & 1), 0) for i in range(8)]
        seed = compress_frame(raw, a, b, shifts, [-value for value in shifts])
        sides = natural_profiles(seed, ENERGIES.tolist())
        oracle = AsymmetricPressure(seed)
        for angle in (.4, math.pi / 4, 1.1):
            tilt = .14
            actual, _ = oracle.evaluate(np.r_[start_params(ENERGIES * tilt, tilt), angle], gradient=False)
            errors.append(abs(actual - support(sides, tilt, angle)))
        balanced = pressure(sides, ENERGIES.tolist())
        errors.append(abs(support(sides, balanced['scale'] / math.sqrt(2), math.pi / 4)
                          - math.sqrt(2) * balanced['rho']))
    assert max(errors) < 1e-9
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    errors = native_controls()
    path = args.archive or ROOT / 'evidence/moving-frame-height1-prefix16384.jsonl.gz'
    payload = gzip.decompress(path.read_bytes())
    if path.name.endswith('.jsonl.gz'):
        records = [json.loads(line) for line in payload.splitlines()]
    else:
        structured = json.loads(payload)
        if not structured.get('completed'):
            raise ValueError('The supplied structured search archive is incomplete.')
        records = structured['records']
    rational_u = [(1, 64), (1, 32), (1, 16), (1, 8), (1, 4), (3, 8), (1, 2),
                  (5, 8), (3, 4), (7, 8), (15, 16), (31, 32), (63, 64)]
    angles = [2 * math.atan(n / d) for n, d in rational_u] + [math.pi / 4]
    tilts = np.geomspace(.002, 1.5, 48)
    grid_angle, grid_tilt = np.meshgrid(angles, tilts)
    aa, bb, tt = np.cos(grid_angle.ravel()), np.sin(grid_angle.ravel()), grid_tilt.ravel()
    physical_tail = np.exp(-tt[:, None] * ENERGIES[None, 1:])
    physical_log = np.log1p(physical_tail.sum(axis=1))
    physical_mean = (physical_tail @ ENERGIES[1:]) / (1 + physical_tail.sum(axis=1))
    entropy = physical_log + tt * physical_mean
    report = {'scope': __doc__, 'completed': False, 'profile_cases': len(records),
              'rational_half_angle_grid': rational_u, 'balanced_angle_also_included': True,
              'tilt_grid': tilts.tolist(), 'native_and_balanced_control_errors': errors,
              'profile_archive_sha256': sha256(path.read_bytes()).hexdigest(),
              'profile_archive_name': path.name,
              'source_sha256': {f'tools/{name}': sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in ('screen_asymmetric_moving_profiles.py', 'asymmetric_crossing_pressure.py',
                                             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
                                             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
                                             'moving_physical_frames.py', 'moving_boundary_crossings.py',
                                             'check_coherent_hot_boundary.py', 'restrict_implicit_reflection.py',
                                             'search_moving_interior_crossings.py', 'search_physical_code_gauges.py')},
              'grid_evaluations': 0, 'negative_grid_points': 0, 'top_grid_cases': []}
    began = time.monotonic()
    top = []
    for begin in range(0, len(records), 16):
        chunk = records[begin:begin + 16]
        sides = np.array([record['orders_infinity_zero'] for record in chunk], float)
        minima = sides.min(axis=2)
        gaps = sides - minima[:, :, None]
        units = (gaps == 0).sum(axis=2)
        logs = []
        for side, weights in enumerate((aa, bb)):
            powers = np.where(gaps[:, side, None, :] > 0,
                              -gaps[:, side, None, :] * (tt / weights)[None, :, None], -np.inf)
            logs.append(np.log(units[:, side, None])
                        + np.log1p(np.exp(powers).sum(axis=2) / units[:, side, None]))
        gap = aa * logs[0] + bb * logs[1] - minima.sum(axis=1)[:, None] * tt - 2 * physical_log
        ranking = gap / (2 * entropy)
        report['negative_grid_points'] += int((gap < 0).sum())
        report['grid_evaluations'] += gap.size
        for offset, record in enumerate(chunk):
            position = int(np.argmin(ranking[offset]))
            recipe = {name: record[name] for name in ('a_mask', 'b_mask', 'interior_index', 'interior_letter', 'a_code', 'b_code')
                      if name in record}
            top.append({'index': begin + offset, 'recipe': recipe,
                        'tilt': float(tt[position]), 'angle': float(grid_angle.ravel()[position]),
                        'gap_at_rate_one': float(gap[offset, position]),
                        'linearized_rate_excess': float(ranking[offset, position])})
        top = sorted(top, key=lambda value: value['linearized_rate_excess'])[:20]
        if (begin + len(chunk)) % 2048 == 0:
            print(json.dumps({'screened_profiles': begin + len(chunk),
                              'best_grid_excess': top[0]['linearized_rate_excess'],
                              'negative_grid_points': report['negative_grid_points']}), flush=True)
    report['top_grid_cases'] = top
    refined = []
    for candidate in top:
        sides = records[candidate['index']]['orders_infinity_zero']
        def objective(values):
            return support(sides, math.exp(float(values[0])), float(values[1]))
        start = [math.log(candidate['tilt']), candidate['angle']]
        initial = objective(start)
        fit = minimize(objective, start, method='Powell',
                       bounds=[(math.log(.0001), math.log(3)), (.001, math.pi / 2 - .001)],
                       options={'maxiter': 80, 'ftol': 1e-10, 'xtol': 1e-8})
        use_initial = initial <= fit.fun
        selected = start if use_initial else fit.x
        refined.append({**candidate, 'support_value': float(min(initial, fit.fun)),
                        'initial_support_value': initial, 'fitted_support_value': float(fit.fun),
                        'used_initial_candidate': use_initial,
                        'tilt': math.exp(float(selected[0])),
                        'angle': float(selected[1]), 'optimizer_success': bool(fit.success),
                        'optimizer_message': str(fit.message), 'evaluations': int(fit.nfev)})
    report.update(completed=True, refined=refined, elapsed_seconds=time.monotonic() - began,
                  best_refined=min(refined, key=lambda value: value['support_value']))
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'completed': True, 'profile_cases': len(records),
                      'best_refined_support': report['best_refined']['support_value'],
                      'negative_grid_points': report['negative_grid_points'],
                      'elapsed_seconds': report['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
