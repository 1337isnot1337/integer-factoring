#!/usr/bin/env python3
"""Cached fixed-gauge screen of four endpoint-vanishing extreme quartics."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from check_rectangular_extra_quartics import PATTERNS
from rectangular_coordinate_cache import cached_coordinates, extra_quartics_from_bases
from screen_rectangular_polynomial_pressure import supplied_starts
from search_rectangular_core_pressure import CoreDefectPressure


def full_parameters(point):
    return point[:11]+[0.]*3+point[11:18]+[0.]*3+point[18:]


def evaluate(M, alpha, beta, full, data, compiled, starts):
    started = time.monotonic()
    oracle = CoreDefectPressure(M, alpha, beta, full_core=full, compiled=compiled, source_data=data)
    attempts, best = [], None
    for index, point in enumerate(starts):
        try:
            rho, _, details = oracle.evaluate(point)
            attempts.append({'start': index, 'rho': rho})
            if best is None or rho < best['rho']:
                best = {'rho': rho, 'parameters': point, 'start': index,
                        **{key: details[key] for key in ('type_x_core', 'type_y_core', 'entropy')}}
        except ArithmeticError as error:
            attempts.append({'start': index, 'error': str(error)})
    return {'best': best, 'attempts': attempts, 'seconds': time.monotonic()-started,
            'exact_polynomial_matrices': [entry[2]['exact_polynomial_matrix'] for entry in compiled],
            'profile_precision': [entry[2]['profile_precision'] for entry in compiled]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--M', type=int, default=8)
    parser.add_argument('--alpha', type=int, default=1)
    parser.add_argument('--beta', type=int, default=2)
    parser.add_argument('--starts', type=Path, action='append', required=True)
    parser.add_argument('--patterns', type=int, nargs='+', default=[3, 0, 1, 2])
    parser.add_argument('--ells', type=int, nargs='+')
    parser.add_argument('--coefficients', type=int, nargs='+', default=[1, 2, 3])
    parser.add_argument('--cache-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--stop-file', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    points = supplied_starts(args.M, args.starts)
    if any(label not in range(4) for label in args.patterns) or any(value not in (1, 2, 3) for value in args.coefficients):
        raise ValueError('Use pattern labels zero to three and nonzero F4 coefficients.')
    settings = []
    for label in args.patterns:
        center = (args.M-1)*(1 if label in (1, 3) else -1)
        ells = args.ells or [0, center, center-args.M//2, center+args.M//2]
        settings.append({'pattern': label, 'ells': ells})
    result = {'status': 'finite fixed-gauge extra-quartic floating screen; no exact pressure certificate',
              'M': args.M, 'alpha': args.alpha, 'beta': args.beta, 'starts': points,
              'settings': settings, 'coefficients': args.coefficients,
              'planned_candidates': sum(len(item['ells'])*len(args.coefficients) for item in settings),
              'baselines': {}, 'candidates': []}
    prior_seconds = 0
    if args.resume and args.output.exists():
        old = json.loads(args.output.read_text())
        for key in ('M', 'alpha', 'beta', 'starts', 'settings', 'coefficients'):
            if old[key] != result[key]:
                raise ValueError(f'Resume parameter differs: {key}.')
        result = old
        prior_seconds = old.get('elapsed_seconds', 0)
    def save():
        result['elapsed_seconds'] = prior_seconds+time.monotonic()-started
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    done = {(item['pattern'], item['coefficient'], item['ell']) for item in result['candidates']}
    prepared = {}
    for setting in settings:
        label, ells = setting['pattern'], setting['ells']
        if all((label, coefficient, ell) in done for coefficient in args.coefficients for ell in ells):
            continue
        full = label == 2
        name = 'full' if full else 'partial'
        starts = [full_parameters(point) for point in points] if full else points
        if name not in prepared:
            data, base, cache = cached_coordinates(args.cache_directory, args.M, args.alpha, args.beta, full_core=full)
            prepared[name] = data, base
            if name not in result['baselines']:
                result['baselines'][name] = evaluate(args.M, args.alpha, args.beta, full, data, base, starts)
                result['baselines'][name]['cache'] = cache
                save()
                print(json.dumps({'basis': name, 'baseline_rho': result['baselines'][name]['best']['rho']}), flush=True)
        data, base = prepared[name]
        modes = [mode*(args.M-1) for mode in PATTERNS[label]]
        _, variation, cache = cached_coordinates(args.cache_directory, args.M, args.alpha, args.beta,
                                                  full_core=full, extra_quartics=[(*modes, [(1, 0)], 0)])
        print(json.dumps({'pattern': label, 'basis_cache_ready': True,
                          'matrix_orders': [len(item[0]) for item in variation]}), flush=True)
        best = min([item['best']['rho'] for item in result['baselines'].values()]+[
            item['best']['rho'] for item in result['candidates'] if item['best']])
        for ell in ells:
            for coefficient in args.coefficients:
                if (label, coefficient, ell) in done:
                    continue
                if args.stop_file and args.stop_file.exists():
                    result['stopped_by_file'] = True
                    save()
                    return
                compiled = extra_quartics_from_bases(base, [(variation, [(coefficient, ell)])])
                record = evaluate(args.M, args.alpha, args.beta, full, data, compiled, starts)
                record.update({'pattern': label, 'modes': modes, 'coefficient': coefficient,
                               'ell': ell, 'full_core': full, 'basis_cache': cache['digest']})
                result['candidates'].append(record)
                save()
                improved = record['best'] is not None and record['best']['rho'] < best-1e-10
                if improved:
                    best = record['best']['rho']
                print(json.dumps({'completed': len(result['candidates']), 'best_rho': best,
                                  'candidate_rho': record['best']['rho'] if record['best'] else None,
                                  'candidate': [label, coefficient, ell],
                                  'elapsed_seconds': result['elapsed_seconds']}), flush=True)
    result['completed'] = len(result['candidates']) == result['planned_candidates']
    save()


if __name__ == '__main__':
    main()
