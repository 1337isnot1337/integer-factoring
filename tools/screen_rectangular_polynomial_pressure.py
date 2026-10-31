#!/usr/bin/env python3
"""Focused fixed-gauge pressure screen of proved polynomial output shears.

Each candidate is I+gamma*w^ell*(1+w)E_rs, hence fixes the ramified
endpoint. Floating pressure values are search evidence, not rate proofs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from check_rectangular_identity_transports import IDENTITY
from check_rectangular_polynomial_shears import polynomial_coordinate_shear
from check_rectangular_refined_repairs import height
from rectangular_coordinate_transforms import transform_coordinate_matrix
from screen_rectangular_refined_pressure import prepare_coordinates, starts_for
from search_rectangular_core_pressure import CoreDefectPressure


DEFAULT_DIRECTIONS = ((0, 1), (0, 4), (0, 5), (0, 10), (0, 11),
                      (3, 8), (8, 3), (10, 11), (9, 8),
                      (10, 1), (11, 0), (8, 1), (2, 0))


def supplied_starts(M, paths):
    neutral = starts_for(M)
    result = []
    for path in paths:
        for point in starts_for(M, path):
            if point not in neutral and point not in result:
                result.append(point)
    if not result:
        raise ValueError('Supply at least one non-neutral nineteen-parameter start.')
    return result


def evaluate(M, alpha, beta, data, compiled, starts):
    started = time.monotonic()
    oracle = CoreDefectPressure(M, alpha, beta, compiled=compiled, source_data=data)
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
    parser.add_argument('--directions', type=Path)
    parser.add_argument('--ells', type=int, nargs='+')
    parser.add_argument('--coefficients', type=int, nargs='+', default=[1, 2, 3])
    parser.add_argument('--cache-directory', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--stop-file', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    starts = supplied_starts(args.M, args.starts)
    directions = json.loads(args.directions.read_text()) if args.directions else list(DEFAULT_DIRECTIONS)
    ells = args.ells or [0, -args.M//2, args.M//2]
    if any(value not in (1, 2, 3) for value in args.coefficients):
        raise ValueError('Coefficients must be nonzero elements of F4.')
    data, bases, cache = prepare_coordinates(args.M, args.alpha, args.beta, args.cache_directory)
    result = {'status': 'finite fixed-gauge polynomial-shear screen; no exact pressure certificate',
              'M': args.M, 'alpha': args.alpha, 'beta': args.beta, 'starts': starts,
              'directions': [list(pair) for pair in directions], 'ells': ells,
              'coefficients': args.coefficients, 'coordinate_cache': cache,
              'planned_candidates': len(directions)*len(ells)*len(args.coefficients),
              'baseline': None, 'candidates': []}
    previous_seconds = 0
    if args.resume and args.output.exists():
        old = json.loads(args.output.read_text())
        for key in ('M', 'alpha', 'beta', 'starts', 'directions', 'ells', 'coefficients'):
            if old[key] != result[key]:
                raise ValueError(f'Resume parameter differs: {key}.')
        result = old
        previous_seconds = result.get('elapsed_seconds', 0)
    def save():
        result['elapsed_seconds'] = previous_seconds+time.monotonic()-started
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    if result['baseline'] is None:
        result['baseline'] = evaluate(args.M, args.alpha, args.beta, data, bases, starts)
        save()
        print(json.dumps({'baseline_rho': result['baseline']['best']['rho']}), flush=True)
    done = {(entry['row'], entry['column'], entry['coefficient'], entry['ell']) for entry in result['candidates']}
    best = min([result['baseline']['best']['rho']]+[entry['best']['rho'] for entry in result['candidates'] if entry['best']])
    for row, column in directions:
        if all((row, column, coefficient, ell) in done for ell in ells for coefficient in args.coefficients):
            continue
        gate = [line[:] for line in IDENTITY]
        gate[row][column] = 1
        steps = [transform_coordinate_matrix(base, source, V=gate) for base, source in zip(bases, data)]
        for ell in ells:
            for coefficient in args.coefficients:
                if (row, column, coefficient, ell) in done:
                    continue
                if args.stop_file and args.stop_file.exists():
                    result['stopped_by_file'] = True
                    save()
                    return
                terms = [(coefficient, ell), (coefficient, ell+1)]
                compiled = [polynomial_coordinate_shear(base, source, row, column, terms, unit_step=step)
                            for base, source, step in zip(bases, data, steps)]
                record = evaluate(args.M, args.alpha, args.beta, data, compiled, starts)
                record.update({'row': row, 'column': column, 'coefficient': coefficient,
                               'ell': ell, 'terms': terms, 'height_change': height(row)-height(column)})
                result['candidates'].append(record)
                save()
                improved = record['best'] is not None and record['best']['rho'] < best-1e-10
                if improved:
                    best = record['best']['rho']
                if improved or len(result['candidates']) % 6 == 0:
                    print(json.dumps({'completed': len(result['candidates']), 'best_rho': best,
                                      'candidate': [row, column, coefficient, ell] if improved else None,
                                      'elapsed_seconds': result['elapsed_seconds']}), flush=True)
    result['completed'] = len(result['candidates']) == result['planned_candidates']
    save()
    print(json.dumps({'completed': len(result['candidates']), 'best_rho': best,
                      'elapsed_seconds': result['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
