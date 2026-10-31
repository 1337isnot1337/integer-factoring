#!/usr/bin/env python3
"""Finite numerical screen of proved secondary-height rectangular repairs.

Every candidate has an explicit admissible additive flag. The reported
floating pressures are search results, not exact rate certificates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from check_rectangular_refined_repairs import height, single_shear
from search_rectangular_core_pressure import CoreDefectPressure


def starts_for(M, source=None):
    starts = [([0.]*18+[t/M]) for t in (0.3, 1., 3.)]
    if source:
        supplied = json.loads(source.read_text())
        if isinstance(supplied, list):
            starts.extend(supplied)
        else:
            point = supplied['best']['parameters']
            if len(point) == 5:
                point = point[:4]+[0.]*14+point[-1:]
            starts.append(point)
            starts.append(point[:-1]+[point[-1]*supplied.get('M', M)/M])
    unique = []
    for point in starts:
        if len(point) != 19:
            raise ValueError('Every start must have nineteen parameters.')
        if point not in unique:
            unique.append(point)
    return unique


def evaluate_candidate(M, alpha, beta, candidate, starts):
    begin = time.monotonic()
    spec = None if candidate is None else single_shear(candidate['row'], candidate['column'], candidate['coefficient'])
    if spec is not None and candidate['row'] == 0:
        spec = (*spec, False)
    oracle = CoreDefectPressure(M, alpha, beta, repair_spec=spec)
    attempts, best = [], None
    for index, point in enumerate(starts):
        try:
            value, _, details = oracle.evaluate(point)
            attempts.append({'start': index, 'rho': value})
            if best is None or value < best['rho']:
                best = {'rho': value, 'parameters': point, 'start': index,
                        **{key: details[key] for key in ('type_x_core', 'type_y_core', 'entropy')}}
        except ArithmeticError as error:
            attempts.append({'start': index, 'error': str(error)})
    record = {'candidate': candidate, 'best': best, 'attempts': attempts,
              'seconds': time.monotonic()-begin,
              'exact_polynomial_matrices': [meta['exact_polynomial_matrix'] for meta in oracle.metadata],
              'program_degree_bounds': [meta['program_degree_bound'] for meta in oracle.metadata]}
    if spec:
        record['repair_spec'] = spec
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--M', type=int, default=4)
    parser.add_argument('--alpha', type=int, default=1)
    parser.add_argument('--beta', type=int, default=2)
    parser.add_argument('--starts', type=Path)
    parser.add_argument('--from-screen', type=Path)
    parser.add_argument('--top', type=int, default=8)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--background-only', action='store_true',
                        help='Separately screen the33 additional output shears that change row zero.')
    args = parser.parse_args()
    begin = time.monotonic()
    starts = starts_for(args.M, args.starts)
    if args.from_screen:
        previous = json.loads(args.from_screen.read_text())
        candidates = [entry['candidate'] for entry in sorted(
            (entry for entry in previous['candidates'] if entry['best'] is not None),
            key=lambda entry: entry['best']['rho'])[:args.top]]
    else:
        candidates = [{'row': row, 'column': column, 'coefficient': coefficient,
                       'height_relation': 'equal' if height(row) == height(column) else 'strict'}
                      for row in ((0,) if args.background_only else range(1, 12)) for column in range(12)
                      if row != column and height(row) <= height(column) for coefficient in (1, 2, 3)]
        candidates.sort(key=lambda candidate: (candidate['height_relation'] != 'equal', candidate['row'], candidate['column'], candidate['coefficient']))
    result = {'status': 'finite fixed-gauge floating screen; no exact rate claimed',
              'M': args.M, 'alpha': args.alpha, 'beta': args.beta,
              'background_row_fixed': not args.background_only,
              'quartic_terms': [[1, 0]], 'starts': starts,
              'planned_candidates': len(candidates), 'baseline': None, 'candidates': []}
    prior_seconds = 0.
    if args.resume and args.output.exists():
        result = json.loads(args.output.read_text())
        if result['M'] != args.M or result['alpha'] != args.alpha or result['beta'] != args.beta or result['starts'] != starts:
            raise ValueError('Resume parameters differ from the checkpoint.')
        prior_seconds = result.get('elapsed_seconds', 0.)
    def save():
        result['elapsed_seconds'] = prior_seconds+time.monotonic()-begin
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    if result['baseline'] is None:
        result['baseline'] = evaluate_candidate(args.M, args.alpha, args.beta, None, starts)
        save()
        print(json.dumps({'baseline_rho': result['baseline']['best']['rho']}), flush=True)
    done = {(entry['candidate']['row'], entry['candidate']['column'], entry['candidate']['coefficient']) for entry in result['candidates']}
    best = min([entry['best']['rho'] for entry in result['candidates'] if entry['best']] + [result['baseline']['best']['rho']])
    for candidate in candidates:
        key = (candidate['row'], candidate['column'], candidate['coefficient'])
        if key in done:
            continue
        record = evaluate_candidate(args.M, args.alpha, args.beta, candidate, starts)
        result['candidates'].append(record)
        save()
        improved = record['best'] is not None and record['best']['rho'] < best-1e-10
        if improved:
            best = record['best']['rho']
        if improved or len(result['candidates']) % 12 == 0:
            print(json.dumps({'completed': len(result['candidates']), 'best_rho': best,
                              'candidate': candidate if improved else None,
                              'elapsed_seconds': result['elapsed_seconds']}), flush=True)
    result['completed'] = len(result['candidates']) == len(candidates)
    save()
    print(json.dumps({'completed': len(result['candidates']), 'best_rho': best,
                      'elapsed_seconds': result['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
