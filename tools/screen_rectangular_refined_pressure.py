#!/usr/bin/env python3
"""Finite numerical screen of proved secondary-height rectangular repairs.

Every candidate has an explicit admissible additive flag. The reported
floating pressures are search results, not exact rate certificates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import pickle
import time

from check_rectangular_defect_core_basis import data_for
from check_rectangular_refined_repairs import height, single_shear
from rectangular_coordinate_transforms import transform_coordinate_matrix
from rectangular_defect_coordinates import coordinate_matrix
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


def prepare_coordinates(M, alpha, beta, cache_directory=None):
    """Build once, or reuse this program's explicitly selected local cache."""
    directory = Path(__file__).resolve().parent
    dependencies = ('check_rectangular_f4_repair.py', 'check_rectangular_defect_support.py',
                    'check_rectangular_defect_core_basis.py', 'rectangular_defect_profiles.py',
                    'rectangular_defect_coordinates.py', 'native_fock_f4.cpp')
    digest = hashlib.sha256()
    for name in dependencies:
        digest.update(name.encode())
        digest.update((directory/name).read_bytes())
    fingerprint = digest.hexdigest()
    cache = None
    if cache_directory:
        cache_directory.mkdir(parents=True, exist_ok=True)
        cache = cache_directory/f'partial-M{M}-a{alpha}-b{beta}-P1-{fingerprint[:16]}.pickle'
        if cache.exists():
            with cache.open('rb') as handle:
                saved = pickle.load(handle)
            if saved['source_sha256'] != fingerprint or saved['parameters'] != [M, alpha, beta]:
                raise ValueError('The local coordinate cache has inconsistent parameters.')
            return saved['data'], saved['compiled'], {'source_sha256': fingerprint, 'reused': True}
    data = [data_for(place, M) for place in ('infinity', 'zero')]
    compiled = []
    for source in data:
        precision = None
        while True:
            result = coordinate_matrix(source, alpha, beta, profile_precision=precision)
            if result[2]['exact_polynomial_matrix']:
                compiled.append(result)
                break
            precision = 2*result[2]['profile_precision']
            if precision > 8192:
                raise ArithmeticError('An exact base numerator needs more than8192 coefficient positions.')
        print(json.dumps({'cache_place': source['place'], 'M': M,
                          'matrix_order': len(result[0]), 'exact': True,
                          'profile_precision': result[2]['profile_precision']}), flush=True)
    if cache:
        temporary = cache.with_suffix('.temporary')
        with temporary.open('wb') as handle:
            pickle.dump({'source_sha256': fingerprint, 'parameters': [M, alpha, beta],
                         'data': data, 'compiled': compiled}, handle, protocol=5)
        temporary.replace(cache)
    return data, compiled, {'source_sha256': fingerprint, 'reused': False}


def evaluate_candidate(M, alpha, beta, candidate, starts, prepared=None):
    begin = time.monotonic()
    spec = None if candidate is None else single_shear(candidate['row'], candidate['column'], candidate['coefficient'])
    if spec is not None and candidate['row'] == 0:
        spec = (*spec, False)
    if prepared is None:
        oracle = CoreDefectPressure(M, alpha, beta, repair_spec=spec)
    else:
        data, baseline = prepared
        compiled = baseline if spec is None else [
            transform_coordinate_matrix(base, source, V=spec[0])
            for base, source in zip(baseline, data)]
        oracle = CoreDefectPressure(M, alpha, beta, compiled=compiled, source_data=data)
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
              'matrix_source': 'cached_physical_core_update' if prepared else 'fresh_CAR_program',
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
    parser.add_argument('--additional-starts', type=Path, action='append', default=[])
    parser.add_argument('--from-screen', type=Path)
    parser.add_argument('--top', type=int, default=8)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--background-only', action='store_true',
                        help='Separately screen the33 additional output shears that change row zero.')
    parser.add_argument('--all-output', action='store_true',
                        help='Screen all264 admissible output shears, including row zero.')
    parser.add_argument('--cache-directory', type=Path)
    parser.add_argument('--rebuild-per-candidate', action='store_true',
                        help='Use fresh CAR matrices instead of the exact coordinate update.')
    args = parser.parse_args()
    begin = time.monotonic()
    starts = starts_for(args.M, args.starts)
    for source in args.additional_starts:
        for point in starts_for(args.M, source):
            if point not in starts:
                starts.append(point)
    if args.background_only and args.all_output:
        raise ValueError('Choose either the background extension or all output shears.')
    if args.from_screen:
        previous = json.loads(args.from_screen.read_text())
        candidates = [entry['candidate'] for entry in sorted(
            (entry for entry in previous['candidates'] if entry['best'] is not None),
            key=lambda entry: entry['best']['rho'])[:args.top]]
    else:
        candidates = [{'row': row, 'column': column, 'coefficient': coefficient,
                       'height_relation': 'equal' if height(row) == height(column) else 'strict'}
                      for row in ((0,) if args.background_only else range(0 if args.all_output else 1, 12)) for column in range(12)
                      if row != column and height(row) <= height(column) for coefficient in (1, 2, 3)]
        candidates.sort(key=lambda candidate: (candidate['height_relation'] != 'equal', candidate['row'], candidate['column'], candidate['coefficient']))
    result = {'status': 'finite fixed-gauge floating screen; no exact rate claimed',
              'M': args.M, 'alpha': args.alpha, 'beta': args.beta,
              'background_row_fixed': not (args.background_only or args.all_output),
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
    prepared = None
    if not args.rebuild_per_candidate:
        data, compiled, cache_information = prepare_coordinates(args.M, args.alpha, args.beta, args.cache_directory)
        prepared = data, compiled
        result['coordinate_cache'] = cache_information
    if result['baseline'] is None:
        result['baseline'] = evaluate_candidate(args.M, args.alpha, args.beta, None, starts, prepared)
        save()
        print(json.dumps({'baseline_rho': result['baseline']['best']['rho']}), flush=True)
    done = {(entry['candidate']['row'], entry['candidate']['column'], entry['candidate']['coefficient']) for entry in result['candidates']}
    best = min([entry['best']['rho'] for entry in result['candidates'] if entry['best']] + [result['baseline']['best']['rho']])
    for candidate in candidates:
        key = (candidate['row'], candidate['column'], candidate['coefficient'])
        if key in done:
            continue
        record = evaluate_candidate(args.M, args.alpha, args.beta, candidate, starts, prepared)
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
