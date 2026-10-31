#!/usr/bin/env python3
"""Cached finite screen of opposite-charge endpoint-vanishing cubic pairs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from check_rectangular_core_monomials import positive_monomials, occupation_charge, monomial_program_terms
from rectangular_coordinate_cache import cached_coordinates, extra_monomials_from_bases
from search_rectangular_core_pressure import CoreDefectPressure


def load_starts(M, paths):
    points = []
    for path in paths:
        source = json.loads(path.read_text())
        if isinstance(source, list):
            values, old_M = source, M
        else:
            values, old_M = [source['best']['parameters']], source.get('M', M)
        for value in values:
            point = list(value)
            if len(point) != 23:
                raise ValueError('Supply full-alphabet seven/nine-mode starts with 23 coordinates.')
            point[-1] *= old_M/M
            points.append(point)
    if not points:
        points = [[0.]*22+[tau/M] for tau in (1., 2., 4.)]
    unique = []
    for point in points:
        if point not in unique:
            unique.append(point)
    return unique


def full_parameters(point):
    return point[:13]+[0.]*3+point[13:22]+[0.]*3+point[22:]


def evaluate(M, gamma, full, data, compiled, starts):
    started = time.monotonic()
    oracle = CoreDefectPressure(M, full_core=full, compiled=compiled,
                               source_data=data, full_alphabet_gamma=gamma)
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
    parser.add_argument('--M', type=int, default=4)
    parser.add_argument('--gamma', type=int, choices=(0, 2, 3), default=2)
    parser.add_argument('--starts', type=Path, action='append', default=[])
    parser.add_argument('--minus', type=int, nargs='+', default=[0, 1])
    parser.add_argument('--plus', type=int, nargs='+', default=list(range(8)))
    parser.add_argument('--ells', type=int, nargs='+')
    parser.add_argument('--coefficients', type=int, nargs='+', default=[1, 2, 3])
    parser.add_argument('--plus-pole', type=int, choices=(0, 1), default=0)
    parser.add_argument('--cache-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--stop-file', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    points = load_starts(args.M, args.starts)
    patterns = positive_monomials(args.M-1, degree=3)
    negative = [pattern for pattern in patterns if occupation_charge(pattern) == -1]
    positive = [pattern for pattern in patterns if occupation_charge(pattern) == 1]
    if any(i not in range(2) for i in args.minus) or any(i not in range(8) for i in args.plus):
        raise ValueError('Choose negative labels 0..1 and positive labels 0..7.')
    if any(coefficient not in (1, 2, 3) for coefficient in args.coefficients):
        raise ValueError('The normalized coefficient product must be nonzero in F4.')
    settings = []
    for p in args.plus:
        for m in args.minus:
            extras = [(*negative[m], [(1, 0)], 0), (*positive[p], [(1, 0)], args.plus_pole)]
            # Both places are validated, including the strict endpoint-pole test.
            inf = monomial_program_terms(args.M, 'infinity', extras)
            zero = monomial_program_terms(args.M, 'zero', extras)
            center = (sum(item[1]+item[3] for item in inf)-sum(item[1] for item in zero))//2
            ells = args.ells if args.ells is not None else [0, center, center-args.M//2, center+args.M//2]
            settings.append({'minus': m, 'plus': p, 'ells': list(dict.fromkeys(ells))})
    result = {'status': 'finite fixed-gauge cubic-pair floating screen; no exact pressure certificate',
              'M': args.M, 'full_alphabet_gamma': args.gamma, 'starts': points,
              'negative_patterns': negative, 'positive_patterns': positive,
              'settings': settings, 'coefficient_products': args.coefficients,
              'plus_pole': args.plus_pole,
              'planned_candidates': sum(len(setting['ells'])*len(args.coefficients) for setting in settings),
              'baselines': {}, 'candidates': []}
    prior_seconds = 0
    if args.resume and args.output.exists():
        old = json.loads(args.output.read_text())
        for key in ('M', 'full_alphabet_gamma', 'starts', 'settings', 'coefficient_products', 'plus_pole'):
            if old[key] != result[key]:
                raise ValueError(f'Resume parameter differs: {key}.')
        result = old
        prior_seconds = old.get('elapsed_seconds', 0)
    def save():
        result['elapsed_seconds'] = prior_seconds+time.monotonic()-started
        temporary = args.output.with_suffix(args.output.suffix+'.tmp')
        temporary.write_text(json.dumps(result, indent=2)+'\n')
        temporary.replace(args.output)
    done = {(item['minus'], item['plus'], item['coefficient'], item['ell']) for item in result['candidates']}
    prepared, basis = {}, {}
    for setting in settings:
        m, p, ells = setting['minus'], setting['plus'], setting['ells']
        if all((m, p, coefficient, ell) in done for coefficient in args.coefficients for ell in ells):
            continue
        full = positive[p][0] == args.M-1
        name = 'full' if full else 'partial'
        starts = [full_parameters(point) for point in points] if full else points
        if name not in prepared:
            data, base, cache = cached_coordinates(args.cache_directory, args.M, full_core=full,
                                                    full_alphabet_gamma=args.gamma)
            prepared[name] = data, base
            if name not in result['baselines']:
                result['baselines'][name] = evaluate(args.M, args.gamma, full, data, base, starts)
                result['baselines'][name]['cache'] = cache
                save()
                print(json.dumps({'basis': name, 'baseline_rho': result['baselines'][name]['best']['rho']}), flush=True)
        data, base = prepared[name]
        variations = []
        for side, label, pattern, pole in (('minus', m, negative[m], 0), ('plus', p, positive[p], args.plus_pole)):
            key = name, side, label, pole
            if key not in basis:
                _, unit, cache = cached_coordinates(args.cache_directory, args.M, full_core=full,
                                                     full_alphabet_gamma=args.gamma,
                                                     extra_monomials=[(*pattern, [(1, 0)], pole)])
                basis[key] = unit, cache['digest']
                print(json.dumps({'basis_ready': key, 'matrix_orders': [len(item[0]) for item in unit]}), flush=True)
            variations.append(basis[key])
        best = min([item['best']['rho'] for item in result['baselines'].values()]+[
            item['best']['rho'] for item in result['candidates'] if item['best']])
        for ell in ells:
            for coefficient in args.coefficients:
                if (m, p, coefficient, ell) in done:
                    continue
                if args.stop_file and args.stop_file.exists():
                    result['stopped_by_file'] = True
                    save()
                    return
                compiled = extra_monomials_from_bases(base, [(variations[0][0], [(1, 0)]),
                                                             (variations[1][0], [(coefficient, ell)])])
                record = evaluate(args.M, args.gamma, full, data, compiled, starts)
                record.update({'minus': m, 'plus': p, 'coefficient': coefficient, 'ell': ell,
                               'full_core': full, 'basis_cache': [item[1] for item in variations]})
                result['candidates'].append(record)
                save()
                if record['best']:
                    best = min(best, record['best']['rho'])
                print(json.dumps({'completed': len(result['candidates']), 'best_rho': best,
                                  'candidate_rho': record['best']['rho'] if record['best'] else None,
                                  'candidate': [m, p, coefficient, ell],
                                  'elapsed_seconds': result['elapsed_seconds']}), flush=True)
    result['completed'] = len(result['candidates']) == result['planned_candidates']
    save()


if __name__ == '__main__':
    main()
