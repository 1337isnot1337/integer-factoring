#!/usr/bin/env python3
"""Exact affine spans of energy-aligned sixteen-letter Laurent frames.

For each of five projective rank-one directions, three F4 scalar frames
span T + l1*(f*T1+f**4*T4) + l2*f**2*T2 + l3*f**3*T3.
Each accepted scalar frame is polynomial, so every such affine sum is
polynomial with the original endpoint and leading coefficient. Complete
weighted profiles are recomputed for the sums.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import signal
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from asymmetric_crossing_pressure import AsymmetricPressure
from dense_coupled_reflection_frames import (
    projective_vectors, dot, rank_one_terms, transformed_from_terms,
    filtration_record, digest, save,
)
from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import physical_bit_energies, tensor
from reflection_laurent_physical_frames import bell_normalize, coupled_admissible_poles, scale
from search_dense_coupled_reflection_frames import profile_rho
from search_physical_code_gauges import start_params
from search_reflection16_coupled_frames import actual_to_scaled, scaled_to_actual


class Deadline(Exception):
    pass


def alarm(signum, frame):
    raise Deadline()


def combine(coefficients, basis):
    return [dot(coefficients, [vector[i] for vector in basis]) for i in range(len(basis[0]))]


def aligned_directions():
    """The two explicit opposite-filtration intersections at energy32."""
    left_basis = [[int(mask >> i & 1) for i in range(16)] for mask in (576, 32872)]
    right_basis = [[int(mask >> i & 1) for i in range(16)] for mask in (576, 5185)]
    result = []
    for coefficients in projective_vectors(2):
        left = combine(coefficients, left_basis)
        rights = [combine(pair, right_basis) for pair in product(range(4), repeat=2) if any(pair)]
        right = next(vector for vector in rights if dot(left, vector) == 0)
        result.append((left, right))
    assert len(result) == 5
    return result


def scaled_terms(terms, scalar):
    result = []
    power = 1
    for entries in terms:
        result.append({index: scale(value, power) for index, value in entries.items()})
        power = F4_PRODUCT[power][scalar]
    return result


def affine_terms(terms, parameters):
    a, b, c = parameters
    return [{index: scale(value, scalar) for index, value in entries.items()
             if scalar} for scalar, entries in zip((1, a, b, c, a), terms)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=2700)
    parser.add_argument('--refine', type=int, default=24)
    parser.add_argument('--trial-seconds', type=float, default=60)
    parser.add_argument('--iterations', type=int, default=100)
    parser.add_argument('--poles', default='1,2,4,8,16,32,64')
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, alarm)
    base = tensor(8)
    _, factors = bell_normalize(base)
    weights = physical_bit_energies(8)
    energy = np.asarray([sum(value for i, value in enumerate(weights) if state >> i & 1)
                         for state in range(16)], float)
    natural = lambda scale: start_params(energy * scale, scale)
    unit = natural(1.0)
    oracle = RectangularOracle(base['values'], 16, 16)
    baseline_orders = oracle.evaluate(unit, [1 / 16] * 16, [1 / 16] * 16, want_spectra=True)[2]
    fitted = minimize_scalar(lambda scale: profile_rho(baseline_orders, energy, scale),
                             bounds=(.005, .4), method='bounded', options={'xatol': 1e-10})
    baseline_scale, baseline_rho = float(fitted.x), float(fitted.fun)
    report = {'scope': 'exact affine strict polynomial paths; bounded numerical complete-profile search only',
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'requested_seconds': args.seconds, 'requested_refinements': args.refine,
              'trial_seconds': args.trial_seconds, 'iterations': args.iterations,
              'alphabet': 16, 'degree': base['degree'],
              'baseline_rho': baseline_rho, 'baseline_scale': baseline_scale,
              'requested_poles': [int(value) for value in args.poles.split(',')],
              'directions': [], 'profiles': [], 'screen': [], 'refinements': [],
              'counters': {'proposed_sums': 0, 'duplicate_tensors': 0, 'unresolved_profiles': 0,
                           'resolved_objectives': 0, 'unresolved_objectives': 0}}
    names = (Path(__file__).name, 'dense_coupled_reflection_frames.py',
             'reflection_laurent_physical_frames.py', 'reflection_dual_rail_gaussian.py',
             'search_dense_coupled_reflection_frames.py', 'asymmetric_crossing_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'search_reflection16_coupled_frames.py')
    report['source_sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in names}
    started = time.monotonic()
    term_cache = []
    for left, right in aligned_directions():
        terms = rank_one_terms(base, left, right, factors)
        assert any(terms[order] for order in range(1, 5))
        filtration = filtration_record(left, right, factors, energy.tolist())
        assert filtration['infinity_active_minimum'] >= 0
        scalar_poles = [coupled_admissible_poles(scaled_terms(terms, scalar), base['degree'])
                        for scalar in (1, 2, 3)]
        common = sorted(set(scalar_poles[0]).intersection(*map(set, scalar_poles[1:])))
        report['directions'].append({'u': left, 'v': right, 'filtration': filtration,
                                     'allowed_poles_by_nonzero_scalar': scalar_poles,
                                     'common_allowed_poles': common})
        term_cache.append(terms)
    seen, profiles = set(), {}
    complete = True
    for direction, terms in enumerate(term_cache):
        for pole in report['requested_poles']:
            if pole not in report['directions'][direction]['common_allowed_poles']:
                continue
            for parameters in product(range(4), repeat=3):
                if not any(parameters):
                    continue
                if time.monotonic() - started > args.seconds * .45:
                    complete = False
                    break
                report['counters']['proposed_sums'] += 1
                changed = transformed_from_terms(base, affine_terms(terms, parameters), pole)
                seed_digest = digest(changed['values'])
                if seed_digest in seen:
                    report['counters']['duplicate_tensors'] += 1
                    continue
                seen.add(seed_digest)
                record = {'direction': direction, 'pole': pole, 'parameters': list(parameters),
                          'seed_sha256': seed_digest, 'degree': changed['degree'],
                          'monomial_content': changed['monomial_content']}
                oracle = RectangularOracle(changed['values'], 16, 16)
                try:
                    raw = oracle.evaluate(unit, [1 / 16] * 16, [1 / 16] * 16, want_spectra=True)[2]
                    key = tuple(round(value) for value in raw)
                    assert max(abs(a - b) for a, b in zip(raw, key)) < 1e-8
                    if key not in profiles:
                        fitted = minimize_scalar(lambda scale: profile_rho(key, energy, scale),
                                                 bounds=(.005, .4), method='bounded',
                                                 options={'xatol': 1e-7, 'maxiter': 35})
                        profiles[key] = len(report['profiles'])
                        report['profiles'].append({'orders': list(key), 'rho': float(fitted.fun),
                                                   'scale': float(fitted.x)})
                    record['profile_index'] = profiles[key]
                    record['natural_rho'] = report['profiles'][record['profile_index']]['rho']
                except (ArithmeticError, ValueError):
                    report['counters']['unresolved_profiles'] += 1
                    record['status'] = 'native profile unresolved; no rank or pressure assigned'
                report['screen'].append(record)
                if len(report['screen']) % 30 == 0:
                    report.update(elapsed_seconds=time.monotonic() - started,
                                  best_natural_rho=min((r['rho'] for r in report['profiles']), default=None))
                    save(args.output, report)
                    print(json.dumps({'screen': len(report['screen']), 'profiles': len(profiles),
                                      'best_rho': report['best_natural_rho'],
                                      'elapsed_seconds': report['elapsed_seconds']}), flush=True)
            if not complete:
                break
        if not complete:
            break
    report.update(complete_requested_screen=complete, screen_seconds=time.monotonic() - started)
    ranked = sorted((r for r in report['screen'] if 'natural_rho' in r), key=lambda r: r['natural_rho'])
    selected, used = [], set()
    for record in ranked:
        if record['profile_index'] not in used:
            used.add(record['profile_index'])
            selected.append(record)
        if len(selected) >= args.refine:
            break
    if args.refine == 0:
        selected = []
    best_support = math.inf
    for case_index, candidate in enumerate(selected):
        remaining = args.seconds - (time.monotonic() - started)
        if remaining <= 0:
            break
        terms = term_cache[candidate['direction']]
        changed = transformed_from_terms(base, affine_terms(terms, candidate['parameters']), candidate['pole'])
        oracle = AsymmetricPressure(changed, 1)
        scale_value = report['profiles'][candidate['profile_index']]['scale']
        first = np.r_[natural(scale_value) / math.sqrt(2), math.pi / 4]
        mode = 'balanced' if case_index % 2 == 0 else 'asymmetric'
        if mode == 'asymmetric':
            first[-1] = .35 if case_index % 4 == 1 else 1.22
        initial = actual_to_scaled(first)
        valid = [math.inf, None, None]

        def objective(variables):
            actual = scaled_to_actual(variables)
            try:
                value, derivative, details = oracle.evaluate(actual)
            except (ArithmeticError, ValueError):
                report['counters']['unresolved_objectives'] += 1
                return 1e6, np.zeros(len(variables))
            report['counters']['resolved_objectives'] += 1
            if value < valid[0]:
                valid[:] = [value, np.asarray(variables).copy(), details]
            gradient = derivative.copy()
            gradient[:-2] *= variables[-2]
            gradient[-2] += float(np.asarray(variables[:-2]) @ derivative[:-2])
            return value, gradient

        bounds = [(-64, 64)] * (len(initial) - 2) + [(.001, .8)]
        bounds.append((math.pi / 4, math.pi / 4) if mode == 'balanced' else (.025, math.pi / 2 - .025))
        begin = time.monotonic()
        record = dict(candidate, mode=mode)
        signal.setitimer(signal.ITIMER_REAL, min(args.trial_seconds, remaining))
        try:
            objective(initial)
            fit = minimize(objective, initial, jac=True, method='SLSQP', bounds=bounds,
                           options={'maxiter': args.iterations, 'ftol': 2e-10})
            record['optimizer'] = {'success': bool(fit.success), 'iterations': int(fit.nit),
                                   'evaluations': int(fit.nfev), 'message': str(fit.message)}
        except Deadline:
            record['optimizer'] = {'message': 'trial deadline; best resolved call retained'}
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        if valid[1] is not None:
            record.update(best_support=valid[0], best_parameters=scaled_to_actual(valid[1]).tolist(),
                          best_scaled_parameters=valid[1].tolist(), details=valid[2])
            best_support = min(best_support, valid[0])
        record['elapsed_seconds'] = time.monotonic() - begin
        report['refinements'].append(record)
        report.update(elapsed_seconds=time.monotonic() - started, best_support_value=best_support)
        save(args.output, report)
        print(json.dumps({'refinements': len(report['refinements']), 'best_support': best_support,
                          'elapsed_seconds': report['elapsed_seconds']}), flush=True)
        if best_support < 1 - 1e-6:
            report['stopping_reason'] = 'numerical threshold candidate needs independent exact replay'
            break
    report.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat(),
                  elapsed_seconds=time.monotonic() - started,
                  stopping_reason=report.get('stopping_reason', 'declared screen/refinement budget reached'))
    save(args.output, report)
    print(json.dumps({key: report[key] for key in ('completed', 'elapsed_seconds', 'stopping_reason')}))


if __name__ == '__main__':
    main()
