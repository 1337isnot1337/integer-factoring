#!/usr/bin/env python3
"""Complete profiles of repeated coherent physical cores at several heights.

Each four-mode core is deletion closed. Their physical tensor product gives
relations at several energy ranges. Scale discovery uses a coarse global
grid before refining every detected local minimum, avoiding a single-start
unimodality assumption. Numerical minima are not exact pressure witnesses.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from itertools import product
import json
from math import exp, log
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

from coherent_hot_compounds import support, validate_maps
from coherent_hot_compounds_fast import profiles
from search_coherent_hot_retractions import _atomic_write


def core_maps(size, coherent):
    words = {8: [0, 1, 2, 3, 4, 5, 8, 9],
             10: [0, 1, 2, 4, 8, 3, 5, 10, 12, 9]}[size]
    embedding = [1 << word for word in words]
    retraction = embedding.copy()
    if coherent:
        retraction[words.index(9)] ^= 1 << 6
    validate_maps(embedding, retraction, 4)
    return embedding, retraction


def product_maps(size, copies, coherent):
    if copies < 1 or 4 * copies > 16:
        raise ValueError('Use one through four four-mode cores.')
    embedding, retraction = core_maps(size, coherent)

    def tensor(vectors):
        answer = []
        for labels in product(range(size), repeat=copies):
            value = 0
            for words in product(*(support(vectors[label]) for label in labels)):
                occupied = sum(word << (4 * block) for block, word in enumerate(words))
                value ^= 1 << occupied
            answer.append(value)
        return answer

    result = tensor(embedding), tensor(retraction)
    validate_maps(*result, 4 * copies)
    return result


def pressure_grid(profile):
    zero = np.asarray(profile['zero_orders_t'], dtype=float)
    infinity = np.asarray(profile['infinity_orders_t'], dtype=float)
    heights = np.asarray(profile['heights'], dtype=float)
    evaluations = 0

    def objective(log_scale):
        nonlocal evaluations
        evaluations += 1
        scale = exp(float(log_scale))
        target = float(logsumexp(-scale * zero) + logsumexp(-scale * infinity))
        return brentq(lambda rho: target - 4 * rho * float(logsumexp(-scale * heights / rho)),
                      1e-7, 2.0, xtol=1e-13)

    grid = np.linspace(log(.001), log(4.), 101)
    values = [objective(point) for point in grid]
    candidates = [{'rho': values[i], 'scale_t': exp(float(grid[i])), 'kind': 'grid_endpoint'}
                  for i in (0, len(grid) - 1)]
    for i in range(1, len(grid) - 1):
        if values[i] <= min(values[i - 1], values[i + 1]):
            fit = minimize_scalar(objective, bounds=(grid[i - 1], grid[i + 1]), method='bounded',
                                  options={'xatol': 1e-10, 'maxiter': 100})
            candidates.append({'rho': float(fit.fun), 'scale_t': exp(float(fit.x)),
                               'kind': 'refined_grid_minimum', 'success': bool(fit.success)})
    return {'best': min(candidates, key=lambda item: item['rho']), 'candidates': candidates,
            'grid': [{'scale_t': exp(float(point)), 'rho': value}
                     for point, value in zip(grid, values)], 'evaluations': evaluations,
            'global_optimality_claimed': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes', default='8,10')
    parser.add_argument('--copies', type=int, default=2)
    parser.add_argument('--coordinate-controls', action='store_true')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    sizes = list(map(int, args.sizes.split(',')))
    if any(size not in (8, 10) for size in sizes):
        raise ValueError('Supported core sizes are eight and ten.')
    names = (Path(__file__).name, 'coherent_hot_compounds.py', 'coherent_hot_compounds_fast.py',
             'coherent_hot_compounds_fast.cpp', 'search_coherent_hot_retractions.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'completed': False, 'sizes': sizes, 'copies': args.copies,
              'coordinate_controls': args.coordinate_controls, 'records': [],
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started = time.monotonic()
    for size in sizes:
        for coherent in ((False, True) if args.coordinate_controls else (True,)):
            p, l = product_maps(size, args.copies, coherent)
            description = {'core_size': size, 'copies': args.copies, 'coherent': coherent,
                           'q': len(p), 'embedding': p, 'retraction': l}
            report['in_progress'] = description
            _atomic_write(args.output, report)
            print(json.dumps({key: value for key, value in description.items()
                              if key not in ('embedding', 'retraction')}), flush=True)
            profile = profiles(p, l, 4 * args.copies)
            fit = pressure_grid(profile)
            report['records'].append({'description': description, 'profile': profile,
                                       'pressure': fit})
            report['in_progress'] = None
            report['elapsed_seconds'] = time.monotonic() - started
            _atomic_write(args.output, report)
            print(json.dumps({'q': len(p), 'coherent': coherent, 'best': fit['best'],
                              'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(completed=True, finished_utc=datetime.now(timezone.utc).isoformat())
    _atomic_write(args.output, report)


if __name__ == '__main__':
    main()
