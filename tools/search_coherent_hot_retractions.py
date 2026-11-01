#!/usr/bin/env python3
"""Complete q32 binary boundary catalogue or bounded q96 coherent walks.

Every trial retains both complete Smith profiles. Numerical pressure roots
guide the search; they are not exact bounds or a factoring improvement.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from itertools import product
import json
from math import exp, log
from pathlib import Path
from random import Random
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

from coherent_hot_compounds import height, profiles, support


def pressure(profile):
    zero = np.asarray(profile['zero_orders_t'], dtype=float)
    infinity = np.asarray(profile['infinity_orders_t'], dtype=float)
    heights = np.asarray(profile['heights'], dtype=float)

    def objective(log_scale):
        scale = exp(log_scale)
        target = float(logsumexp(-scale * zero) + logsumexp(-scale * infinity))
        equation = lambda rho: target - 4 * rho * float(logsumexp(-scale * heights / rho))
        return brentq(equation, 1e-7, 2.0, xtol=1e-13)

    fit = minimize_scalar(objective, bounds=(log(.01), log(2)), method='bounded',
                          options={'xatol': 1e-10, 'maxiter': 80})
    return {'rho': float(fit.fun), 'scale_t': exp(float(fit.x)),
            'completed': bool(fit.success), 'queries': int(fit.nfev)}


def _lift(vector, words):
    return sum(1 << word for index, word in enumerate(words) if vector >> index & 1)


def q32_catalogue():
    """All rank-three idempotents on the four two-particle height-16 words."""
    lower = sorted((w for w in range(256) if height(w) < 16), key=lambda w: (height(w), w))
    boundary = (24, 36, 66, 129)
    assert len(lower) == 28
    base = [1 << w for w in lower + [15]]
    for kernel in range(1, 16):
        for covector in range(1, 16):
            if (kernel & covector).bit_count() % 2 != 1:
                continue
            pivot = (covector & -covector).bit_length() - 1
            keep = [i for i in range(4) if i != pivot]
            right = [(1 << i) ^ (((covector >> i) & 1) << pivot) for i in keep]
            left = [(1 << i) ^ (covector if kernel >> i & 1 else 0) for i in keep]
            assert all((a & b).bit_count() % 2 == int(i == j)
                       for i, a in enumerate(left) for j, b in enumerate(right))
            yield {'kernel': kernel, 'covector': covector, 'right': right, 'left': left,
                   'boundary_words': list(boundary), 'fixed_boundary_word': 15}, \
                base + [_lift(v, boundary) for v in right], base + [_lift(v, boundary) for v in left], 8
    # The other occupation-preserving rank-four height-16 choice retains all
    # four two-particle states and omits the four-particle state.
    vectors = [1 << w for w in lower + list(boundary)]
    yield {'all_two_particle_boundary': True, 'boundary_words': list(boundary)}, vectors, vectors.copy(), 8


def q96_start():
    words = sorted(range(1 << 12), key=lambda w: (height(w), w))
    lower = [w for w in words if height(w) < 25]
    boundary = [w for w in words if height(w) == 25 and w.bit_count() == 3]
    assert len(lower) == 92 and len(boundary) == 10
    base = [1 << w for w in lower + [31]]
    return base, boundary, [1, 2, 4], [1, 2, 4]


def mutate_projector(right, left, dimension, random):
    right, left = list(right), list(left)
    slot = random.randrange(len(right))
    for _ in range(100):
        vector = random.randrange(1, 1 << dimension)
        change_right = bool(random.randrange(2))
        if change_right:
            projected = vector
            for p, l in zip(right, left):
                if (l & vector).bit_count() % 2:
                    projected ^= p
            if projected:
                right[slot] ^= projected
                break
        else:
            projected = vector
            for p, l in zip(right, left):
                if (p & vector).bit_count() % 2:
                    projected ^= l
            if projected:
                left[slot] ^= projected
                break
    else:
        raise ArithmeticError('No nonzero projector mutation was found.')
    assert all((a & b).bit_count() % 2 == int(i == j)
               for i, a in enumerate(left) for j, b in enumerate(right))
    return right, left


def _atomic_write(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--family', choices=('q32', 'q96'), default='q32')
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--seed', type=int, default=731026)
    parser.add_argument('--max-trials', type=int, default=1000)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.seconds <= 0 or args.max_trials <= 0:
        raise ValueError('Use a positive time and trial budget.')
    names = (Path(__file__).name, 'coherent_hot_compounds.py', 'coherent_hot_compounds.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'family': args.family, 'seed': args.seed,
              'seconds_budget': args.seconds, 'trial_budget': args.max_trials,
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'records': [], 'completed': False,
              'unresolved': 0, 'in_progress': None}
    began = time.monotonic()
    random = Random(args.seed)
    queue = iter(q32_catalogue()) if args.family == 'q32' else None
    if args.family == 'q96':
        base, boundary, current_right, current_left = q96_start()
        current_rho = float('inf')
    best = None
    exhausted = False
    while len(report['records']) < args.max_trials and time.monotonic() - began < args.seconds:
        index = len(report['records'])
        if queue is not None:
            try:
                description, embedding, retraction, modes = next(queue)
            except StopIteration:
                exhausted = True
                break
        else:
            right, left = ((current_right.copy(), current_left.copy()) if index == 0 else
                           mutate_projector(current_right, current_left, len(boundary), random))
            description = {'right': right, 'left': left, 'boundary_words': boundary,
                           'fixed_boundary_word': 31}
            embedding = base + [_lift(v, boundary) for v in right]
            retraction = base + [_lift(v, boundary) for v in left]
            modes = 12
        report['in_progress'] = {'index': index, 'description': description,
                                 'embedding': embedding, 'retraction': retraction, 'modes': modes}
        _atomic_write(args.output, report)
        start = time.monotonic()
        record = {'index': index, 'description': description, 'completed': False}
        try:
            result = profiles(embedding, retraction, modes)
            fit = pressure(result)
            record.update(completed=True, embedding=embedding, retraction=retraction, modes=modes,
                          heights=result['heights'], charges=result['charges'], pressure=fit,
                          zero_orders_t=result['zero_orders_t'], infinity_orders_t=result['infinity_orders_t'],
                          block_sizes=[b['size'] for b in result['cuts'][0]['blocks']],
                          precisions=[[b['precision_t'] for b in cut['blocks']] for cut in result['cuts']])
            if best is None or fit['rho'] < best['pressure']['rho']:
                best = record
                report['best_record_index'] = index
                print(json.dumps({'family': args.family, 'record': index, 'rho': fit['rho'],
                                  'elapsed': time.monotonic() - began}), flush=True)
            if queue is None:
                temperature = 1e-6
                if fit['rho'] <= current_rho or random.random() < exp(min(0., (current_rho - fit['rho']) / temperature)):
                    current_right, current_left, current_rho = right, left, fit['rho']
        except (ValueError, ArithmeticError) as error:
            record['error'] = f'{type(error).__name__}: {error}'
            report['unresolved'] += 1
        record['seconds'] = time.monotonic() - start
        report['records'].append(record)
        report['in_progress'] = None
        report['elapsed_seconds'] = time.monotonic() - began
        _atomic_write(args.output, report)
    report.update(completed=True, finite_list_exhausted=exhausted,
                  elapsed_seconds=time.monotonic() - began,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    _atomic_write(args.output, report)


if __name__ == '__main__':
    main()
