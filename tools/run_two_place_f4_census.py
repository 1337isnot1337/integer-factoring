#!/usr/bin/env python3
"""Registered, resumable exact F4 pencil census with independent class replays.

Every q=(1+t)r with deg(r)<=10 is enumerated. Profile weights are integers.
Supporting-line optimization is numerical and never establishes a bound.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
from hashlib import sha256
import json
import math
from pathlib import Path
import secrets
import subprocess
import time

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp

from check_joint_energy_conjugations import cut_matrix
from weighted_smith_oracle_f4 import WeightedSmithOracle

TOTAL = 1 << 22


def q_planes(code):
    result = []
    for plane in range(2):
        r = sum(((code >> (2*i + plane)) & 1) << i for i in range(11))
        result.append(r ^ (r << 1))
    return result


def read_input(path):
    fields = path.read_text().split()
    if len(fields) != 1 + 4 * 625 + 33 or fields[0] != 'PENCIL1':
        raise ValueError('Unexpected registered pencil input.')
    entries = [tuple(map(int, fields[1 + 4*i:1 + 4*(i+1)])) for i in range(625)]
    parameters = list(map(float, fields[2501:]))
    if any(value != round(value) for value in parameters) or parameters[-1] != 1:
        raise ValueError('The registered profile weights must be integral with scale one.')
    return entries, parameters


def replay_class(record, entries):
    lo, hi = record['q_planes']
    values = [(a ^ (lo if k else 0) ^ (lo << 4 if l else 0),
               b ^ (hi if k else 0) ^ (hi << 4 if l else 0)) for a, b, k, l in entries]
    half = [0, 5, 3, 1, 4]
    weights = [half[i // 5] for i in range(25)]
    infinity, _ = WeightedSmithOracle(cut_matrix(values, True), infinity=True).profile(weights, weights)
    energies = [half[i % 5] + half[i // 5] for i in range(25)]
    zero, _ = WeightedSmithOracle(cut_matrix(values, False)).profile(energies, [0] * 25)
    exact = list(infinity + zero)
    if exact != record['orders']:
        raise ArithmeticError('An independent fraction-free complete profile disagreed.')
    return {'status': 'pass', 'complete_integer_orders': 50,
            'method': 'Packed-F4 fraction-free exact polynomial elimination; no series truncation.'}


def asymmetric_fit(orders):
    profiles = [np.asarray(orders[:25], dtype=float), np.asarray(orders[25:], dtype=float)]
    minima = [min(values) for values in profiles]
    centered = [values - minimum for values, minimum in zip(profiles, minima)]
    energies = np.asarray([[0, 10, 6, 2, 8]] * 2, dtype=float)

    def value(parameters):
        angle, log_scale = parameters
        a, b, t = math.cos(angle), math.sin(angle), math.exp(log_scale)
        phi = (a * logsumexp(-t * centered[0] / a)
               + b * logsumexp(-t * centered[1] / b) - t * sum(minima))
        if phi <= 1e-12:
            return 2.

        def equation(rate):
            return phi - rate * float(logsumexp(-t * energies / rate, axis=1).sum())

        return brentq(equation, 1e-7, 8., xtol=2e-12)

    best = None
    for initial in ((.25, -1.), (math.pi / 4, -1.), (1.3, -1.)):
        fit = minimize(value, initial, method='Nelder-Mead',
                       bounds=((.01, math.pi / 2 - .01), (-8., math.log(4.))),
                       options={'maxiter': 400, 'xatol': 1e-8, 'fatol': 1e-10})
        result = {'R': float(fit.fun), 'angle': float(fit.x[0]),
                  'scale': math.exp(fit.x[1]), 'success': bool(fit.success),
                  'scope': 'Finite floating search; not a global optimum or an exact pressure certificate.'}
        if best is None or result['R'] < best['R']:
            best = result
    return best


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--state', type=Path, required=True)
    parser.add_argument('--chunk', type=int, default=32768)
    parser.add_argument('--stop', type=int, default=TOTAL,
                        help='A bounded prefix checkpoint; full coverage still requires 2^22 positions.')
    args = parser.parse_args()
    if not 0 < args.chunk <= TOTAL or not 0 < args.stop <= TOTAL:
        parser.error('Use positive chunk and stop bounds at most 2^22.')
    args.state.mkdir(parents=True, exist_ok=True)
    lock = (args.state / 'runner.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    root = Path(__file__).resolve().parents[1]
    source_paths = [Path(__file__), root / 'tools/scan_two_place_f4_pencil.cpp',
                    root / 'tools/native_weighted_smith_f4.cpp',
                    root / 'tools/weighted_smith_oracle_f4.py',
                    root / 'tools/check_joint_energy_conjugations.py']
    fingerprints = {str(path.relative_to(root)): sha256(path.read_bytes()).hexdigest()
                    for path in source_paths}
    manifest_path = args.state / 'manifest.json'
    identity = {'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
                'binary_sha256': sha256(args.binary.read_bytes()).hexdigest(),
                'source_sha256': fingerprints}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if any(manifest[key] != value for key, value in identity.items()):
            raise ValueError('The resumed source, binary or mathematical input changed.')
    else:
        seed = secrets.randbits(64)
        selector = __import__('random').Random(seed)
        manifest = {**identity, 'created_utc': datetime.now(timezone.utc).isoformat(),
                    'os_random_seed': seed, 'multiplier': selector.randrange(TOTAL // 2) * 2 + 1,
                    'offset': selector.randrange(TOTAL), 'total_registered': TOTAL,
                    'family': 'q=(1+t)r, r=sum_(j=0)^10 c_j t^j, c_j in F4',
                    'coverage': 'An odd affine permutation modulo 2^22 is bijective; every coefficient vector occurs exactly once.',
                    'construction_source_commit': '8a3ea3e6aeacba7d0a53764de5579e3707bb2ca7',
                    'scope': 'Complete finite F4 family at inherited integer heights. Numerical asymmetric fits are separate; no all-field or factoring-exponent claim.'}
        atomic_json(manifest_path, manifest)
    entries, _ = read_input(args.input)
    input_bytes = args.input.read_bytes()
    multiplier, offset = manifest['multiplier'], manifest['offset']
    inverse = pow(multiplier, -1, TOTAL)
    classes = {}
    completed = unresolved = endpoint_failures = chunks = 0
    start_time = time.monotonic()

    def absorb(chunk):
        nonlocal completed, unresolved, endpoint_failures, chunks
        if chunk['start'] != completed or chunk['stop'] <= completed:
            raise ValueError('Chunk coverage is not contiguous.')
        if chunk['multiplier'] != multiplier or chunk['offset'] != offset:
            raise ValueError('A chunk used a different candidate permutation.')
        if chunk['visited'] != chunk['stop'] - chunk['start']:
            raise ValueError('The chunk visit count disagrees with its interval.')
        if sum(row['count'] for row in chunk['classes']) + chunk['unresolved'] + chunk['endpoint_failures'] != chunk['visited']:
            raise ValueError('The chunk outcome ledger does not cover every visited candidate.')
        for record in chunk['classes']:
            code = record['representative']
            position = (inverse * (code - offset)) % TOTAL
            if not chunk['start'] <= position < chunk['stop'] or q_planes(code) != record['q_planes']:
                raise ValueError('A profile representative is outside its registered chunk.')
            key = tuple(record['orders'])
            if len(key) != 50 or sum(key[25:]) != 444:
                raise ValueError('A complete profile or determinant invariant is missing.')
            if key not in classes:
                replay = replay_class(record, entries)
                fit = asymmetric_fit(key)
                classes[key] = {**record, 'count': 0, 'independent_replay': replay,
                                'asymmetric_fit': fit}
                if fit['R'] < .999999:
                    atomic_json(args.state / 'candidate-alert.json', classes[key])
            classes[key]['count'] += record['count']
        completed = chunk['stop']
        unresolved += chunk['unresolved']
        endpoint_failures += chunk['endpoint_failures']
        chunks += 1

    for path in sorted(args.state.glob('chunk-*.json')):
        absorb(json.loads(path.read_text()))
    while completed < args.stop:
        stop = min(completed + args.chunk, args.stop)
        result = subprocess.run([str(args.binary), str(completed), str(stop),
                                 str(multiplier), str(offset)], input=input_bytes,
                                capture_output=True, check=True)
        chunk = json.loads(result.stdout)
        path = args.state / f'chunk-{completed:07d}-{stop:07d}.json'
        atomic_json(path, chunk)
        absorb(chunk)
        best = min(classes.values(), key=lambda row: row['asymmetric_fit']['R'])
        summary = {'status': 'complete' if completed == TOTAL else 'incomplete',
                   'visited': completed, 'registered': TOTAL, 'chunks': chunks,
                   'unresolved': unresolved, 'endpoint_failures': endpoint_failures,
                   'distinct_complete_profiles': len(classes),
                   'independently_replayed_class_orders': 50 * len(classes),
                   'seconds_this_invocation': time.monotonic() - start_time,
                   'best_numerical_support': best,
                   'classes': list(classes.values()), 'scope': manifest['scope']}
        atomic_json(args.state / 'summary.json', summary)
        print(json.dumps({key: summary[key] for key in
                          ('status', 'visited', 'registered', 'unresolved', 'distinct_complete_profiles',
                           'independently_replayed_class_orders', 'seconds_this_invocation')}
                         | {'best_numerical_R': best['asymmetric_fit']['R']}), flush=True)


if __name__ == '__main__':
    main()
