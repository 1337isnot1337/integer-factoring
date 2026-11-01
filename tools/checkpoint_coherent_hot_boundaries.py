#!/usr/bin/env python3
"""Archive the complete boundary-idempotent grid and replay its four winners."""
from __future__ import annotations
import argparse
import gzip
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp
from check_coherent_hot_boundary import full_maps, logical, natural_profiles
from coherent_hot_boundary import compress
from hot_code_charge_profiles import hot_energy
from physical_split_crossings import field_dot, retract
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path('evidence/coherent-hot-boundary-search.json.gz'))
    parser.add_argument('--output', type=Path, default=Path('evidence/coherent-hot-boundary-search-replay.json'))
    parser.add_argument('--build-from', type=Path, nargs=3, metavar=('SMALL', 'PREFIX', 'LAST'))
    args = parser.parse_args()
    if args.build_from:
        runs = []
        for role, path in zip(('small', 'prefix', 'last'), args.build_from):
            raw = path.read_bytes()
            runs.append({'role': role, 'original_path': str(path), 'sha256': sha256(raw).hexdigest(),
                         'bytes': len(raw), 'raw_json': raw.decode()})
        args.archive.write_bytes(gzip.compress(json.dumps({'runs': runs}, separators=(',', ':')).encode(), mtime=0))
    compressed = args.archive.read_bytes()
    envelope = json.loads(gzip.decompress(compressed))
    assert len(envelope['runs']) == 3 and {run['role'] for run in envelope['runs']} == {'small', 'prefix', 'last'}
    reports, bindings = {}, []
    for run in envelope['runs']:
        raw = run['raw_json'].encode()
        assert len(raw) == run['bytes'] and sha256(raw).hexdigest() == run['sha256']
        reports[run['role']] = json.loads(raw)
        bindings.append({key: run[key] for key in ('role', 'sha256', 'bytes')})
        for name, expected in reports[run['role']]['source_sha256'].items():
            assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected, name
    small, prefix, last = (reports[role] for role in ('small', 'prefix', 'last'))
    assert small['completed'] and small['finite_list_exhausted'] and len(small['records']) == 40
    assert not prefix['completed'] and len(prefix['records']) == 355
    assert last['completed'] and len(last['records']) == 1
    assert last['prefix_sha256'] == next(row['sha256'] for row in bindings if row['role'] == 'prefix')
    assert last['prefix_saved_cases'] == 355 and last['combined_saved_cases'] == 356
    assert last['prefix_elapsed_seconds'] == prefix['elapsed_seconds']
    rows = small['records'] + prefix['records'] + last['records']
    labels = sorted(range(256), key=lambda value: (hot_energy(value), value))
    expected, winners = set(), {}
    for q in (8, 10, 12, 16):
        level = hot_energy(labels[q - 1])
        boundary = [x for x in labels if hot_energy(x) == level]
        for right in product(range(4), repeat=len(boundary)):
            if not any(right) or next(value for value in right if value) != 1:
                continue
            for left in product(range(4), repeat=len(boundary)):
                if field_dot(left, right) == 1:
                    expected.add((q, left, right))
    observed = [(r['q'], tuple(r['left']), tuple(r['right'])) for r in rows]
    assert len(observed) == len(set(observed)) == 396 and set(observed) == expected
    maximum_pressure_error = 0.
    for row in rows:
        q, level = row['q'], row['energy']
        assert row['lower_letters'] == [x for x in labels if hot_energy(x) < level]
        assert row['boundary_letters'] == [x for x in labels if hot_energy(x) == level]
        assert row['energies'] == list(map(hot_energy, row['lower_letters'])) + [level]
        sides = row['orders_infinity_zero']
        assert len(sides) == 2 and all(len(side) == q * q for side in sides)
        assert all(side == sorted(side) and min(side) == 0 and all(type(x) is int for x in side) for side in sides)
        scale = row['pressure']['scale']
        logsum = sum(float(logsumexp(-scale * np.asarray(side))) for side in sides)
        energies = np.asarray(row['energies'])
        root = brentq(lambda R: logsum - 2 * R * float(logsumexp(-scale * energies / R)), 1e-10, 4, xtol=1e-13)
        error = abs(root / 2 - row['pressure']['rho'])
        assert error < 2e-12 and row['pressure']['completed']
        maximum_pressure_error = max(maximum_pressure_error, error)
        if q not in winners or row['pressure']['rho'] < winners[q]['pressure']['rho']:
            winners[q] = row
    replayed = []
    coefficient_count = channel_count = 0
    for q, row in sorted(winners.items()):
        lower, boundary = row['lower_letters'], row['boundary_letters']
        raw_seed = restricted_tensor(16, list(map(logical, lower + boundary)))
        left, right = [row['left']], [row['right']]
        seed = compress(raw_seed, len(lower), left, right)
        left_full, right_full = full_maps(len(lower), len(boundary), left, right)
        reference = retract(raw_seed, left_full, right_full)
        assert seed['values'] == reference['values']
        coefficient_count += len(seed['values'])
        content = min((part & -part).bit_length() - 1 for value in seed['values'] for part in value if part)
        seed['values'] = [(a >> content, b >> content) for a, b in seed['values']]
        degree = max(max(a.bit_length(), b.bit_length()) - 1 for a, b in seed['values'])
        assert content == row['monomial_content'] and degree == row['primitive_degree']
        assert digest(seed['values']) == row['tensor_sha256']
        literal = natural_profiles(seed, row['energies'])
        shift = 128 - content
        expected_sides = [[value - shift for value in row['orders_infinity_zero'][0]],
                          [value + shift for value in row['orders_infinity_zero'][1]]]
        assert literal == expected_sides
        channel_count += 2 * q * q
        replayed.append({'q': q, 'left': row['left'], 'right': row['right'],
                         'degree': degree, 'tensor_sha256': row['tensor_sha256'],
                         'natural_profiles_infinity_zero': literal, 'pressure': row['pressure']})
    replay_sources = (Path(__file__).name, 'check_coherent_hot_boundary.py',
                      'reflection_dual_rail_gaussian.py', 'check_schur_gaussian_reduction.py',
                      'selected_gaussian_coefficients.cpp', 'native_weighted_smith_f4.py',
                      'reflection_laurent_physical_frames.py')
    report = {'scope': 'Complete finite boundary-idempotent coverage; all saved pressures replayed and four literal winners independently reconstructed. Numerical minima are not optimality proofs.',
              'archive_sha256': sha256(compressed).hexdigest(), 'archive_bytes': len(compressed),
              'original_byte_bindings': bindings,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in replay_sources},
              'counts': {'idempotents': len(rows), 'saved_pressure_replays': len(rows),
                         'literal_winners': len(replayed), 'independent_coefficients': coefficient_count,
                         'complete_literal_channels': channel_count},
              'coverage': {str(q): sum(row['q'] == q for row in rows) for q in sorted(winners)},
              'interruption': {'prefix_records': 355, 'explicit_continuation_records': 1,
                               'prefix_saved_seconds': prefix['elapsed_seconds'],
                               'continuation_seconds': last['elapsed_seconds'],
                               'scope': 'The original prefix remains marked incomplete. The exact union, not a claimed uninterrupted run, exhausts the deterministic grid.'},
              'small_run': {key: small[key] for key in ('started_utc', 'finished_utc', 'elapsed_seconds')},
              'maximum_saved_pressure_error': maximum_pressure_error,
              'winner_replays': replayed, 'passed': True}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'counts': report['counts'], 'archive_bytes': len(compressed)}))


if __name__ == '__main__':
    main()
