#!/usr/bin/env python3
"""Independent minors, polynomial frames, and saved bounded-type gauge replay.

The natural integer-height profiles have exact polynomial-frame proofs.
Optional full-gauge states are independently reassembled, but their real
height comparisons and entropy/KKT checks remain floating point.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

from balanced_particle_extensions import construct, coordinates, extension_cases, matrices
from check_aligned_gaussian_complements import expand_square, minor, multiply, strict_matrix
from check_f4_recentered_walks import certify_profile, same_cost
from f4_recentered_lattices import pack
from native_matrix_smith_f4 import MatrixOracle
from replay_source_provenance import resolve_source_hashes
from search_f4_recentered_crossings import digest, root_for_profiles


ROOT = Path(__file__).resolve().parent.parent


def certify_tensor(seed):
    a, b = seed['dimensions'][:2]
    assert seed['dimensions'] == [a, b, b, a]
    for key, size in (('code_x', a), ('code_y', b)):
        code = seed[key]
        assert len(code) == len(set(code)) == size
        assert all(0 <= mask < 256 and mask.bit_count() == 4 for mask in code)
    assert len(seed['values']) == a*a*b*b
    common = seed['entry_gcd_in_z']
    assert common > 0 and common.bit_count() % 2 == 1
    matrix = strict_matrix()
    for index, (value, high) in enumerate(seed['values']):
        x, y, u, v = coordinates(index, a, b)
        assert high == 0 and value.bit_count() % 2 == int(x == v and y == u)
        actual = minor(matrix, seed['code_x'][v] | (seed['code_y'][y] << 8),
                       seed['code_y'][u] | (seed['code_x'][x] << 8))
        assert actual == expand_square(multiply(common, value))
    assert max(value.bit_length()-1 for value, _ in seed['values']) == seed['degree']
    for key, energy in (('code_x', 'energies_x'), ('code_y', 'energies_y')):
        grades = [sum(2*i+1 for i in range(8) if (mask >> i) & 1) for mask in seed[key]]
        assert list(seed[energy]) == [value-min(grades) for value in grades]
    return a*a*b*b


def local_integer_matrix(seed, side):
    a, b = seed['dimensions'][:2]
    n, degree = a*b, seed['degree']
    hx = [value//2 for value in seed['energies_x']]
    hy = [value//2 for value in seed['energies_y']]
    result = [[0]*n for _ in range(n)]
    for index, pair in enumerate(seed['values']):
        x, y, u, v = coordinates(index, a, b)
        if side == 0:
            pair = tuple(int(format(value, '0'+str(degree+1)+'b')[::-1], 2) for value in pair)
            row, column, weight = b*x+u, a*y+v, hy[u]+hx[v]
        else:
            row, column, weight = b*x+y, a*u+v, hx[x]+hy[y]
        result[row][column] = pack(pair) << (2*weight)
    return result, degree if side == 0 else 0


def natural_replay(path):
    record = json.loads(path.read_text())
    provenance, unresolved = resolve_source_hashes(record['source_sha256'])
    assert not unresolved, unresolved
    expected = list(extension_cases())
    rows = record['results']
    assert len(rows) == record['complete_tensor_queries'] == 196
    assert record['completed_requested_schedule'] and record['completed_full_extension_schedule']
    minors = profiles = 0
    for index, (row, (family, added, xcode, ycode)) in enumerate(zip(rows, expected)):
        assert (row['family'], row['added_mask']) == (family, added)
        assert tuple(row['code_x']) == xcode and tuple(row['code_y']) == ycode
        seed = construct(xcode, ycode)
        assert digest(seed['values']) == row['digest']
        assert seed['entry_gcd_in_z'] == row['entry_gcd_in_z']
        assert seed['degree'] == row['degree']
        minors += certify_tensor(seed)
        for side in (0, 1):
            matrix, scalar = local_integer_matrix(seed, side)
            certify_profile(matrix, scalar, row['profiles'][side])
            profiles += 1
        same_cost(root_for_profiles(row['profiles'], seed['energies_x'], seed['energies_y']), row)
        if (index+1) % 20 == 0:
            print(json.dumps({'event': 'natural_replay', 'tensors': index+1,
                              'literal_minors': minors, 'exact_profiles': profiles}), flush=True)
    assert minors == record['endpoint_entries_checked']
    assert abs(min(row['rho'] for row in rows)-record['best_recorded_rho']) < 1e-10
    return {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(),
            'source_provenance': provenance, 'complete_extension_tensors': len(rows),
            'literal_complete_gaussian_minors': minors,
            'complete_integer_profiles_certified_by_regular_frames': profiles,
            'best_recorded_rho': record['best_recorded_rho']}


def close(actual, expected, tolerance=3e-7):
    if isinstance(actual, (tuple, list)):
        assert len(actual) == len(expected)
        for first, second in zip(actual, expected):
            close(first, second, tolerance)
    else:
        assert abs(actual-expected) <= tolerance*max(1., abs(actual), abs(expected)), (actual, expected)


def type_kkt(energy, probabilities, lower, rate):
    assert len(energy) == len(probabilities) == len(lower)
    assert abs(sum(probabilities)-1) < 1e-10
    assert all(p >= bound-1e-12 and p >= 0 for p, bound in zip(probabilities, lower))
    potentials = [g+rate*math.log(p) if p > 0 else -math.inf for g, p in zip(energy, probabilities)]
    free = [i for i, (p, bound) in enumerate(zip(probabilities, lower)) if p > bound+1e-11]
    assert free
    level = potentials[free[0]]
    for i in free:
        close(potentials[i], level)
    # Coordinates whose mass is tiny but above a zero bound are also free;
    # a zero numerical mass is only acceptable as harmless exponential underflow.
    for i, (p, bound) in enumerate(zip(probabilities, lower)):
        if i in free or bound == 0:
            if p > 0 and bound == 0:
                close(potentials[i], level)
            continue
        assert potentials[i] >= level-3e-7*max(1., abs(level))


def gauge_state(seed, state, lower_x, lower_y):
    from scipy.special import logsumexp
    a, b = seed['dimensions'][:2]
    alphabets = (a, b, b, a)
    params = state['parameters']
    rows, offset = [], 0
    for size in alphabets*2:
        rows.append((0.,)+tuple(params[offset:offset+size-1]))
        offset += size-1
    assert len(params) == offset+1 and params[-1] > 0
    first, second = matrices(seed)
    p = MatrixOracle(first, a*b).profile(
        row_weights=[rows[0][x]+rows[2][u] for x in range(a) for u in range(b)],
        column_weights=[rows[1][y]+rows[3][v] for y in range(b) for v in range(a)],
        scale=params[-1], infinity=True)['orders']
    q = MatrixOracle(second, a*b).profile(
        row_weights=[rows[4][x]+rows[5][y] for x in range(a) for y in range(b)],
        column_weights=[rows[6][u]+rows[7][v] for u in range(b) for v in range(a)],
        scale=params[-1])['orders']
    gx = [sum(rows[leg][i] for leg in (0, 3, 4, 7)) for i in range(a)]
    gy = [sum(rows[leg][i] for leg in (1, 2, 5, 6)) for i in range(b)]
    close(p+q, state['orders_times_log_s'])
    close(gx, state['energy_x_times_log_s'])
    close(gy, state['energy_y_times_log_s'])
    px, py = state['type_x'], state['type_y']
    rate = 2*state['rho']
    type_kkt(gx, px, lower_x, rate)
    type_kkt(gy, py, lower_y, rate)
    entropy = -sum(value*math.log(value) for value in px+py if value > 0)
    numerator = float(logsumexp([-value for value in p])+logsumexp([-value for value in q]))
    numerator += sum(value*energy for value, energy in zip(px, gx))
    numerator += sum(value*energy for value, energy in zip(py, gy))
    close(numerator/(2*entropy), state['rho'])
    return 2


def gauge_replay(path):
    record = json.loads(path.read_text())
    assert record['completed_requested_schedule']
    provenance, unresolved = resolve_source_hashes(record['source_sha256'])
    assert not unresolved, unresolved
    for name, expected in record['input_sha256'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest() == expected
    minors = profiles = states = 0
    for row in record['candidates']:
        seed = row['literal_seed']
        assert digest(seed['values']) == row['source_row']['digest']
        minors += certify_tensor(seed)
        for run in row['runs']:
            for key, code in (('lower_x', seed['code_x']), ('lower_y', seed['code_y'])):
                assert run[key] == [run['minimum_added_mass'] if mask == row['source_row']['added_mask'] else 0.
                                    for mask in code]
            if run['best'] is None:
                continue
            profiles += gauge_state(seed, run['best'], run['lower_x'], run['lower_y'])
            states += 1
        close(min(run['best']['rho'] for run in row['runs'] if run['best'] is not None), row['best']['rho'])
    return {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(),
            'source_provenance': provenance, 'literal_complete_gaussian_minors': minors,
            'saved_floating_gauge_states': states,
            'independent_complete_weighted_matrix_profiles': profiles,
            'best_recorded_rho': record['best_recorded_rho'],
            'scope': 'Floating profile replay, type lower bounds, entropy quotient and sufficient Gibbs KKT conditions; no exact gauge optimum or pressure-sign certificate.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--natural', type=Path)
    parser.add_argument('--gauges', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not args.natural and not args.gauges:
        parser.error('Supply a natural search record, a gauge record, or both.')
    started = time.monotonic()
    report = {'status': 'passed', 'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.natural:
        report['natural'] = natural_replay(args.natural)
    if args.gauges:
        report['gauges'] = gauge_replay(args.gauges)
    report['elapsed_seconds'] = time.monotonic()-started
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
