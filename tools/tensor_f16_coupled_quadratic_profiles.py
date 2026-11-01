#!/usr/bin/env python3
"""F16 coupled height-eight/height-twelve screen over 851 Frobenius classes.

These classes cover all 3375 normalized nonzero triples. The initial queue
uses the binary basis 1, omega, z, omega*z. Complete profile multiplicities are
retained; orthogonal infinity profiles are explicitly marked as derived.
Numerical rho values are discovery data, not factoring certificates.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

import tensor_hot_compounds_f16 as field
from search_coherent_hot_retractions import pressure


def canonical_frobenius(coefficients):
    orbit, current = [], tuple(coefficients)
    for _ in range(4):
        orbit.append(current)
        current = tuple(field.mul(x, x) for x in current)
    assert current == tuple(coefficients)
    return min(orbit)


def orthogonal_profile(code, left):
    p, modes = code['embedding'], code['modes']
    heights, charges = field.validate_maps(p, left, modes)
    q, orders = len(p), []
    started = time.monotonic()
    for charge in sorted({a + b for a in charges for b in charges}):
        pairs = [(i, j) for i in range(q) for j in range(q) if charges[i] + charges[j] == charge]
        rows = [field.tensor_vector(left[i], left[j], modes) for i, j in pairs]
        columns = [field.tensor_vector(p[j], p[i], modes) for i, j in pairs]
        block = field.compound_block(modes, rows, columns)
        orders.extend(block['orders_t'])
    orders.sort()
    assert len(orders) == q * q
    return {'q': q, 'modes': modes, 'heights': heights, 'charges': charges,
            'zero_orders_t': orders, 'infinity_orders_t': orders[:],
            'computed_cuts': ['zero'], 'computed_smith_orders': q * q,
            'infinity_profile_source': 'orthogonal equal-profile theorem; not separately computed',
            'seconds': time.monotonic() - started}


def packed_profile(profile):
    return {'q': profile['q'], 'zero_counts': sorted(Counter(profile['zero_orders_t']).items()),
            'infinity_counts': sorted(Counter(profile['infinity_orders_t']).items()),
            'computed_cuts': profile.get('computed_cuts', ['zero', 'infinity']),
            'computed_smith_orders': profile.get('computed_smith_orders', 2 * profile['q'] ** 2),
            'infinity_profile_source': profile.get('infinity_profile_source', 'directly computed'),
            'seconds': profile['seconds']}


def map_hash(vectors):
    return sha256(json.dumps([sorted(v.items()) for v in vectors], separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-classes', type=int, default=851)
    parser.add_argument('--seconds', type=float, default=10800)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--output', type=Path, default=Path('evidence/tensor-f16-coupled-quadratic-profiles.json'))
    args = parser.parse_args()
    all_classes = {canonical_frobenius(c) for c in product(range(1, 16), repeat=3)}
    assert len(all_classes) == 851
    priority = {canonical_frobenius(c) for c in product((1, 2, 4, 8), repeat=3)}
    queue = sorted(priority) + sorted(all_classes - priority)
    names = [Path(__file__).name, 'tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.cpp',
             'tensor_orthogonal_scaling.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'search_coherent_hot_retractions.py']
    hashes = {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names}
    if args.resume:
        report = json.loads(args.output.read_text())
        assert report['source_sha256'] == hashes
        report['in_progress'] = None
    else:
        report = {'scope': __doc__, 'field': 'F4[z]/(z^2+z+omega)', 'modes': 6, 'height_cap': 36,
                  'all_nonzero_triples': 3375, 'frobenius_profile_classes': 851,
                  'initial_basis_classes': len(priority), 'source_sha256': hashes,
                  'records': [], 'completed_classes': [], 'best_record_index': None,
                  'completed': False, 'in_progress': None}
    completed = {tuple(c) for c in report['completed_classes']}
    reference_path = Path('evidence/tensor-coupled-quadratic-profiles.json')
    known = json.loads(reference_path.read_text())
    reference = {(tuple(r['coefficients']), r['retraction_kind']): r for r in known['records']}
    report['f4_reference_sha256'] = sha256(reference_path.read_bytes()).hexdigest()
    previous = report.get('elapsed_seconds', 0)
    started = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        report['elapsed_seconds'] = previous + time.monotonic() - started
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':')) + '\n')
        temporary.replace(args.output)

    for coefficients in queue:
        if coefficients in completed:
            continue
        if len(completed) >= args.max_classes or time.monotonic() - started >= args.seconds:
            break
        a, b, c = coefficients
        code = field.differential_kernel([{9: 1, 6: a}, {33: 1, 18: b, 12: c}], 6, 36)
        assert code['q'] == 28
        if 'physical_heights' not in report:
            report['physical_heights'], report['physical_charges'] = code['heights'], code['charges']
        assert code['heights'] == report['physical_heights'] and code['charges'] == report['physical_charges']
        orthogonal = field.orthogonal_retraction(code)
        for kind, left in [('coordinate', code['retraction'])] + ([('orthogonal', orthogonal)] if orthogonal else []):
            if any(tuple(r['coefficients']) == coefficients and r['retraction_kind'] == kind for r in report['records']):
                continue
            report['in_progress'] = {'coefficients': coefficients, 'retraction_kind': kind}
            save()
            profile = (field.profiles(code['embedding'], left, 6) if kind == 'coordinate'
                       else orthogonal_profile(code, left))
            if (coefficients, kind) in reference:
                old = reference[coefficients, kind]['profile']
                assert profile['zero_orders_t'] == old['zero_orders_t']
                assert profile['infinity_orders_t'] == old['infinity_orders_t']
                report['embedded_f4_complete_profile_comparisons'] = report.get('embedded_f4_complete_profile_comparisons', 0) + 2
            fit = pressure(profile)
            report['records'].append({'coefficients': coefficients, 'retraction_kind': kind,
                                      'embedding_sha256': map_hash(code['embedding']), 'retraction_sha256': map_hash(left),
                                      'profile': packed_profile(profile), 'pressure': fit})
            best = report['best_record_index']
            if best is None or fit['rho'] < report['records'][best]['pressure']['rho']:
                report['best_record_index'] = len(report['records']) - 1
                report['best_maps'] = {'coefficients': coefficients, 'retraction_kind': kind, 'code': code, 'retraction': left}
                print(json.dumps({'new_best': coefficients, 'kind': kind, 'rho': fit['rho'],
                                  'seconds': profile['seconds'], 'classes_done': len(completed)}), flush=True)
            report['in_progress'] = None
            save()
        completed.add(coefficients)
        report['completed_classes'].append(coefficients)
        if len(completed) % 25 == 0:
            print(json.dumps({'classes_done': len(completed), 'records': len(report['records']),
                              'best_rho': report['records'][report['best_record_index']]['pressure']['rho'],
                              'elapsed_seconds': previous + time.monotonic() - started}), flush=True)
        save()
    report['completed'] = len(completed) == 851
    report['stopped_reason'] = 'finite catalogue exhausted' if report['completed'] else 'class or time budget reached'
    save()
    print(json.dumps({'completed': report['completed'], 'classes_done': len(completed),
                      'records': len(report['records']), 'output': str(args.output)}), flush=True)


if __name__ == '__main__':
    main()
