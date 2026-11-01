#!/usr/bin/env python3
"""All-scale product-pressure domination of six height-20 orthogonal maps.

The exact certificate compares complete saved profiles. --recompute-control
reconstructs both coordinate-control profiles; --recompute-all additionally
reconstructs all six orthogonal profiles. No floating pressure is used.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from coherent_differential_kernels_f4 import differential_kernel
from coherent_hot_compounds_f4 import profiles
from orthogonal_differential_pairs_f4 import orthogonal_retraction
from tensor_orthogonal_pair_profiles import normalized_scalings


def product_polynomial(profile):
    q = profile['q']
    assert len(profile['zero_orders_t']) == len(profile['infinity_orders_t']) == q * q
    result = Counter()
    for a, x in Counter(profile['zero_orders_t']).items():
        for b, y in Counter(profile['infinity_orders_t']).items():
            result[a + b] += x * y
    return result


def factor_difference(P, Q):
    difference = [P[i] - Q[i] for i in range(1 + max(max(P), max(Q)))]
    shift = next(i for i, x in enumerate(difference) if x)
    coefficients = difference[shift:]
    one_minus_power = 0
    while sum(coefficients) == 0:
        quotient, partial = [], 0
        for x in coefficients[:-1]:
            partial += x
            quotient.append(partial)
        assert partial == -coefficients[-1]
        coefficients = quotient
        one_minus_power += 1
    assert shift == 44 and one_minus_power == 1 and len(coefficients) == 163
    assert all(x > 0 for x in coefficients)
    # Reconstruct the original integer difference coefficient by coefficient.
    reconstructed = coefficients[:]
    for _ in range(one_minus_power):
        next_poly = [0] * (len(reconstructed) + 1)
        for i, x in enumerate(reconstructed):
            next_poly[i] += x
            next_poly[i + 1] -= x
        reconstructed = next_poly
    reconstructed = [0] * shift + reconstructed
    assert reconstructed == difference
    return {'monomial_power': shift, 'one_minus_power': one_minus_power,
            'quotient_degree': len(coefficients) - 1, 'positive_coefficients': coefficients}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profiles', type=Path, default=Path('evidence/tensor-orthogonal-pair-profiles-small.json'))
    ap.add_argument('--output', type=Path, default=Path('evidence/tensor-orthogonal-pressure-dominance.json'))
    ap.add_argument('--recompute-control', action='store_true')
    ap.add_argument('--recompute-all', action='store_true')
    args = ap.parse_args()
    saved = json.loads(args.profiles.read_text())
    assert saved['completed'] and saved['in_progress'] is None
    trials = [record for record in saved['records'] if record['height_cap'] == 20]
    expected = list(normalized_scalings())
    assert len(trials) == 6 and {tuple(t['scaling']) for t in trials} == {s for s, _, _ in expected}
    pairs = [(0, 7), (1, 6), (2, 5), (3, 4)]
    C = [[1, 1, 1, 1], [0, 1, 2, 3]]
    operators = [{(1 << a) | (1 << b): row[i] for i, (a, b) in enumerate(pairs) if row[i]}
                 for row in C]
    code = differential_kernel(operators, 12, 20)
    if args.recompute_control or args.recompute_all:
        print('Recomputing both coordinate-control profiles.', flush=True)
        control = profiles(code['embedding'], code['retraction'], 12)
    else:
        if not args.output.exists():
            raise FileNotFoundError('Use --recompute-control to build the initial certificate.')
        control = json.loads(args.output.read_text())['coordinate_control_profile']
    assert code['q'] == control['q'] == 48
    assert Counter(control['heights']) == Counter(code['heights'])
    Q = product_polynomial(control)
    certificates = []
    recomputed_orders = 2 * 48 * 48 if args.recompute_control or args.recompute_all else 0
    for trial in trials:
        profile = trial['profile']
        if args.recompute_all:
            stored = trial['code']
            code2 = differential_kernel([{int(k): v for k, v in op.items()} for op in stored['operators']], 12, 20)
            left = orthogonal_retraction(code2)
            fresh = profiles(code2['embedding'], left, 12)
            for cut in ('zero_orders_t', 'infinity_orders_t'):
                assert fresh[cut] == profile[cut]
            recomputed_orders += 2 * 48 * 48
        assert profile['zero_orders_t'] == profile['infinity_orders_t']
        assert Counter(profile['heights']) == Counter(control['heights'])
        certificate = factor_difference(product_polynomial(profile), Q)
        certificates.append({'scaling': trial['scaling'], **certificate})
    report = {'scope': __doc__, 'passed': True,
              'coordinate_control_coefficients': C, 'coordinate_control_profile': control,
              'orthogonal_maps': len(trials), 'physical_letters': 48,
              'exact_differences': certificates,
              'input_profile_artifact': args.profiles.as_posix(),
              'input_profile_sha256': sha256(args.profiles.read_bytes()).hexdigest(),
              'recomputed_smith_orders_this_run': recomputed_orders,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest()
                                for n in (Path(__file__).name, 'tensor_orthogonal_pair_profiles.py',
                                          'coherent_differential_kernels_f4.py', 'coherent_hot_compounds_f4.py',
                                          'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
                                          'native_rectangular_smith_f4.cpp')}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'orthogonal_maps': 6, 'quotient_degree': 162,
                      'all_coefficients_strictly_positive': True,
                      'recomputed_smith_orders_this_run': recomputed_orders}))


if __name__ == '__main__':
    main()
