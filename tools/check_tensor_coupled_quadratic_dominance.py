#!/usr/bin/env python3
"""Exact all-scale comparison of the completed coupled-quadratic F4 screen.

This replays integer polynomial identities from saved complete profiles;
it does not recompute their Smith orders.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path


def profile_product(profile):
    P = Counter()
    for a, x in Counter(profile['zero_orders_t']).items():
        for b, y in Counter(profile['infinity_orders_t']).items():
            P[a + b] += x * y
    return P


def positive_difference(P, Q):
    original = [P[i] - Q[i] for i in range(1 + max(max(P), max(Q)))]
    while original and not original[-1]:
        original.pop()
    if not original:
        return {'comparison': 'equal'}
    shift = next(i for i, coefficient in enumerate(original) if coefficient)
    coefficients = original[shift:]
    one_minus_power = 0
    while sum(coefficients) == 0:
        quotient, partial = [], 0
        for coefficient in coefficients[:-1]:
            partial += coefficient
            quotient.append(partial)
        assert partial == -coefficients[-1]
        coefficients = quotient
        one_minus_power += 1
    degree = len(coefficients) - 1
    if min(coefficients) >= 0:
        basis, nonnegative = 'power', coefficients[:]
    else:
        basis = 'bernstein'
        nonnegative = [sum(coefficients[i] * comb(degree - i, j - i) for i in range(j + 1))
                       for j in range(degree + 1)]
        rebuilt = [0] * (degree + 1)
        for j, value in enumerate(nonnegative):
            for k in range(degree - j + 1):
                rebuilt[j + k] += value * (-1 if k & 1 else 1) * comb(degree - j, k)
        assert rebuilt == coefficients
    assert min(nonnegative) >= 0 and max(nonnegative) > 0
    rebuilt = coefficients[:]
    for _ in range(one_minus_power):
        new = [0] * (len(rebuilt) + 1)
        for i, value in enumerate(rebuilt):
            new[i] += value
            new[i + 1] -= value
        rebuilt = new
    assert [0] * shift + rebuilt == original
    return {'comparison': 'strictly_larger', 'monomial_power': shift,
            'one_minus_power': one_minus_power, 'quotient_degree': degree,
            'power_coefficients': coefficients, 'nonnegative_basis': basis,
            'nonnegative_coefficients': nonnegative}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profiles', type=Path, default=Path('evidence/tensor-coupled-quadratic-profiles.json'))
    ap.add_argument('--output', type=Path, default=Path('evidence/tensor-coupled-quadratic-dominance.json'))
    args = ap.parse_args()
    data = json.loads(args.profiles.read_text())
    assert data['completed'] and data['in_progress'] is None
    assert data['modes'] == 6 and data['height_cap'] == 36
    records = data['records']
    coordinate = [r for r in records if r['retraction_kind'] == 'coordinate']
    orthogonal = [r for r in records if r['retraction_kind'] == 'orthogonal']
    assert {tuple(r['coefficients']) for r in coordinate} == set(product((1, 2, 3), repeat=3))
    assert len(coordinate) == 27 and len(orthogonal) == 4
    assert {tuple(r['coefficients']) for r in orthogonal} == {(2, 1, 2), (2, 1, 3), (3, 1, 2), (3, 1, 3)}
    control = next(r for r in coordinate if r['coefficients'] == [1, 1, 1])
    physical = Counter(control['profile']['heights'])
    reference = profile_product(control['profile'])
    certificates, methods = [], Counter()
    computed_orders = 0
    for record in records:
        p = record['profile']
        assert p['q'] == 28 and len(p['zero_orders_t']) == len(p['infinity_orders_t']) == 784
        assert Counter(p['heights']) == physical
        computed_orders += p.get('computed_smith_orders', 2 * 784)
        certificate = positive_difference(profile_product(p), reference)
        methods[certificate.get('nonnegative_basis', 'equal')] += 1
        certificates.append({'coefficients': record['coefficients'],
                             'retraction_kind': record['retraction_kind'], **certificate})
    assert methods == {'equal': 1, 'power': 28, 'bernstein': 2}
    assert computed_orders == 45472
    result = {'scope': __doc__, 'passed': True, 'profile_records': len(records),
              'comparison_counts': dict(methods), 'native_smith_orders_in_source_run': computed_orders,
              'derived_infinity_orders_in_source_run': 4 * 784,
              'native_smith_orders_recomputed_by_this_checker': 0,
              'input_profile_path': args.profiles.as_posix(),
              'input_profile_sha256': sha256(args.profiles.read_bytes()).hexdigest(),
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(), 'certificates': certificates}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'certificates'}))


if __name__ == '__main__':
    main()
