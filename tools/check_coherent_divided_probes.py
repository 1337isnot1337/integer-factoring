#!/usr/bin/env python3
"""Replay retained divided-power maps, fixed fits and 48 exact rational points."""
import argparse
from hashlib import sha256
import json
from math import cos, sin
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from certify_hot_charge_pressure import homogeneous_sum, integer_digest
from coherent_divided_powers import paired_top
from coherent_inverse_systems import inverse_system
from screen_unequal_coherent_hot_codes import support

ROOT = Path(__file__).resolve().parents[1]


def exact_replay(report, point):
    p, q, h = point['weights']
    n, d = point['base']
    assert p*p+q*q == h*h and 0 < n < d
    profile = report['profile']
    ex, ey = ([2*x for x in profile[key]] for key in ('heights_x', 'heights_y'))
    pi, di = homogeneous_sum([q*h*x for x in profile['infinity_orders_t']], n, d)
    pz, dz = homogeneous_sum([p*h*x for x in profile['zero_orders_t']], n, d)
    qx, dx = homogeneous_sum([p*q*x for x in ex], n, d)
    qy, dy = homogeneous_sum([p*q*x for x in ey], n, d)
    left, right = pi**p*pz**q, (qx*qy)**h
    difference = p*di+q*dz-h*(dx+dy)
    if difference >= 0:
        right *= d**difference
    else:
        left *= d**(-difference)
    assert point['strict_support_below_one'] == (left < right) == False
    assert point['integer_bits'] == [left.bit_length(), right.bit_length()]
    assert point['left_sha256'] == integer_digest(left)
    assert point['right_sha256'] == integer_digest(right)
    assert point['absolute_difference_sha256'] == integer_digest(abs(left-right))


def certificate():
    path = ROOT/'evidence/coherent-divided-power-probes.json'
    data = json.loads(path.read_text())
    records, points, order_count = [data['equal']]+data['unequal'], 0, 0
    for report in records:
        for name, digest in report['source_sha256'].items():
            assert sha256((ROOT/'tools'/name).read_bytes()).hexdigest() == digest, name
        powers = report.get('powers', [report.get('power')])
        maps = [inverse_system([sum(1 << word for word in paired_top(report['modes'], r))], report['modes'])
                for r in powers]
        assert maps == report.get('codes', [report.get('code')])
        profile = report['profile']
        size = profile['q']**2 if 'q' in profile else profile['q_x']*profile['q_y']
        for name in ('zero_orders_t', 'infinity_orders_t'):
            assert len(profile[name]) == size and profile[name] == sorted(profile[name])
        for cut in profile['cuts']:
            assert sorted(value for block in cut['blocks'] for value in block['orders_t']) == cut['orders_t']
            for block in cut['blocks']:
                assert block['precision_attempts'][-1]['status'] == 0
                assert block['endpoint_entries'] == block['size']**2
        order_count += 2*size
        if 'q' in profile:
            rho, t = report['pressure']['rho'], report['pressure']['scale_t']
            residual = sum(float(logsumexp(-t*np.asarray(profile[key])))
                           for key in ('zero_orders_t', 'infinity_orders_t'))
            residual -= 4*rho*float(logsumexp(-t*np.asarray(profile['heights'])/rho))
            assert abs(residual) < 1e-9
        else:
            best = report['best_numerical']
            value = support([profile['infinity_orders_t'], profile['zero_orders_t']],
                            [2*x for x in profile['heights_x']], [2*x for x in profile['heights_y']],
                            best['scale'], cos(best['angle']), sin(best['angle']), 1)
            assert abs(value-best['support']) < 1e-9
            for point in report['exact_rational_tests']:
                exact_replay(report, point)
                points += 1
    names = (Path(__file__).name, 'coherent_divided_powers.py', 'coherent_inverse_systems.py',
             'certify_hot_charge_pressure.py', 'screen_unequal_coherent_hot_codes.py')
    return dict(retained_pairs=len(records), exact_rational_points_replayed=points,
                recorded_complete_orders=order_count, large_native_profiles_recomputed=False,
                input_sha256=sha256(path.read_bytes()).hexdigest(),
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names},
                scope='Maps, complete-profile receipts, fixed numerical points and all retained exact integer comparisons replayed. Large native evaluations and global optimization are not repeated.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/coherent-divided-probes-replay.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({key: value for key, value in result.items() if key != 'source_sha256'}))


if __name__ == '__main__':
    main()
