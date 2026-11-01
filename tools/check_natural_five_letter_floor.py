#!/usr/bin/env python3
"""Exact global cell-cap floor for paired energies proportional to 0,2,6,8,10.

Uses existing cell/shared-minor inequalities. This is a balanced criterion
obstruction, not a stronger factoring algorithm or an asymmetric exclusion.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial
import json
import lzma
from pathlib import Path
import time

from charged_gaussian_floor import Dyadic
from check_cell_minor_caps import log_at_least_one, partition, point, relaxed_orders
from check_quadratic_crossings import bernstein_certificate

ENERGIES = tuple(map(F, (0, 2, 6, 8, 10)))
R = F(10, 7)  # Twice the claimed balanced floor rho=5/7.
SMALL = F(1, 10000)
ROOT = (F(0), F(5), SMALL, F(7))


def corner(arithmetic, box):
    ul, uh, vl, vh = box
    orders, _, _ = relaxed_orders(uh, [vh*e for e in ENERGIES], [vh*e for e in ENERGIES])
    raw = arithmetic.add(point(arithmetic, ul), arithmetic.multiply(
        point(arithmetic, 2), log_at_least_one(arithmetic, partition(arithmetic, orders))))
    physical = arithmetic.multiply(point(arithmetic, 2*R), log_at_least_one(
        arithmetic, partition(arithmetic, [vl*e/R for e in ENERGIES])))
    return arithmetic.subtract(raw, physical)


def tails(arithmetic):
    # Existing all-x>=1 shared-minor certificate U5(x)>=1029/100.
    polynomial = [0]*10
    polynomial[9], polynomial[5], polynomial[4], polynomial[0] = 100, -1029, 400, 2000
    transformed = [sum(polynomial[i]*comb(i, j)*2**j for i in range(j, 10)) for j in range(10)]
    bernstein = bernstein_certificate(transformed)
    assert 100*3**4 > 1029
    log_t = log_at_least_one(arithmetic, point(arithmetic, F(1029, 100)))
    log_five = log_at_least_one(arithmetic, point(arithmetic, 5))
    small_margin = arithmetic.subtract(arithmetic.subtract(
        arithmetic.multiply(point(arithmetic, 2), log_t),
        arithmetic.multiply(point(arithmetic, 2*R), log_five)), point(arithmetic, 40*SMALL))
    large_u_margin = arithmetic.subtract(point(arithmetic, 5),
        arithmetic.multiply(point(arithmetic, 2*R), log_five))
    large_v_comparison = sum(F(14, 5)**k/factorial(k) for k in range(5))-8*R
    assert small_margin[0] > 0 and large_u_margin[0] > 0 and large_v_comparison > 0
    assert F(5, 7)**2 > F(1, 2)
    return {'shared_min_orbit_polynomial': polynomial, 'shared_min_bernstein': bernstein,
            'small_v_threshold': str(SMALL), 'small_v_lower_margin_dyadic': small_margin,
            'large_u_threshold': '5', 'large_u_lower_margin_dyadic': large_u_margin,
            'large_v_threshold': '7', 'large_v_taylor_comparison': str(large_v_comparison)}


def source_and_scope_controls(arithmetic):
    from check_joint_energy_conjugations import cut_matrix
    from joint_energy_conjugations import model
    from search_f4_recentered_crossings import five_letter_source
    from weighted_smith_oracle_f4 import WeightedSmithOracle

    source, _, ex, ey = five_letter_source()
    assert sorted(ex) == sorted(ey) == list(ENERGIES)
    weights = [ex[i//5]//2 for i in range(25)]
    infinity, _ = WeightedSmithOracle(cut_matrix(source, True), infinity=True).profile(weights, weights)
    zero, _ = WeightedSmithOracle(cut_matrix(source, False)).profile(model()['energies'], [0]*25)
    assert zero[0]+infinity[0] == 0
    prefix = [sum(zero[:k])+sum(infinity[:k]) for k in range(1, 6)]
    assert prefix == [sum(ENERGIES[:k]) for k in range(1, 6)]
    average = [(a-zero[0]+b-infinity[0])/F(2) for a, b in zip(zero, infinity)]
    ceiling, _, _ = relaxed_orders(F(0), ENERGIES, ENERGIES)
    assert all(sum(average[:k]) <= sum(ceiling[:k]) for k in range(1, 26))
    # The known negative five-letter target has a different paired energy
    # shape. It must remain outside the theorem's scope.
    other = tuple(map(F, (0, 0, 2, F(8, 3), 3)))
    orders, _, _ = relaxed_orders(F(1), other, other)
    negative = arithmetic.subtract(arithmetic.add(point(arithmetic, 1), arithmetic.multiply(
        point(arithmetic, 2), log_at_least_one(arithmetic, partition(arithmetic, orders)))),
        arithmetic.multiply(point(arithmetic, 2*R), log_at_least_one(
            arithmetic, partition(arithmetic, [v/R for v in other]))))
    assert negative[1] < 0
    return {'independent_full_fraction_free_orders': 50,
            'source_complete_zero_orders': zero, 'source_complete_infinity_orders': infinity,
            'first_five_combined_prefixes': prefix,
            'reference_average_normalized_orders': list(map(str, average)),
            'd0_cell_cap_ceiling': list(map(str, ceiling)),
            'different_energy_shape_negative_control': list(map(str, other)),
            'different_energy_shape_pressure_interval': negative}


def split(box, axis):
    ul, uh, vl, vh = box
    if axis == 'u':
        middle = (ul+uh)/2
        return (ul, middle, vl, vh), (middle, uh, vl, vh)
    middle = (vl+vh)/2
    return (ul, uh, vl, middle), (ul, uh, middle, vh)


def verify_cover(arithmetic, leaves, recheck=True):
    directions, endpoints, kraft = {}, set(), F(0)
    minimum = None
    for leaf in leaves:
        box, bits = ROOT, ''
        for code in leaf['path']:
            assert code in '0123'
            axis = 'u' if code in '01' else 'v'
            assert bits not in endpoints
            if bits in directions:
                assert directions[bits] == axis
            else:
                directions[bits] = axis
            side = int(code) % 2
            box = split(box, axis)[side]
            bits += str(side)
        assert bits not in endpoints and bits not in directions
        endpoints.add(bits)
        kraft += F(1, 1 << len(bits))
        if recheck:
            bound = corner(arithmetic, box)
            assert bound[0] > 0 and bound[0] == leaf['lower_dyadic']
        assert leaf['lower_dyadic'] > 0
        minimum = leaf['lower_dyadic'] if minimum is None else min(minimum, leaf['lower_dyadic'])
    assert kraft == 1
    return minimum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--max-leaves', type=int, default=250000)
    args = parser.parse_args()
    arithmetic, start = Dyadic(80), time.monotonic()
    tail_controls = tails(arithmetic)
    source_controls = source_and_scope_controls(arithmetic)
    if args.replay:
        data = args.replay.read_bytes()
        report = json.loads(lzma.decompress(data) if args.replay.suffix == '.xz' else data)
        assert report['status'] == 'pass' and report['precision'] == arithmetic.precision
        assert report['tails'] == json.loads(json.dumps(tail_controls))
        assert report['source_and_scope_controls'] == json.loads(json.dumps(source_controls))
        verify_cover(arithmetic, report['leaves'])
        print(json.dumps({'status': 'pass', 'replayed_compact_leaves': len(report['leaves'])}))
        return
    queue, leaves, visits, max_depth = [(ROOT, '')], [], 0, 0
    while queue:
        box, path = queue.pop()
        visits += 1
        bound = corner(arithmetic, box)
        if bound[0] > 0:
            leaves.append({'path': path, 'lower_dyadic': bound[0]})
            max_depth = max(max_depth, len(path))
        else:
            assert len(path) < 60 and len(leaves)+len(queue) < args.max_leaves, 'INCOMPLETE box budget'
            ul, uh, vl, vh = box
            axis = 'u' if uh-ul > 4*(vh-vl) else 'v'
            left, right = split(box, axis)
            base = 0 if axis == 'u' else 2
            queue.append((right, path+str(base+1)))
            queue.append((left, path+str(base)))
        if visits % 10000 == 0:
            print(json.dumps({'visited_boxes': visits, 'certified_leaves': len(leaves),
                              'pending_boxes': len(queue), 'seconds': time.monotonic()-start}), flush=True)
    minimum = verify_cover(arithmetic, leaves, recheck=False)
    root = Path(__file__).resolve().parents[1]
    names = ['tools/check_natural_five_letter_floor.py', 'tools/check_cell_minor_caps.py',
             'tools/charged_gaussian_floor.py', 'tools/check_quadratic_crossings.py']
    names += ['tools/check_joint_energy_conjugations.py', 'tools/joint_energy_conjugations.py',
              'tools/search_f4_recentered_crossings.py', 'tools/weighted_smith_oracle_f4.py']
    report = {'status': 'pass', 'balanced_rho_floor': '5/7', 'paired_energy_shape': list(map(str, ENERGIES)),
              'precision': arithmetic.precision, 'compact_root': list(map(str, ROOT)),
              'leaves': leaves, 'compact_leaves': len(leaves), 'maximum_depth': max_depth,
              'minimum_compact_lower_dyadic': minimum, 'tails': tail_controls,
              'source_and_scope_controls': source_controls,
              'source_sha256': {n: sha256((root/n).read_bytes()).hexdigest() for n in names},
              'seconds': time.monotonic()-start,
              'scope': 'Every scalar Laurent-polynomial strict crossing over any field with both '
              'paired physical energies proportional by one positive common scale to (0,2,6,8,10), '
              'arbitrary independent height splitting, degree and regular local frames. Balanced '
              'two-place Smith criterion only; other energy shapes, asymmetric criteria and auxiliary '
              'modules remain open. No stronger factoring bound.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({key: report[key] for key in
                      ('status', 'balanced_rho_floor', 'compact_leaves', 'maximum_depth', 'seconds')}))


if __name__ == '__main__':
    main()
