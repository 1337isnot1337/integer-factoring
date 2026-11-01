#!/usr/bin/env python3
"""Replay complete unequal coherent tensors and exact supporting-line tests."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import time
import numpy as np

from certify_hot_charge_pressure import homogeneous_sum, integer_digest
from check_coherent_hot_boundary import full_maps, logical
from hot_code_charge_profiles import hot_energy
from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import retract
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from screen_unequal_coherent_hot_codes import natural_parameters, support


def exact_comparison(row):
    record = row['exact_point']
    p, q, h = record['weights']
    n, d = record['base']
    assert p*p + q*q == h*h and 0 < n < d
    infinity, zero = row['regular_profiles_infinity_zero']
    pi, di = homogeneous_sum([q*h*x for x in infinity], n, d)
    pz, dz = homogeneous_sum([p*h*x for x in zero], n, d)
    qx, dx = homogeneous_sum([p*q*x for x in row['energies_x']], n, d)
    qy, dy = homogeneous_sum([p*q*x for x in row['energies_y']], n, d)
    left, right = pi**p * pz**q, (qx*qy)**h
    difference = p*di + q*dz - h*(dx+dy)
    if difference >= 0:
        right *= d**difference
    else:
        left *= d**(-difference)
    assert record['strict_support_below_one'] == (left < right)
    assert record['integer_bits'] == [left.bit_length(), right.bit_length()]
    assert record['left_sha256'] == integer_digest(left)
    assert record['right_sha256'] == integer_digest(right)
    assert record['absolute_difference_sha256'] == integer_digest(abs(left-right))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('evidence/unequal-coherent-hot-screen.json'))
    parser.add_argument('--output', type=Path, default=Path('evidence/unequal-coherent-hot-screen-replay.json'))
    args = parser.parse_args()
    began = time.monotonic()
    report = json.loads(args.input.read_text())
    assert report['completed'] and len(report['records']) == 18
    for name, expected in report['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected
    seen, endpoint_count, channel_count, cost_count = set(), 0, 0, 0
    for row in report['records']:
        key = row['source_size'], row['x_kind'], row['y_kind']
        assert key not in seen
        seen.add(key)
        lower, boundary = row['lower_letters'], row['boundary_letters']
        raw = restricted_tensor(16, list(map(logical, lower+boundary)))
        lx, rx = full_maps(len(lower), len(boundary), row['x_left'], row['x_right'])
        ly, ry = full_maps(len(lower), len(boundary), row['y_left'], row['y_right'])
        seed = retract(raw, lx, rx, ly, ry)
        content = min((part & -part).bit_length()-1 for value in seed['values'] for part in value if part)
        assert content == row['scalar_content']
        seed['values'] = [(a>>content, b>>content) for a,b in seed['values']]
        seed['degree'] = max(max(a.bit_length(), b.bit_length())-1 for a,b in seed['values'])
        assert digest(seed['values']) == row['tensor_sha256'] and seed['degree'] == row['degree']
        ex, ey = row['energies_x'], row['energies_y']
        assert ex == list(map(hot_energy,lower))+[hot_energy(boundary[0])]*len(row['x_left'])
        assert ey == list(map(hot_energy,lower))+[hot_energy(boundary[0])]*len(row['y_left'])
        qx, qy = seed['q_x'], seed['q_u']
        assert (qx,qy) == (row['q_x'],row['q_y']) and qx != qy
        values = RectangularOracle(seed['values'],qx,qy).prefixes(natural_parameters(ex,ey))['orders']
        assert len(values) == 2*qx*qy
        sides = [sorted(values[:qx*qy]),sorted(values[qx*qy:])]
        assert [min(side) for side in sides] == row['opposite_shifts']
        assert sum(row['opposite_shifts']) == 0
        assert [[v-min(side) for v in side] for side in sides] == row['regular_profiles_infinity_zero']
        assert len(row['branch_trials']) == 6 and all(t['completed'] for t in row['branch_trials'])
        for trial in row['branch_trials']:
            value = support(row['regular_profiles_infinity_zero'],ex,ey,trial['scale'],*trial['weights'])
            assert abs(value-trial['support']) < 1e-10
            cost_count += 1
        assert row['best_numerical_support'] == min(row['branch_trials'],key=lambda t:t['support'])
        exact_comparison(row)
        endpoint_count += len(seed['values'])
        channel_count += len(values)
    result = {'passed':True,'scope':'Complete exact tensors/profiles and retained rational tests; numerical scale optima are not exclusions.',
              'input_sha256':sha256(args.input.read_bytes()).hexdigest(),'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
              'tensors':len(seen),'endpoints':endpoint_count,'integer_profile_channels':channel_count,
              'numerical_costs':cost_count,'exact_rational_comparisons':len(seen),
              'strict_tests_passed':sum(r['exact_point']['strict_support_below_one'] for r in report['records']),
              'seconds':time.monotonic()-began}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
