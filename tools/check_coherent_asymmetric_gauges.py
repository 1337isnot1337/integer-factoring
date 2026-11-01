#!/usr/bin/env python3
"""Literal profiles, shared limits, and derivatives for independent cut gauges."""
import argparse
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import numpy as np
from check_coherent_shared_gauges import literal_profiles
from coherent_asymmetric_gauges import AsymmetricGaugeOracle
from coherent_graded_quotients import parameter_images, quotient_maps
from coherent_hot_compounds import support
from coherent_inverse_systems import inverse_system
from coherent_weighted_gauges import SharedGaugeOracle
from search_coherent_core_products import core_maps


def controls():
    inverse = inverse_system([(1 << 6) ^ (1 << 9)], 4)
    quotient = quotient_maps(parameter_images(3, 5, 1), 3)
    cases = [('inverse_q6', inverse['embedding'], inverse['retraction'], 4),
             ('coherent_q8', *core_maps(8, True), 4),
             ('nonlinear_q8', quotient['embedding'], quotient['retraction'], 5)]
    rng = np.random.default_rng(2026110131)
    counts, records = Counter(), []
    for name, embedding, retraction, modes in cases:
        oracle = AsymmetricGaugeOracle(embedding, retraction, modes)
        shared = SharedGaugeOracle(embedding, retraction, modes)
        first = []
        for p, l, h in zip(embedding, retraction, oracle.heights.astype(int)):
            shifts = {sum(1 if a % 2 == 0 else -1 for a in range(modes)
                          if word >> a & 1) for word in support(p | l)}
            assert len(shifts) == 1
            first.append(int(h + 8 * shifts.pop()))
        second = [x + int(h) for x, h in zip(first, oracle.heights)]
        point = np.array([.5] + [first[j] / 32 for j in oracle.free]
                         + [second[j] / 32 for j in oracle.free])
        actual = oracle.evaluate(point, details=True)
        assert min(actual['physical_energies'][j] for j in oracle.free) > 0
        profiles = []
        for place, nums, key in ((0, first, 'zero_orders'), (1, second, 'infinity_orders')):
            exact, _ = literal_profiles(embedding, retraction, modes, nums)
            scaled = [32 * value for value in actual[key]]
            assert all(value == round(value) for value in scaled)
            assert exact[place] == sorted(map(round, scaled))
            profiles.append(exact[place])
            counts['complete_literal_cut_profiles'] += 1
            counts['complete_smith_orders'] += len(exact[place])
        tied = np.array([.5] + [first[j] / 32 for j in oracle.free] * 2)
        joined, common = oracle.evaluate(tied, details=True), shared.evaluate(tied[:oracle.q], details=True)
        assert abs(joined['rho'] - common['rho']) < 1e-12
        assert np.max(np.abs(np.array(joined['physical_type']) - np.array(common['physical_type']))) < 1e-12
        assert abs(joined['gradient'][0] - common['gradient'][0]) < 1e-12
        assert np.max(np.abs(joined['gradient'][1:oracle.q] + joined['gradient'][oracle.q:]
                             - common['gradient'][1:])) < 1e-11
        counts['shared_gauge_limit_controls'] += 1
        interior = point.copy()
        interior[1:] += rng.uniform(-.0001, .0001, len(point) - 1)
        at = oracle.evaluate(interior)
        assert at['minimum_cone_slack'] > 0
        errors = []
        for _ in range(6):
            direction = rng.normal(size=len(point)); direction /= np.linalg.norm(direction)
            step = 1e-6
            numerical = (oracle.evaluate(interior + step * direction)['rho']
                         - oracle.evaluate(interior - step * direction)['rho']) / (2 * step)
            error = abs(numerical - float(at['gradient'] @ direction))
            assert error < 1e-6, (name, error)
            errors.append(error); counts['directional_gradient_controls'] += 1
        records.append({'name': name, 'modes': modes, 'zero_numerators_over_32': first,
                        'infinity_numerators_over_32': second, 'scale': .5,
                        'exact_ramified_profiles': profiles,
                        'minimum_cone_slack': actual['minimum_cone_slack'],
                        'numerical_rho': actual['rho'], 'gradient_absolute_errors': errors})
    return dict(counts), records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    counts, records = controls()
    names = (Path(__file__).name, 'coherent_asymmetric_gauges.py',
             'coherent_weighted_gauges.py', 'check_coherent_shared_gauges.py')
    report = {'status': 'pass', 'counts': counts, 'records': records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'rational_profiles_exact': True,
              'gradient_checks_numerical': True, 'factoring_exponent_improved': False}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'pass', 'counts': counts}, indent=2))


if __name__ == '__main__':
    main()
