#!/usr/bin/env python3
"""Independent finite-field and full-tensor controls for Hasse frame fibers."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import random
import time

import numpy as np

from affine_place_profile_oracle import AffinePlaceOracle
from check_bell_jet_crt_crossings import hasse
from check_dense_coupled_reflection_frames import coefficients, multiply
from entrywise_reflection_constraints import field_dot, hasse_row
from finite_jet_smith_fibers import FiniteJetFiber, insert
from finite_precision_smith_frames import build_lattice
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import digest
from search_finite_jet_smith_fibers import mixed_fixed_gauge_score
from search_fused_mixed_weighted import MixedPressure


def pack(values):
    return tuple(sum(((value >> plane) & 1) << index for index, value in enumerate(values))
                 for plane in range(2))


def independent_dot(row, values):
    result = 0
    for coefficient, value in zip(row, values):
        result ^= multiply(coefficient, value)
    return result


def cost_replay(seed, lattice, values):
    profiles = lattice.profiles(values)
    oracle = AffinePlaceOracle(seed)
    finite = oracle.log_profile_with_gradient(0, (2, 1), [0] * 16)[3]
    score = mixed_fixed_gauge_score(profiles[1], finite, [value + 24 for value in profiles[0]])
    t, alpha = score['outer_log_scale'], score['row_allocations'][0]
    probabilities = np.asarray(score['paired_type'])
    baseline = {'row_places': [0, 2], 'column_places': [None],
                'row_allocations': [alpha, 1 - alpha], 'column_allocations': [1.],
                'outer_log_scale': t}
    pressure = MixedPressure(oracle, baseline, 'full', False, True, time.monotonic() + 60)
    parameters = np.zeros(pressure.dimension)
    parameters[pressure.scale_index] = t
    h, zero = np.asarray([0, 6, 2, 8]), np.zeros(4)
    at_zero = np.concatenate((h, h, zero, zero)) * (t / alpha)
    at_infinity = np.concatenate((zero, zero, h, h)) * t
    parameters[pressure.height_slices[0][0]] = at_zero[pressure.active_heights]
    parameters[pressure.height_slices[1][0]] = at_infinity[pressure.active_heights]
    parameters[pressure.type_slice] = np.tile(np.log(probabilities[1:] / probabilities[0]), 2)
    parameters[pressure.allocation_slices[0]] = math.log(alpha / (1 - alpha))
    pressure.evaluate(parameters, .001)
    assert abs(score['rho'] - pressure.best['rho']) < 2e-9
    return {'fast_rho': score['rho'], 'complete_oracle_rho': pressure.best['rho'],
            'absolute_difference': abs(score['rho'] - pressure.best['rho'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(312529)
    systems, assignments = [], 0
    for trial in range(24):
        width = 2 + trial % 4
        equations = [([rng.randrange(4) for _ in range(width)], rng.randrange(4))
                     for _ in range(1 + trial % 7)]
        expected = {values for values in product(range(4), repeat=width)
                    if all(independent_dot(row, values) == rhs for row, rhs in equations)}
        assignments += 4 ** width
        pivots, feasible = {}, True
        for row, rhs in equations:
            if not insert(pivots, pack(row), rhs):
                feasible = False
                break
        assert feasible == bool(expected)
        if feasible:
            actual = {values for values in product(range(4), repeat=width)
                      if all(field_dot(row, pack(values)) == rhs for row, rhs in pivots.values())}
            assert actual == expected and len(actual) == 4 ** (width - len(pivots))
        systems.append({'variables': width, 'equations': len(equations), 'solutions': len(expected)})
    base = tensor(4)
    h, zero = [0, 6, 2, 8], [0] * 4
    lattice = build_lattice(base['values'], 4, 24, [zero, zero, h, h, h, h, zero, zero], 8)
    floors = tuple(max(0, value - 6) for value in lattice.orders)
    fiber = FiniteJetFiber(lattice, floors)
    assert fiber.variables == 394
    full_base_values = [(entry, 2, 0, field_dot(hasse_row(2, 0, 24), value))
                        for entry, value in enumerate(base['values'])]
    partial_zero = [(entry, 2, 0, 0) for entry, bits in enumerate(lattice.labels) if bits[0] != bits[2]]
    cases = [(lattice.orders, []), (lattice.orders, full_base_values),
             (floors, partial_zero),
             (floors, full_base_values[:64] + [(entry, 3, 1, field_dot(hasse_row(3, 1, 24), base['values'][entry]))
                                              for entry in range(32)])]
    records, finite_entries, endpoint_entries = [], 0, 0
    for targets, requests in cases:
        equations = [fiber.jet_equation(*request) for request in requests]
        pivots = fiber.solve(equations, targets)
        assert pivots is not None
        for sample in range(2):
            values = fiber.sample(pivots, None if sample == 0 else rng)
            profiles = lattice.check_witness(values, targets)
            for entry, point, order, target in requests:
                assert hasse(coefficients(values[entry]), point, order) == target
                finite_entries += 1
            for entry, value in enumerate(values):
                x, y, u, v = lattice.labels[entry]
                assert hasse(coefficients(value), 1, 0) == int(x == v and y == u)
                endpoint_entries += 1
            seed = dict(base, values=values)
            record = {'requests': len(requests), 'sample': sample, 'solution_dimension': fiber.variables - len(pivots),
                      'tensor_sha256': digest(values), 'complete_frame_profiles': profiles}
            if len(records) in (0, 4):
                record['fixed_gauge_cost_replay'] = cost_replay(seed, lattice, values)
            records.append(record)
    inconsistent = [fiber.jet_equation(0, 2, 0, 0), fiber.jet_equation(0, 2, 0, 1)]
    assert fiber.solve(inconsistent) is None
    names = (Path(__file__).name, 'finite_jet_smith_fibers.py', 'finite_precision_smith_frames.py',
             'full_alphabet_two_cut_lattices.py', 'search_finite_jet_smith_fibers.py',
             'check_bell_jet_crt_crossings.py', 'check_dense_coupled_reflection_frames.py',
             'search_fused_mixed_weighted.py', 'affine_place_profile_oracle.py')
    report = {'status': 'passed',
              'scope': 'Exact F4 equation reduction, complete fixed-frame replay, independent Hasse/endpoint controls, and complete-oracle pressure agreement; no global optimization claim.',
              'small_systems': systems, 'brute_force_assignments': assignments,
              'fiber_variables': fiber.variables, 'full_tensor_controls': records,
              'independent_Hasse_entries': finite_entries, 'independent_endpoint_entries': endpoint_entries,
              'inconsistent_target_rejected': True,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('small_systems', 'full_tensor_controls', 'source_sha256')}))


if __name__ == '__main__':
    main()
