#!/usr/bin/env python3
"""Controls for the separate stable Gibbs evaluator and compatible adapter."""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import time

from balanced_particle_extensions import REFERENCE_CODE, construct
from closed_auxiliary_blocks import tensor
from native_rectangular_smith_f4 import RectangularOracle
from search_rectangular_f4_pressure import paired_value_gradient
from stable_paired_gibbs import paired_value_gradient_with_fallback, stable_paired_result


def decimal_root(orders, energies):
    """Independent200-digit root of complete profile data; numerical control only."""
    with localcontext() as context:
        context.prec = 200
        profiles = [[Decimal.from_float(float(value)) for value in row] for row in orders]
        physical = [[Decimal.from_float(float(value)) for value in row] for row in energies]

        def log_partition(values):
            minimum = min(values)
            first = values.index(minimum)
            tail = sum((-(value-minimum)).exp() for index, value in enumerate(values) if index != first)
            return -minimum+(1+tail).ln()

        constant = sum(log_partition(row) for row in profiles)+sum(min(row) for row in physical)
        physical = [[value-min(row) for value in row] for row in physical]

        def equation(rate):
            return constant-rate*sum(log_partition([value/rate for value in row]) for row in physical)

        low, high = Decimal('.01'), Decimal(8)
        assert equation(low) > 0 > equation(high)
        for _ in range(160):
            middle = (low+high)/2
            if equation(middle) > 0:
                low = middle
            else:
                high = middle
        return float((low+high)/4)


def raw_data(oracle, parameters, cut):
    _, _, flat = oracle.evaluate(parameters, [1/oracle.A]*oracle.A,
                                  [1/oracle.B]*oracle.B, cut, True)
    energies = [[0.]*oracle.A, [0.]*oracle.B]
    for side in (0, 1):
        for leg in range(4):
            energy = energies[0 if leg in (0, 3) else 1]
            for label in range(1, len(energy)):
                energy[label] += parameters[side*oracle.stride+oracle.offsets[leg]+label-1]
    return (flat[:oracle.n],flat[oracle.n:]), energies


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--boundary', type=Path, default=Path('evidence/paired-gibbs-boundary-diagnostics.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    boundary = json.loads(args.boundary.read_text())
    assert boundary['status'] == 'passed'
    base = tensor(3, 0, 0)
    oracle = RectangularOracle(base['values'], 8, 8)
    recovered = 0
    boundary_rows = []
    directional_errors = []
    for item in boundary['cases']:
        result = stable_paired_result(oracle, item['parameters'], item['cut'])
        low, high = map(Fraction, item['certified_rho_bracket'])
        assert float(low) < result['rho'] < float(high)
        assert abs(result['rho']-float(item['high_precision_rho_diagnostic'])) < 2e-12
        assert result['gradient'] is not None and all(math.isfinite(value) for value in result['gradient'])
        compatible = paired_value_gradient_with_fallback(oracle, item['parameters'], item['cut'])
        if item['floating_guard_result']:
            recovered += 1
        assert abs(compatible[0]-result['rho']) < 1e-12
        if abs(item['normalized_heights_each_leg'][1]) >= 6:
            h = 1/1024
            sign = 1 if item['kind'] == 'positive-nonvacuum' else -1
            direction = [0.]*oracle.parameters
            for side in (0, 1):
                for leg in range(4):
                    for label in range(1, 8):
                        if item['kind'] == 'positive-nonvacuum' or label == 1:
                            direction[side*oracle.stride+oracle.offsets[leg]+label-1] = sign*item['scale']
            plus = [value+h*delta for value,delta in zip(item['parameters'],direction)]
            minus = [value-h*delta for value,delta in zip(item['parameters'],direction)]
            p = decimal_root(*raw_data(oracle,plus,item['cut']))
            m = decimal_root(*raw_data(oracle,minus,item['cut']))
            actual = (p-m)/(2*h)
            predicted = math.fsum(value*delta for value,delta in zip(result['gradient'],direction))
            directional_errors.append(abs(actual-predicted))
            assert abs(actual-predicted) < 2e-7, (actual,predicted)
        boundary_rows.append({'kind': item['kind'], 'normalized_height': item['normalized_heights_each_leg'][1],
                              'cut': item['cut'], 'rho': result['rho'],
                              'log_paired_entropy': result['log_paired_entropy'],
                              'log_zero_rate_constant': result['log_zero_rate_constant'],
                              'delta_numerator': result['delta_numerator'],
                              'gradient_status': result['gradient_status'],
                              'affine_consistency': result['exact_affine_consistency']})
    rng = random.Random(2026110107)
    ordinary = root_errors = gradient_errors = finite_differences = 0
    maximum_root_error = maximum_gradient_error = maximum_finite_error = 0.
    nonsmooth_directions = 0
    seeds = [construct(REFERENCE_CODE), construct(REFERENCE_CODE,REFERENCE_CODE+(85,)),
             tensor(3,0x318c,0xc33c)]
    for seed in seeds:
        a,b = seed['dimensions'][:2]
        oracle = RectangularOracle(seed['values'],a,b)
        for trial in range(8):
            parameters = [rng.randrange(-192,193)/256 for _ in range(oracle.parameters-1)]+[rng.randrange(32,129)/256]
            for cut in (0,1):
                stable = stable_paired_result(oracle,parameters,cut)
                ordinary_value,ordinary_gradient,px,py = paired_value_gradient(oracle,parameters,cut)
                error = abs(stable['rho']-ordinary_value)
                maximum_root_error = max(maximum_root_error,error)
                assert error < 2e-10
                error = max(abs(first-second) for first,second in zip(stable['gradient'],ordinary_gradient))
                maximum_gradient_error = max(maximum_gradient_error,error)
                assert error < 2e-8
                ordinary += 1
                coordinate = (7*trial+cut) % oracle.parameters
                step = 1/2**18
                plus,minus = parameters[:],parameters[:]
                plus[coordinate] += step
                minus[coordinate] -= step
                base_value = oracle.types(parameters,cut)[0]
                p = oracle.types(plus,cut)[0]
                m = oracle.types(minus,cut)[0]
                forward,backward = (p-base_value)/step,(base_value-m)/step
                if abs(forward-backward) > 3e-5:
                    nonsmooth_directions += 1
                    continue
                error = abs((p-m)/(2*step)-stable['gradient'][coordinate])
                maximum_finite_error = max(maximum_finite_error,error)
                assert error < 3e-6, (error,coordinate,parameters)
                finite_differences += 1
    assert ordinary == 48 and finite_differences >= 24
    report = {'status': 'passed',
              'scope': 'Exact affine consistency checks, stable floating roots within independently certified boundary brackets, and finite numerical gradient controls. Native global profile attainment is not proved by the affine reconstruction.',
              'boundary_record': str(args.boundary), 'boundary_record_sha256': sha256(args.boundary.read_bytes()).hexdigest(),
              'certified_boundary_brackets_matched': len(boundary_rows), 'diagnostic_guard_examples_recovered': recovered,
              'successful_low_entropy_diagnostics_replayed': len(boundary_rows)-recovered,
              'boundary_directional_decimal_controls': len(directional_errors),
              'maximum_boundary_directional_error': max(directional_errors),
              'ordinary_root_and_gradient_controls': ordinary,
              'maximum_ordinary_root_error': maximum_root_error,
              'maximum_ordinary_gradient_error': maximum_gradient_error,
              'ordinary_smooth_finite_difference_controls': finite_differences,
              'ordinary_nonsmooth_directions_skipped': nonsmooth_directions,
              'maximum_ordinary_finite_difference_error': maximum_finite_error,
              'boundary_cases': boundary_rows,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
                  (Path(__file__).name,'stable_paired_gibbs.py','native_rectangular_smith_f4.py',
                   'native_rectangular_smith_f4.cpp','search_rectangular_f4_pressure.py')},
              'elapsed_seconds': time.monotonic()-started}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:value for key,value in report.items() if key != 'boundary_cases'}),flush=True)


if __name__ == '__main__':
    main()
