#!/usr/bin/env python3
"""Replay comparable extreme gauges behind the floating Gibbs-root guard.

The historical failures did not retain their parameters. This diagnostic
does not reconstruct those events. It constructs separate comparable
gauges within the same bounds, certifies complete integer profiles, and
brackets their true Gibbs ratios with outward integer interval arithmetic.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import time

from charged_gaussian_floor import Dyadic
from check_closed_auxiliary_blocks import literal_entry
from check_f4_recentered_walks import certify_profile
from closed_auxiliary_blocks import tensor
from f4_recentered_lattices import pack
from native_rectangular_smith_f4 import RectangularOracle
from paired_gibbs_interval import pressure
from search_rectangular_f4_pressure import logsumexp, paired_value_gradient


def local_matrix(seed, heights, side, cut):
    q, degree = seed['dimensions'][0], seed['degree']
    n = q*q
    clear = max(0, -4*min(heights))
    result = [[0]*n for _ in range(n)]
    actual_cut = cut if side == 0 else 1-cut
    for index, pair in enumerate(seed['values']):
        x, y, u, v = ((index//q**leg) % q for leg in range(4))
        if actual_cut == 0:
            row, column = x+q*y, u+q*v
        else:
            row, column = x+q*u, y+q*v
        weight = sum(heights[label] for label in (x, y, u, v))+clear
        assert weight >= 0
        if side == 0:
            pair = tuple(int(format(value, '0'+str(degree+1)+'b')[::-1], 2) for value in pair)
        result[row][column] = pack(pair) << (2*weight)
    return result, clear+(degree if side == 0 else 0)


def high_precision_root(profiles, energies, scale):
    from decimal import Decimal, localcontext
    with localcontext() as context:
        context.prec = 200
        scaled = [[Decimal(scale)*value for value in row] for row in profiles]
        gx = [Decimal(scale)*value for value in energies]

        def log_partition(values):
            minimum = min(values)
            counts = Counter(values)
            tail = Decimal(counts.pop(minimum)-1)
            tail += sum(count*(-(value-minimum)).exp() for value, count in counts.items())
            return -minimum+(1+tail).ln()

        minimum = min(gx)
        shifted = [value-minimum for value in gx]
        constant = sum(log_partition(row) for row in scaled)+2*minimum

        def equation(rate):
            return constant-2*rate*log_partition([value/rate for value in shifted])

        low, high = Decimal('.01'), Decimal(8)
        assert equation(low) > 0 > equation(high)
        for _ in range(170):
            middle = (low+high)/2
            if equation(middle) > 0:
                low = middle
            else:
                high = middle
        rho = (low+high)/4
        return rho, constant


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--historical', type=Path, default=Path('evidence/closed-auxiliary-eight-symbol-search.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    historical = json.loads(args.historical.read_text())
    data = tensor(3, 0, 0)
    q = 8
    for index, (value, high) in enumerate(data['values']):
        labels = tuple((index//q**leg) % q for leg in range(4))
        assert high == 0 and value == literal_entry(3, labels, 0, 0)
        assert value.bit_count() % 2 == int(labels[0] == labels[3] and labels[1] == labels[2])
    oracle = RectangularOracle(data['values'], q, q)
    scale = 8
    cases = []
    exact_profiles = 0
    arithmetic = Dyadic(768)
    for height in (3, 6, 12):
        for kind in ('positive-nonvacuum', 'negative-preferred-one'):
            h = (0,)+(height,)*7 if kind == 'positive-nonvacuum' else (0, -height)+(0,)*6
            energies = [4*value for value in h]
            parameters = [scale*value for _ in range(8) for value in h[1:]]+[float(scale)]
            for cut in (0, 1):
                _, _, orders = oracle.evaluate(parameters, [1/q]*q, [1/q]*q, cut, True)
                assert all(abs(value/scale-round(value/scale)) < 1e-10 for value in orders)
                profiles = [tuple(round(value/scale) for value in orders[:q*q]),
                            tuple(round(value/scale) for value in orders[q*q:])]
                for side in (0, 1):
                    matrix, scalar = local_matrix(data, h, side, cut)
                    certify_profile(matrix, scalar, profiles[side])
                    exact_profiles += 1
                try:
                    actual = paired_value_gradient(oracle, parameters, cut)
                    error = None
                    floating_rho = actual[0]
                except (ArithmeticError, ValueError, OverflowError) as exception:
                    error = type(exception).__name__+': '+str(exception)
                    floating_rho = None
                gx = [scale*value for value in energies]
                constant = 2*min(gx)+logsumexp([-value for value in orders[:q*q]])+logsumexp([-value for value in orders[q*q:]])
                rho, exact_constant_numeric = high_precision_root(profiles, energies, scale)
                integer = int(rho*10**6)
                low, high = Fraction(integer-1,10**6), Fraction(integer+2,10**6)
                lo_pressure = pressure(arithmetic, profiles, energies, energies, arithmetic.rational(scale), low)
                hi_pressure = pressure(arithmetic, profiles, energies, energies, arithmetic.rational(scale), high)
                assert lo_pressure[0] > 0 > hi_pressure[1]
                assert low > Fraction(7,10)
                cases.append({'kind': kind, 'normalized_heights_each_leg': h, 'scale': scale,
                              'cut': cut, 'parameters': parameters, 'complete_unscaled_profiles': profiles,
                              'normalized_physical_energies': energies,
                              'floating_guard_result': error, 'floating_rho_if_returned': floating_rho,
                              'floating_zero_rate_constant': constant,
                              'high_precision_zero_rate_constant_diagnostic': str(exact_constant_numeric),
                              'high_precision_rho_diagnostic': str(rho),
                              'certified_rho_bracket': [str(low), str(high)],
                              'lower_rho_pressure_interval_numerators': list(map(str,lo_pressure)),
                              'upper_rho_pressure_interval_numerators': list(map(str,hi_pressure))})
                print(json.dumps({'event': 'boundary_case', 'kind': kind, 'height': height,
                                  'cut': cut, 'guard': error, 'floating_constant': constant,
                                  'rho_bracket': [str(low),str(high)]}), flush=True)
    report = {'status': 'passed',
              'scope': 'Separate comparable extreme-gauge examples with exact profile certificates and positive true Gibbs-root brackets; the historical rejected parameter vectors were not saved and are not reconstructed.',
              'historical_record': str(args.historical),
              'historical_record_sha256': sha256(args.historical.read_bytes()).hexdigest(),
              'historical_native_profile_queries': historical['native_profile_queries'],
              'historical_rejected_evaluations': historical['rejected_evaluations'],
              'historical_failures_individually_classified': False,
              'source_seed': {'block_length': 3, 'mask_b': 0, 'mask_c': 0, 'dimensions': [8]*4},
              'literal_auxiliary_path_and_endpoint_entries': 4096,
              'exact_complete_matrix_profiles': exact_profiles,
              'dyadic_denominator': str(arithmetic.scale),
              'observed_diagnostic_guard_results': dict(Counter(row['floating_guard_result'] or 'returned' for row in cases)),
              'cases': cases,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
                  (Path(__file__).name, 'closed_auxiliary_blocks.py', 'native_rectangular_smith_f4.py',
                   'native_rectangular_smith_f4.cpp', 'search_rectangular_f4_pressure.py',
                   'paired_gibbs_interval.py', 'charged_gaussian_floor.py')},
              'elapsed_seconds': time.monotonic()-started}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'cases': len(cases), 'exact_profiles': exact_profiles,
                      'guard_results': report['observed_diagnostic_guard_results'],
                      'seconds': report['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
