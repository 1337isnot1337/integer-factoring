#!/usr/bin/env python3
"""Compare complete spectator-convolved small/full-core spectra at M8."""
from collections import Counter
import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from check_rectangular_defect_core_basis import data_for
from check_rectangular_defect_coordinates import coordinate_weights
from native_matrix_smith_f4 import MatrixOracle
from rectangular_defect_coordinates import coordinate_matrix


def integer_quarter(value):
    result = round(4*value)
    assert abs(result-4*value) < 1e-9
    return result


def spectrum(M, place, full_core, gauges):
    started = time.monotonic()
    data = data_for(place, M, full_core=full_core)
    matrix, shifts, metadata = coordinate_matrix(data)
    rw, cw = coordinate_weights(metadata, gauges)
    orders = MatrixOracle(matrix, len(matrix), shifts, precision=metadata['profile_precision']).profile(rw, cw)['orders']
    if not metadata['exact_polynomial_matrix']:
        known = metadata['known_below']+min(rw)+min(cw)
        assert max(orders) < known-1e-9
    histogram = Counter(integer_quarter(value) for value in orders)
    inactive = [i for i in range(2*M) if i not in data['active_indices']]
    for index, order in zip(inactive, data['spectator_orders']):
        sign = 1 if place == 'infinity' or index < M else -1
        shift = integer_quarter(order+sign*(gauges[0]+gauges[1]))
        next_histogram = histogram.copy()
        for value, count in histogram.items():
            next_histogram[value+shift] += count
        histogram = next_histogram
    assert sum(histogram.values()) == 3*(1 << (2*M-2))
    record = {'full_core': full_core, 'active_modes': len(data['active_indices']),
              'active_rank': len(matrix), 'spectator_modes': len(inactive),
              'full_rank': sum(histogram.values()), 'seconds': time.monotonic()-started,
              'profile_precision': metadata['profile_precision'],
              'exact_polynomial_matrix': metadata['exact_polynomial_matrix'],
              'active_orders_times_four': [integer_quarter(value) for value in orders]}
    return histogram, record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--M', type=int, default=8)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = []
    for place, charges in (('infinity', (0.75, -1.25)), ('zero', (-0.5, 1.5))):
        gauges = (*charges, [0, 0.5, -1, 1.5], [0, -0.75, 1.25], [0, 0.25, -0.5])
        small, small_record = spectrum(args.M, place, False, gauges)
        full, full_record = spectrum(args.M, place, True, gauges)
        assert small == full, (place, list((small-full).items())[:8], list((full-small).items())[:8])
        encoded = sorted(small.items())
        cases.append({'place': place, 'M': args.M,
                      'lambda': gauges[0], 'mu': gauges[1], 'v': gauges[2], 'y': gauges[3], 'u': gauges[4], 'x': [0, 0, 0, 0],
                      'all_orders_compared': sum(small.values()),
                      'spectrum_times_four': encoded,
                      'spectrum_sha256': sha256(json.dumps(encoded, separators=(',', ':')).encode()).hexdigest(),
                      'small_support': small_record, 'full_support': full_record})
        print(json.dumps({'place': place, 'orders_matched': sum(small.values())}), flush=True)
    result = {'status': 'complete weighted-spectrum agreement of independent active supports', 'cases': cases}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
