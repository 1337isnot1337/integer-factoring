#!/usr/bin/env python3
"""Independent finite derivatives of the normalized root on a literal q12 chart."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from normalized_support_gap import normalized_gap
from q12_reflection_graphs import natural, tensor
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    c = [0]*192
    c[0], c[51], c[97], c[148] = 2, 3, 1, 2
    model = ScaledAsymmetricPressure(tensor(c))
    rng = np.random.default_rng(121220261101)
    coordinates = sorted({11*k+j for k in range(8) for j in (0, 10)} | {88, 89})
    records, errors = [], []
    for y in (-1.2, 0., .8):
        point = np.r_[natural()[:-1]+rng.normal(0, .08, 88), math.log(.06), y]
        rate, gradient, _ = model.evaluate(point)
        value, derivative = normalized_gap(rate, point, gradient)
        a, b = math.cos(math.atan(math.exp(y))), math.sin(math.atan(math.exp(y)))
        assert abs(value-(rate-1)/(a*b)) < 1e-12
        for coordinate in coordinates:
            delta = np.zeros(90)
            delta[coordinate] = 1e-6
            plus = normalized_gap(model.evaluate(point+delta, False)[0], point+delta)
            minus = normalized_gap(model.evaluate(point-delta, False)[0], point-delta)
            error = abs(float((plus-minus)/2e-6-derivative[coordinate]))
            assert error < 2e-5, (y, coordinate, error)
            errors.append(error)
        assert normalized_gap(1.-1e-4, point) < 0 < normalized_gap(1.+1e-4, point)
        records.append({'point': point.tolist(), 'rate': rate, 'objective': value,
                        'gradient': derivative.tolist()})
    names = ['normalized_support_gap.py', 'check_normalized_support_gap.py',
             'search_q12_reflection_graphs.py', 'q12_reflection_graphs.py',
             'scaled_asymmetric_crossing_pressure.py', 'asymmetric_crossing_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp']
    report = {'status': 'pass', 'finite_derivative_comparisons': len(errors),
              'maximum_error': max(errors), 'records': records,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'scope': 'Numerical objective derivatives and threshold-sign controls on an actual nonbinary '
                       'q12 tensor. No exact support or new factoring bound.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'finite_derivative_comparisons', 'maximum_error')}))


if __name__ == '__main__':
    main()
