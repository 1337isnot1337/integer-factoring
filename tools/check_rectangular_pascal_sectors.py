#!/usr/bin/env python3
"""Independent rectangular entries, endpoints, complete orders and derivatives."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_f4_pascal_coordinate_sectors import independent_entries, endpoint_expected
from f4_pascal_coordinate_sectors import F4SectorPascalLetterPressure
from f4_pascal_flags import change_flag
from permuted_pascal_kernels import permuted_pascal_kernel
from rectangular_pascal_sectors import CoordinateTensorBuilder, RectangularPascalPressure
from weighted_smith_oracle_f4 import WeightedSmithOracle
from windowed_f4_smith import pack


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng, controls = np.random.default_rng(0xABBA42), []
    entries, orders, derivatives, equal_orders = 0, 0, 0, 0
    for images, ab, ax in (([0,1],[0,3],[0,1,2,3]),
                           ([0,2,1,3],[0,1,3],[0,1,2,3,7]),
                           ([0,1,2,3],[0,1,3,7,15],[0,2,5])):
        original = permuted_pascal_kernel(images)
        flag = [[int(i == j) for j in range(len(images))] for i in range(len(images))]
        flag[0][1] = 2
        seed = change_flag(original, flag)
        builder = CoordinateTensorBuilder(seed)
        tensor = builder.tensor(ab, ax)
        independent = independent_entries(seed, sorted(set(ab+ax)))
        rectangular = {}
        for labels in product(ab, ax, ax, ab):
            pair = tensor['values'].get(labels, (0,0))
            assert pack(pair) == independent[labels]
            assert (pair[0].bit_count() % 2 | (pair[1].bit_count() % 2) << 1) == endpoint_expected(seed, labels)
            rectangular[labels] = independent[labels]
        entries += len(rectangular)
        model = RectangularPascalPressure(seed, ab, ax, builder=builder)
        parameters = model.natural(1.)
        parameters[:-1] += rng.integers(-3,4,len(parameters)-1)
        value = model.evaluate(parameters)
        ib, ix = {w:i for i,w in enumerate(ab)}, {w:i for i,w in enumerate(ax)}
        a, b = len(ab), len(ax)
        for place, (_, legs) in enumerate(model.cuts):
            matrix = [[0]*(a*b) for _ in range(a*b)]
            for (wb, wx, wv, wy), entry in rectangular.items():
                r = ib[wb]+a*ix[wv if place == 0 else wx]
                c = ix[wx if place == 0 else wv]+b*ib[wy]
                matrix[r][c] = entry
            h = model.heights(parameters)
            row = [int(h[place][legs[0]][i]+h[place][legs[1]][j]) for j in range(b) for i in range(a)]
            col = [int(h[place][legs[2]][i]+h[place][legs[3]][j]) for j in range(a) for i in range(b)]
            exact, _ = WeightedSmithOracle(matrix, infinity=bool(place)).profile(row, col)
            assert list(exact) == model.last['complete_original_profiles'][place]['orders']
            orders += len(exact)
        probe = model.natural(.035)
        probe[:-1] += rng.normal(size=len(probe)-1)*.0001
        _, analytic = model.evaluate(probe, True)
        errors = []
        for coordinate in sorted(set([0,a,a+b,model.size,2*model.size]+rng.integers(0,2*model.size,5).tolist())):
            delta = np.zeros_like(probe)
            delta[coordinate] = 1e-7
            finite = (model.evaluate(probe+delta)-model.evaluate(probe-delta))/(2e-7)
            errors.append(abs(float(finite-analytic[coordinate])))
            assert errors[-1] < 2e-5, (coordinate, errors[-1])
            derivatives += 1
        controls.append({'images':images, 'alphabet_b':ab, 'alphabet_x':ax,
                         'independent_orders_per_cut':a*b, 'maximum_gradient_error':max(errors),
                         'integer_gauge_rho':float(value), 'nonbinary_entries':model.construction['nonbinary_original_entries']})
        # Equal alphabets must reproduce the established implementation,
        # including all gradient components and the scalar-scale derivative.
        equal = sorted(set(ab+ax))[:5]
        previous = F4SectorPascalLetterPressure(seed, equal)
        current = RectangularPascalPressure(seed, equal, equal, builder=builder)
        probe = current.natural(.035)
        probe[:-1] += rng.normal(size=len(probe)-1)*.0001
        v0, g0 = previous.evaluate(probe, True)
        v1, g1 = current.evaluate(probe, True)
        assert abs(v0-v1) < 1e-12 and np.allclose(g0,g1,atol=1e-9,rtol=1e-9)
        for c0,c1 in zip(previous.last['complete_original_profiles'],current.last['complete_original_profiles']):
            assert np.allclose(c0['orders'],c1['orders'],atol=1e-12,rtol=1e-12)
            equal_orders += len(c0['orders'])
    root = Path(__file__).resolve().parents[1]
    names = ['tools/rectangular_pascal_sectors.py','tools/check_rectangular_pascal_sectors.py',
             'tools/check_f4_pascal_coordinate_sectors.py','tools/pascal_coordinate_sectors.py',
             'tools/f4_pascal_coordinate_sectors.py','tools/weighted_pascal_flags.py',
             'tools/native_matrix_smith_f4.py','tools/native_matrix_smith_f4.cpp',
             'tools/native_rectangular_smith_f4.cpp','tools/weighted_smith_oracle_f4.py']
    report = {'status':'pass', 'independent_entries_and_endpoints':entries,
              'independent_nonbinary_complete_orders':orders, 'finite_derivatives':derivatives,
              'equal_alphabet_regression_orders':equal_orders, 'controls':controls,
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope':'Independent paired coordinate alphabets with complete original Fock profiles. '
                      'These finite controls do not certify a stronger contraction or factoring bound.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('controls','source_sha256')}))


if __name__ == '__main__':
    main()
