#!/usr/bin/env python3
"""Independent coefficient controls for finite-precision F4 Smith frames."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from affine_place_profile_oracle import AffinePlaceOracle, bitplanes
from bell_jet_crt_crossings import construct
from check_dense_coupled_reflection_frames import coefficients, multiply
from f16_identity_operator_profiles import inverse_matrix
from f4_regular_smith_frames import regular_frame, bounded_frame
from physical_split_crossings import identity_split
from reflection_laurent_physical_frames import digest
from search_entrywise_reflection_jets import bell_values


def add(first, second):
    result = dict(first)
    for exponent, coefficient in second.items():
        value = result.get(exponent, 0) ^ coefficient
        if value:
            result[exponent] = value
        else:
            result.pop(exponent, None)
    return result


def product(first, second, precision=None):
    result = {}
    for a, c in first.items():
        for b, d in second.items():
            exponent = a + b
            if precision is None or exponent < precision:
                result = add(result, {exponent: multiply(c, d)})
    return result


def matrix_product(first, second, precision=None):
    result = [[{} for _ in second[0]] for _ in first]
    for i in range(len(first)):
        for j in range(len(second[0])):
            for k in range(len(second)):
                result[i][j] = add(result[i][j], product(first[i][k], second[k][j], precision))
    return result


def replay(source, frame, precision):
    orders, left, right, diagonal = frame
    first = [[coefficients(value) for value in row] for row in left]
    second = [[coefficients(value) for value in row] for row in right]
    for matrix in (first, second):
        constant = [[value.get(0, 0) for value in row] for row in matrix]
        inverse_matrix(constant)
    observed = matrix_product(matrix_product(first, source, precision), second, precision)
    assert observed == [[coefficients(value) for value in row] for row in diagonal]
    assert all(not observed[i][j] for i in range(len(source)) for j in range(len(source)) if i != j)
    assert tuple(min(observed[i][i]) for i in range(len(source))) == orders
    return len(source) ** 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(312511)
    records, compared = [], 0
    for trial in range(24):
        size = 1 + trial % 5
        expected = sorted(rng.randrange(13) for _ in range(size))
        units = []
        for _ in range(2):
            matrix = []
            for i in range(size):
                row = []
                for j in range(size):
                    value = {0: 1} if i == j else {}
                    for exponent in range(1, 4):
                        coefficient = rng.randrange(4)
                        if coefficient:
                            value[exponent] = coefficient
                    row.append(value)
                matrix.append(row)
            units.append(matrix)
        diagonal = [[{expected[i]: rng.randrange(1, 4)} if i == j else {}
                     for j in range(size)] for i in range(size)]
        source = matrix_product(matrix_product(units[0], diagonal), units[1])
        packed = [[bitplanes(value) for value in row] for row in source]
        precision = max(expected) + 4
        frame = regular_frame(packed, precision)
        assert list(frame[0]) == expected
        compared += replay(source, frame, precision)
        assert bounded_frame(packed, 4)[0] == frame[0]
        records.append({'size': size, 'orders': expected, 'precision': precision})
    identity = identity_split(4)
    alternating = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
    crt_records = []
    for label, matrix in (('identity', identity), ('alternating', alternating)):
        carrier = bell_values(identity, 1)
        for entry in range(256):
            x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
            carrier[entry] ^= 3 * int(x == v and y == u)
        seed = construct(matrix, identity, carrier)
        oracle = AffinePlaceOracle(seed)
        for cut, point in ((0, (2, 1)), (1, (3, 1))):
            source = oracle.local_matrix(cut, point)
            packed = [[bitplanes(value) for value in row] for row in source]
            frame = regular_frame(packed, 40)
            assert frame[0] == oracle.profile(cut, point, [0] * 16)
            compared += replay(source, frame, 40)
            crt_records.append({'relative_form': label, 'cut': cut, 'point': list(point),
                                'tensor_sha256': digest(seed['values']), 'orders': list(frame[0])})
    insufficient = False
    try:
        regular_frame([[(1 << 15, 0)]], 15)
    except ValueError:
        insufficient = True
    assert insufficient
    names = (Path(__file__).name, 'f4_regular_smith_frames.py',
             'check_dense_coupled_reflection_frames.py', 'affine_place_profile_oracle.py',
             'bell_jet_crt_crossings.py', 'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')
    report = {'status': 'passed',
              'scope': 'Exact local regular frames and complete Smith orders to sufficient precision; no tensor pressure claim.',
              'constructed_profile_controls': len(records), 'independent_coefficient_matrix_entries': compared,
              'nonsingular_constant_frame_checks': 2 * (len(records) + len(crt_records)),
              'complete_CRT_cut_controls': crt_records, 'random_controls': records,
              'insufficient_precision_rejected': insufficient,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('complete_CRT_cut_controls', 'random_controls', 'source_sha256')}))


if __name__ == '__main__':
    main()
