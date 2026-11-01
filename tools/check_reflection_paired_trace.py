#!/usr/bin/env python3
"""Independent exact controls for the uniform reflection paired trace."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from reflection_dual_rail_gaussian import laurent_kernel, tensor
from reflection_laurent_physical_frames import add, scale
from search_reflection_physical_bases import change_basis, generate_basis_pair
from laurent_series_coupled_frames import construct


def times(a, b):
    result = set()
    for i in a:
        for j in b:
            exponent = i + j
            if exponent in result:
                result.remove(exponent)
            else:
                result.add(exponent)
    return result


def schur_controls(interval):
    kernel = [[set(entry) for entry in row] for row in laurent_kernel(interval)]
    m = interval // 2
    a, b = {0, -2 * interval}, {0, 2 * interval}
    ab = times(a, b)
    count = 0
    for first in (1, 2):
        for second in (1, 2):
            for i in range(m):
                for j in range(m):
                    row, col = first * m + i, second * m + j
                    value = times(ab, kernel[row][col])
                    for k in range(m):
                        value ^= times(b, times(kernel[row][k], kernel[m - 1 - k][col]))
                        value ^= times(a, times(kernel[row][3 * m + k], kernel[4 * m - 1 - k][col]))
                    assert value == (ab if first != second and i == j else set())
                    count += 1
    scalar = {0}
    for _ in range(m // 2):
        scalar = times(scalar, ab)
    scalar = {e + interval * interval // 2 for e in scalar}
    assert scalar == {0, interval * interval}
    for i in range(m):
        for j in range(m):
            assert kernel[i][3 * m + j] == ({0} if i == j else set())
    return count


def paired_trace(seed, interval):
    q = seed['q_x']
    expected = (1 | (1 << (interval * interval)), 0)
    contractions = []
    for pair in ((0, 3), (1, 2)):
        other = [leg for leg in range(4) if leg not in pair]
        values = []
        for i in range(q):
            for j in range(q):
                coordinates = [0] * 4
                coordinates[other[0]], coordinates[other[1]] = i, j
                value = (0, 0)
                for k in range(q):
                    coordinates[pair[0]] = coordinates[pair[1]] = k
                    index = sum(c * q ** leg for leg, c in enumerate(coordinates))
                    value = add(value, seed['values'][index])
                assert value == (expected if i == j else (0, 0))
                values.append(value)
        contractions.append(sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest())
    return {'alphabet': q, 'paired_entries': 2 * q * q,
            'contracted_tensor_entries': 2 * q ** 3, 'trace_sha256': contractions}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    schur = [{'interval': m, 'cleared_schur_entries': schur_controls(m)} for m in (4, 8, 16)]
    originals = [dict(interval=m, **paired_trace(tensor(m), m)) for m in (4, 8)]
    rng = random.Random(312412)
    base = tensor(4)
    transformed = []
    for index in (1, 7, 37, 117, 138, 171):
        label, left, right = generate_basis_pair(index, 4, rng, [])
        changed = change_basis(base, left, right)
        transformed.append(changed)
        paired_trace(changed, 4)
    for sample in range(8):
        coefficients = {}
        for pole in (1, 2):
            g = [[0] * 4 for _ in range(4)]
            for i in range(2):
                for j in range(2, 4):
                    g[i][j] = rng.randrange(4)
            coefficients[pole] = g
        changed = construct(base, coefficients)
        transformed.append(changed)
        paired_trace(changed, 4)
    for sample in range(16):
        values = list(base['values'])
        for changed in rng.sample(transformed, 4):
            coefficient = rng.randrange(1, 4)
            values = [add(value, scale(add(a, b), coefficient))
                      for value, a, b in zip(values, base['values'], changed['values'])]
        paired_trace(dict(base, values=values), 4)
    names = (Path(__file__).name, 'reflection_dual_rail_gaussian.py', 'laurent_series_coupled_frames.py',
             'reflection_laurent_physical_frames.py', 'search_reflection_physical_bases.py')
    report = {'status': 'passed', 'scope': 'Uniform proof with exact cleared-Schur, complete-tensor, and affine-frame controls.',
              'schur_controls': schur, 'original_tensors': originals,
              'constant_basis_frames': 6, 'independent_laurent_series_frames': 8, 'affine_combinations': 16,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names}}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ('status', 'scope', 'schur_controls', 'constant_basis_frames',
                                          'independent_laurent_series_frames', 'affine_combinations')}))


if __name__ == '__main__':
    main()
