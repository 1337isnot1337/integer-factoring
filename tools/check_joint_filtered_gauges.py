#!/usr/bin/env python3
"""Opposite-filtration rigidity and exact entry-pole controls for q5 gauges."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json

from check_joint_endpoint_loops import operator_coefficients
from joint_energy_conjugations import IDENTITY, N, Space, binary_product, kernel, model, scale

ROOT = Path(__file__).resolve().parents[1]


def add(first, second):
    return first[0] ^ second[0], first[1] ^ second[1]


def transvection(source, r, c, scalar):
    value = [row[:] for row in source]
    value[r] = [add(a, scale(b, scalar)) for a, b in zip(value[r], value[c])]
    for i in range(N):
        value[i][c] = add(value[i][c], scale(value[i][r], scalar))
    return value


def certificate():
    data = model()
    heights = data['energies']
    c0 = data['weighted'][-8]
    square = binary_product(c0, c0)
    assert binary_product(square, square) == IDENTITY
    positions = [(r, c) for r in range(N) for c in range(N) if heights[r] >= heights[c]]
    energy_positions = [i for i, (r, c) in enumerate(positions) if heights[r] == heights[c]]
    assert len(positions) == 347 and len(energy_positions) == 69
    coefficients = operator_coefficients(data)
    equations = [[0]*(N*N) for _ in range(10)]
    for variable, (r, c) in enumerate(positions):
        elementary = 1 << (N*r+c)
        for degree, coefficient in enumerate(coefficients):
            value = binary_product(elementary, coefficient) ^ binary_product(coefficient, elementary)
            while value:
                bit = value & -value
                position = bit.bit_length()-1
                power = 8+heights[position % N]-degree
                if power < 0:
                    assert power >= -10
                    equations[power+10][position] ^= 1 << variable
                value ^= bit
    tangent = kernel([row for block in equations for row in block if row], len(positions))
    assert len(tangent) == 69
    assert len(Space(tangent+[1 << i for i in energy_positions]).basis) == 69
    gap_counts, controls = Counter(), 0
    for r, c in positions:
        if heights[r] == heights[c]:
            continue
        gap = heights[r]-heights[c]
        for scalar in (1, 2, 3):
            transformed = transvection(data['operator'], r, c, scalar)
            minimum = min(heights[j]-max(a.bit_length(), b.bit_length())+1
                          for row in transformed for j, (a, b) in enumerate(row) if a | b)
            assert minimum == -8-gap
            for i, row in enumerate(transformed):
                for j, (a, b) in enumerate(row):
                    assert ((a.bit_count() % 2) | ((b.bit_count() % 2) << 1)) == int(i == j)
            controls += 1
            gap_counts[gap] += 1
    assert controls == 834
    return dict(joint_operator_leading_fourth_power_identity=True,
                filtered_matrix_dimension=347, energy_block_dimension=69,
                no_new_pole_tangent_dimension=69, entry_pole_controls=controls,
                strict_transvection_gap_counts={str(k): v for k, v in sorted(gap_counts.items())},
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('check_joint_filtered_gauges.py', 'check_joint_endpoint_loops.py',
                                            'joint_energy_conjugations.py', 'search_f4_recentered_crossings.py')},
                scope='Every characteristic-two field and every invertible constant gauge preserving the stated energy filtration; proof uses the regular invertible joint operator before opposite-cut realignment.',
                conclusion='The whole zero Smith profile is preserved. The smallest infinity order is exactly -8-max(H_r-H_c: (U^-1)_rc nonzero). Retaining the original pole forces the gauge to preserve exact energy blocks.',
                limit='Gauges spending a worse pole budget, or preserving the zero profile without preserving this filtration, are not excluded.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-filtered-gauges.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}))


if __name__ == '__main__':
    main()
