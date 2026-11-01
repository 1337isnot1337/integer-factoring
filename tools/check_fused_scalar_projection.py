#!/usr/bin/env python3
"""Exact binary rank-span obstruction for a two-copy fused coefficient envelope."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
from itertools import product
import json
import os
from pathlib import Path
import shlex
import struct
import subprocess
from tempfile import TemporaryDirectory

from f4_lattice_linear import f16_multiply
from joint_energy_conjugations import Space

ROOT = Path(__file__).resolve().parents[1]
DUAL = 1572866
TARGET = 1118481
SPARSE = (
    ((0, 8), (12, 7), (24, 3), (28, 7), (30, 11), (40, 7), (42, 6),
     (52, 7), (56, 3), (58, 15), (60, 1), (68, 7), (70, 11), (80, 8)),
    ((0, 4), (10, 11), (12, 4), (22, 15), (24, 4), (28, 4), (30, 10),
     (38, 6), (40, 5), (42, 6), (50, 11), (52, 4), (56, 4), (58, 15),
     (68, 4), (70, 10), (80, 4)),
    ((0, 13), (10, 10), (12, 3), (20, 1), (22, 15), (24, 7), (28, 3),
     (38, 6), (40, 3), (50, 10), (52, 3), (56, 7), (68, 3), (80, 13)),
)


def coefficient_matrices():
    row, column = [], []
    for sparse in SPARSE:
        matrix = [[0]*9 for _ in range(9)]
        for position, value in sparse:
            matrix[position//9][position % 9] = value
        other = [[0]*9 for _ in range(9)]
        for x, y, u, v in product(range(3), repeat=4):
            other[x+3*u][v+3*y] = matrix[3*y+x][3*u+v]
        row.append(matrix)
        column.append(other)
    for side in (row, column):
        assert all(side[0][r][c] ^ side[1][r][c] ^ side[2][r][c] == int(r == c)
                   for r, c in product(range(9), repeat=2))
    return row, column


def projected_tensor(first, second, scalar):
    value = 0
    for r, s, c, t in product(range(9), repeat=4):
        entry = f16_multiply(scalar, f16_multiply(first[r][c], second[s][t]))
        value |= (entry & 1) << (81*(9*r+s)+9*c+t)
    return value


def coefficient_space(ordered=False):
    cuts = coefficient_matrices()
    rows, columns, labels = [], [], []
    for i in range(3):
        for j in range(i, 3):
            for scalar in (1, 2, 4, 8):
                for side, vectors in zip(cuts, (rows, columns)):
                    value = projected_tensor(side[i], side[j], scalar)
                    if i != j:
                        value ^= projected_tensor(side[j], side[i], scalar)
                    vectors.append(value)
                labels.append([i, j, scalar, 'symmetric'])
    if ordered:
        for i, j, scalar in product(range(3), range(3), (1, 2, 4, 8)):
            for side, vectors in zip(cuts, (rows, columns)):
                vectors.append(projected_tensor(side[i], side[j], scalar))
            labels.append([i, j, scalar, 'ordered'])
    space = Space(rows)
    indices = [rows.index(v) for v in space.basis]
    identity = sum(1 << (81*r+r) for r in range(81))
    assert space.binary_coordinates(identity) == TARGET
    values = ([rows[i] for i in indices], [columns[i] for i in indices])
    for side in values:
        assert combine(side, TARGET) == identity
    return values, [labels[i] for i in indices]


def combine(basis, coordinates):
    value = 0
    while coordinates:
        bit = coordinates & -coordinates
        value ^= basis[bit.bit_length()-1]
        coordinates ^= bit
    return value


def partitions():
    result = []
    for opposite in (False, True):
        groups = {}
        for r in range(81):
            a, b = divmod(r, 9)
            key = ((a % 3-a//3, b % 3-b//3) if opposite else
                   (a % 3+a//3, b % 3+b//3))
            groups.setdefault(key, []).append(r)
        result.append(sorted(groups.values(), key=lambda group: (-len(group), group)))
    assert [len(g) for g in result[0]] == [len(g) for g in result[1]]
    return result


def write_blocks(path, cuts, groups):
    dimension = len(cuts[0])
    with path.open('wb') as output:
        output.write(struct.pack('<IIQQ', dimension, len(groups[0]), TARGET, DUAL))
        output.write(struct.pack('<25I', *(len(g) for g in groups[0])))
        for basis, partition in zip(cuts, groups):
            for value in basis:
                packed = []
                for block in partition:
                    for r in block:
                        row = sum(((value >> (81*r+c)) & 1) << j
                                  for j, c in enumerate(block))
                        # No off-block entries: block ranks add exactly.
                        assert row.bit_count() == ((value >> (81*r)) & ((1 << 81)-1)).bit_count()
                        packed.append(row)
                output.write(struct.pack('<81H', *packed))


def dense_binary_rank(value):
    # Independent full 81 by 81 calculation, with highest-column row pivots.
    pivots = {}
    for r in range(81):
        row = (value >> (81*r)) & ((1 << 81)-1)
        while row:
            p = row.bit_length()-1
            if p in pivots:
                row ^= pivots[p]
            else:
                pivots[p] = row
                break
    return len(pivots)


def certificate():
    cuts, labels = coefficient_space()
    assert len(labels) == 21 and (DUAL & TARGET).bit_count() % 2 == 1
    groups = partitions()
    with TemporaryDirectory(prefix='factoring-fused-scalar-') as directory:
        directory = Path(directory)
        data, native = directory/'blocks.bin', directory/'screen'
        write_blocks(data, cuts, groups)
        source = ROOT/'tools/screen_fused_scalar_projection.cpp'
        compiler = shlex.split(os.environ.get('CXX', 'c++'))
        subprocess.run(compiler+['-O3', '-std=c++17', str(source), '-o', str(native)], check=True)
        reports, cheap = {}, []
        for mode in ('full', 'dual'):
            output = directory/f'{mode}.csv'
            run = subprocess.run([str(native), str(data), str(output), mode],
                                 check=True, text=True, capture_output=True)
            report = json.loads(run.stdout)
            assert report['dimension'] == 21 and report['target'] == TARGET
            assert report['states'] == 1 << (21-int(mode == 'dual'))
            reports[mode] = report
            records = [[int(row[k]) for k in ('coordinates', 'row_rank', 'column_rank')]
                       for row in csv.DictReader(output.open())]
            assert report['cheap'] == len(records)
            if mode == 'full':
                cheap = sorted(records)
            else:
                assert not records
    assert len(cheap) == 17
    assert len(Space([v for v, _, _ in cheap]).basis) == 6
    for coordinate, rr, rc in cheap:
        assert (DUAL & coordinate).bit_count() % 2 == 0 and rr*rc <= 501
        assert [dense_binary_rank(combine(basis, coordinate)) for basis in cuts] == [rr, rc]
    # sqrt(2) < 99/70 and 502^(1/2) > 9^(99/70) > 9^sqrt(2).
    assert 99**2 > 2*70**2 and 502**35 > 9**99
    names = ('check_fused_scalar_projection.py', 'screen_fused_scalar_projection.cpp',
             'f4_lattice_linear.py', 'joint_energy_conjugations.py')
    return dict(dimension=21, physical_alphabet=9, scalar_projection='low bit; chi(1)=1',
                target_coordinates=TARGET, dual_functional=DUAL, basis_labels=labels,
                block_sizes=[len(g) for g in groups[0]],
                rank_product_threshold=501, cheap_pieces=cheap, cheap_span_dimension=6,
                exhaustive_reports=reports, independent_full_matrix_ranks=34,
                threshold_integer_comparison=dict(lower_rank_product=502,
                                                  sqrt2_upper_numerator=99,
                                                  sqrt2_upper_denominator=70),
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in names},
                scope='Binary symmetric two-copy scalar-projected coefficient envelope of the literal q3 theta4 fused crossing; every vector and every binary linear decomposition in this space.',
                conclusion='Every crossing decomposition contains a piece with rank product at least 502, so its rank-sum cost exceeds 9^sqrt(2).',
                limit='The full ordered 33-dimensional envelope, other fields, higher powers and physical type restrictions are not settled by this certificate.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/fused-scalar-projection.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('source_sha256', 'basis_labels', 'cheap_pieces')}))


if __name__ == '__main__':
    main()
