#!/usr/bin/env python3
"""Exact tangent-space certificate for formal joint leading-rank rigidity."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from joint_energy_conjugations import N, binary_product, kernel, model

ROOT = Path(__file__).resolve().parents[1]


def realign(matrix):
    result = 0
    for x in range(5):
        for y in range(5):
            for u in range(5):
                for v in range(5):
                    if matrix >> (N*(x+5*y)+v+5*u) & 1:
                        result ^= 1 << (N*(x+5*u)+y+5*v)
    return result


def certificate():
    data = model()
    leading, positions = data['weighted'][-8], data['positions']
    matrix = realign(leading)
    rows = [(matrix >> (N*i)) & ((1 << N)-1) for i in range(N)]
    assert rows[0] & 1 and binary_rank(rows) == 1
    a, b = [row & 1 for row in rows], rows[0]
    assert all(row == (b if bit else 0) for row, bit in zip(rows, a))
    equations = [0]*(N*N)
    for variable, (i, j) in enumerate(positions):
        elementary = 1 << (N*i+j)
        derivative = realign(binary_product(elementary, leading) ^ binary_product(leading, elementary))
        delta = [(derivative >> (N*r)) & ((1 << N)-1) for r in range(N)]
        for r in range(N):
            for c in range(N):
                value = (((delta[r] >> c) & 1)
                         ^ (((b >> c) & 1)*(delta[r] & 1))
                         ^ (a[r]*((delta[0] >> c) & 1))
                         ^ (a[r]*((b >> c) & 1)*(delta[0] & 1)))
                if value:
                    equations[N*r+c] ^= 1 << variable
    tangent = kernel(equations, len(positions))
    physical = [sum(1 << k for k, (i, j) in enumerate(positions)
                    if i == j and (i % 5 if side == 0 else i//5) == letter)
                for side in (0, 1) for letter in range(5)]
    centralizer = [sum(((value >> (N*i+j)) & 1) << k for k, (i, j) in enumerate(positions))
                   for value in data['leading'].basis]
    assert len(tangent) == 16
    assert binary_rank(physical) == 9 and binary_rank(centralizer) == 10
    assert binary_rank(physical+centralizer) == 16
    assert binary_rank(tangent+physical+centralizer) == 16
    assert all(not ((row & vector).bit_count() % 2)
               for vector in physical+centralizer for row in equations)
    return dict(energy_block_dimension=69, tangent_equation_rank=53, tangent_dimension=16,
                physical_diagonal_dimension=9, leading_commutant_dimension=10,
                sum_dimension=16, intersection_dimension=3,
                tangent_basis=tangent, physical_generators=physical,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('check_joint_leading_tangent.py', 'joint_energy_conjugations.py',
                                            'search_f4_recentered_crossings.py')},
                scope='All characteristic-two coefficient fields; tangent identity at the reference, with formal lifting proved in the note. No classification of disconnected constant solutions.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-leading-tangent.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('tangent_basis', 'physical_generators', 'source_sha256')}))


if __name__ == '__main__':
    main()
