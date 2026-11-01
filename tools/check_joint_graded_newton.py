#!/usr/bin/env python3
"""Energy-refined Newton ceiling for every rational joint-energy gauge."""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import reduce
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_joint_energy_conjugations import cut_matrix
from f4_lattice_linear import f16_multiply
from joint_energy_conjugations import model
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]
MUL = [[f16_multiply(a, b) for b in range(16)] for a in range(16)]
INV = [0]+[next(b for b in range(1, 16) if MUL[a][b] == 1) for a in range(1, 16)]
SCALARS = (14, 1, 6, 14, 6, 15, 5, 8, 3, 12, 13)
STRIDE = 1024


def scale(planes, scalar):
    output = [0]*4
    for i, value in enumerate(planes):
        for j in range(4):
            if MUL[1 << i][scalar] >> j & 1:
                output[j] ^= value
    return output


def characteristic_polynomial(matrix, blocks, diagonal):
    # det(I+z D R) in characteristic two: a subset permanent is a determinant.
    # Bit STRIDE*k+j is t^j z^k, in four independent F16 coefficient planes.
    result = [1, 0, 0, 0]
    for block in blocks:
        size = len(block)
        states = [[0]*4 for _ in range(1 << size)]
        states[0] = [1, 0, 0, 0]
        for used in range((1 << size)-1):
            r = used.bit_count()
            value = states[used]
            if not any(value):
                continue
            scaled = scale(value, diagonal[block[r]])
            for j, c in enumerate(block):
                if used >> j & 1:
                    continue
                term = value[:] if r == j else [0]*4
                entry = matrix[block[r]][c]
                while entry:
                    bit = entry & -entry
                    power = STRIDE+bit.bit_length()-1
                    for plane in range(4):
                        term[plane] ^= scaled[plane] << power
                    entry ^= bit
                target = states[used | (1 << j)]
                for plane in range(4):
                    target[plane] ^= term[plane]
        other = states[-1]
        support = other[0] | other[1] | other[2] | other[3]
        combined = [0]*4
        while support:
            bit = support & -support
            position = bit.bit_length()-1
            scalar = sum(((p >> position) & 1) << i for i, p in enumerate(other))
            scaled = scale(result, scalar)
            for plane in range(4):
                combined[plane] ^= scaled[plane] << position
            support ^= bit
        result = combined
    mask = (1 << STRIDE)-1
    return [[(plane >> (STRIDE*k)) & mask for plane in result] for k in range(26)]


def evaluate_binary(polynomial, point):
    value = 0
    for i in range(polynomial.bit_length()-1, -1, -1):
        value = MUL[value][point] ^ ((polynomial >> i) & 1)
    return value


def determinant(matrix):
    # Independent finite-field Gaussian elimination, without subset expansion.
    matrix = [row[:] for row in matrix]
    result = 1
    for c in range(len(matrix)):
        p = next((r for r in range(c, len(matrix)) if matrix[r][c]), None)
        if p is None:
            return 0
        matrix[c], matrix[p] = matrix[p], matrix[c]
        pivot = matrix[c][c]
        result = MUL[result][pivot]
        for r in range(c+1, len(matrix)):
            scalar = MUL[matrix[r][c]][INV[pivot]]
            for j in range(c+1, len(matrix)):
                matrix[r][j] ^= MUL[scalar][matrix[c][j]]
    return result


def certificate():
    data = model()
    assert all(not pair[1] for row in data['operator'] for pair in row)
    matrix = [[pair[0] << data['energies'][r] for pair in row]
              for r, row in enumerate(data['operator'])]
    assert 25*max(value.bit_length()-1 for row in matrix for value in row if value) < STRIDE
    reach = [[bool(value) or r == c for c, value in enumerate(row)]
             for r, row in enumerate(matrix)]
    for k, r, c in product(range(25), repeat=3):
        reach[r][c] |= reach[r][k] and reach[k][c]
    remaining, blocks = set(range(25)), []
    while remaining:
        r = min(remaining)
        block = sorted(c for c in remaining if reach[r][c] and reach[c][r])
        blocks.append(block)
        remaining.difference_update(block)
    assert sorted(map(len, blocks)) == [12, 13]
    diagonal = [0]*25
    for block, scalar in zip(data['blocks'], SCALARS):
        for r in block:
            diagonal[r] = scalar
    coefficients = characteristic_polynomial(matrix, blocks, diagonal)
    valuations = []
    for planes in coefficients:
        support = planes[0] | planes[1] | planes[2] | planes[3]
        assert support
        valuations.append((support & -support).bit_length()-1)
    assert valuations == [0, 8, 20, 28, 40, 56, 68, 88, 100, 120, 132, 152, 168,
                          188, 204, 224, 244, 264, 284, 304, 328, 348, 372, 392, 416, 444]
    for t in range(16):
        numeric = [[MUL[diagonal[r]][evaluate_binary(value, t)] for value in row]
                   for r, row in enumerate(matrix)]
        evaluated = []
        for planes in coefficients:
            value = 0
            for i, polynomial in enumerate(planes):
                value ^= MUL[1 << i][evaluate_binary(polynomial, t)]
            evaluated.append(value)
        for z in range(16):
            value = 0
            for scalar in reversed(evaluated):
                value = MUL[value][z] ^ scalar
            literal = [[MUL[z][entry] ^ int(r == c) for c, entry in enumerate(row)]
                       for r, row in enumerate(numeric)]
            assert determinant(literal) == value
    hull = []
    for i, value in enumerate(valuations):
        while len(hull) > 1:
            a, b = hull[-2:]
            if (b[1]-a[1])*(i-b[0]) < (value-b[1])*(b[0]-a[0]):
                break
            hull.pop()
        hull.append([i, value])
    slopes = []
    for a, b in zip(hull, hull[1:]):
        slope = Fraction(b[1]-a[1], b[0]-a[0])
        slopes.extend([[slope.numerator, slope.denominator]]*(b[0]-a[0]))
    zero, minors = WeightedSmithOracle(cut_matrix(data['values'], False)).profile(
        data['energies'], [0]*25)
    # The top characteristic coefficient agrees with a separate fraction-free determinant.
    determinant_bits = minors[-1].polynomial
    assert all(not ((determinant_bits >> (2*i+1)) & 1)
               for i in range((determinant_bits.bit_length()+1)//2))
    binary = sum(((determinant_bits >> (2*i)) & 1) << i
                 for i in range((determinant_bits.bit_length()+1)//2))
    assert binary
    weighted = scale([binary << sum(data['energies']), 0, 0, 0],
                     reduce(lambda a, b: MUL[a][b], diagonal, 1))
    assert weighted == coefficients[-1]
    names = ('check_joint_graded_newton.py', 'joint_energy_conjugations.py',
             'search_f4_recentered_crossings.py', 'f4_lattice_linear.py')
    return dict(energy_scalars=list(SCALARS), support_components=blocks,
                characteristic_coefficients=coefficients, coefficient_valuations=valuations,
                Newton_vertices=hull, Newton_slopes=slopes, original_zero_orders=list(zero),
                independent_F16_determinants=256, fraction_free_determinant_agrees=True,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in names},
                scope='Every characteristic-two field and every invertible rational gauge preserving the joint-energy blocks, without polynomiality or endpoint assumptions.',
                conclusion='Every zero Smith prefix lies below the listed energy-refined Newton polygon; total order is 444.',
                limit='The polygon need not be attainable, and no control of the opposite-cut infinity profile follows.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-graded-newton.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('source_sha256', 'characteristic_coefficients')}))


if __name__ == '__main__':
    main()
