#!/usr/bin/env python3
"""Exact all-field pressure obstruction for the complete two-place pencil."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from math import ceil, comb
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from check_joint_endpoint_loops import integer_product, operator_coefficients, regular_matrices
from check_joint_energy_conjugations import polynomial
from check_joint_graded_newton import MUL, determinant, evaluate_binary
from check_joint_zero_prefix import two_place_pencil
from joint_energy_conjugations import N, binary_product, model
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]
STRIDE = 512


def symbolic_minors(source, change):
    # Each coefficient of any determinant has local-variable degree < STRIDE.
    maximum = max(value.bit_length()-1 for matrix in (source, change)
                  for row in matrix for value in row if value)
    assert N*maximum < STRIDE
    encoded = [[polynomial((a ^ (b << STRIDE), 0)) for a, b in zip(row, other)]
               for row, other in zip(source, change)]
    orders, minors = WeightedSmithOracle(encoded).profile([0]*N, [0]*N)
    records = []
    prefix = 0
    for order, minor in zip(orders, minors):
        prefix += order
        value = minor.polynomial
        assert not any((value >> (2*i+1)) & 1 for i in range((value.bit_length()+1)//2))
        binary = sum(((value >> (2*i)) & 1) << i for i in range((value.bit_length()+1)//2))
        coefficients = []
        while binary:
            coefficients.append(binary & ((1 << STRIDE)-1))
            binary >>= STRIDE
        valuations = [(a & -a).bit_length()-1 if a else None for a in coefficients]
        assert valuations[0] == prefix
        records.append(dict(size=minor.size, rows=minor.row_mask, columns=minor.column_mask,
                            coefficients=coefficients, valuations=valuations))
    return list(orders), records


def minor_controls(source, change, records):
    controls = 0
    for t in (2, 3, 4, 8):
        base = [[evaluate_binary(value, t) for value in row] for row in source]
        delta = [[evaluate_binary(value, t) for value in row] for row in change]
        for parameter in (1, 2, 4, 8):
            matrix = [[a ^ MUL[parameter][b] for a, b in zip(row, other)]
                      for row, other in zip(base, delta)]
            for record in records:
                rows = [i for i in range(N) if record['rows'] >> i & 1]
                columns = [i for i in range(N) if record['columns'] >> i & 1]
                value = 0
                for coefficient in reversed(record['coefficients']):
                    value = MUL[value][parameter] ^ evaluate_binary(coefficient, t)
                assert determinant([[matrix[r][c] for c in columns] for r in rows]) == value
                controls += 1
    return controls


def lower_hull(points):
    caps = {}
    for k, value in points:
        caps[k] = min(value, caps.get(k, value))
    vertices = []
    for k, value in sorted(caps.items()):
        while len(vertices) > 1:
            a, b = vertices[-2:]
            if (b[1]-a[1])*(k-b[0]) < (value-b[1])*(b[0]-a[0]):
                break
            vertices.pop()
        vertices.append([k, value])
    slopes = []
    for a, b in zip(vertices, vertices[1:]):
        slopes.extend([Fraction(b[1]-a[1], b[0]-a[0])]*(b[0]-a[0]))
    return vertices, slopes


def pressure_certificate(zero, infinity):
    # Rounding slopes up weakens the partition lower bound and removes fractions.
    integer_orders = [[ceil(value) for value in orders] for orders in (zero, infinity)]
    factors = []
    for orders in integer_orders:
        factor = [0]*(8*max(orders)+1)
        for order in orders:
            factor[8*order] += 1
        factors.append(factor)
    product = integer_product(*factors)
    physical = [0]*56
    for exponent in (0, 11, 33, 44, 55):
        physical[exponent] = 1
    square = integer_product(physical, physical)
    cube = integer_product(square, physical)
    difference = [0]*max(len(product), len(cube))
    for i, value in enumerate(product):
        difference[i] += 11*value
    for factor, multiplier in ((square, 1), (cube, 10)):
        for i, value in enumerate(factor):
            difference[i] -= multiplier*value
    while not difference[-1]:
        difference.pop()
    order = next(i for i, value in enumerate(difference) if value)
    assert order == 8
    quotient = difference[order:]
    certificate = None
    for elevation in (0, 16, 32, 64):
        degree = len(quotient)-1+elevation
        coefficients = [sum(value*comb(degree-i, k-i)
                            for i, value in enumerate(quotient[:k+1]) if value)
                        for k in range(degree+1)]
        if min(coefficients) > 0:
            certificate = dict(degree=degree, minimum=min(coefficients),
                               positive_coefficients=len(coefficients),
                               coefficient_sha256=sha256(json.dumps(coefficients).encode()).hexdigest())
            break
    assert certificate is not None and certificate['minimum'] == 22
    return dict(zero_integer_orders=integer_orders[0], infinity_integer_orders=integer_orders[1],
                difference_polynomial=difference, zero_factor=order, Bernstein=certificate)


def certificate():
    data = model()
    pencil = two_place_pencil(data, operator_coefficients(data))
    generator = pencil['generator']
    regular = [0]*29
    zero_source = []
    for r, row in enumerate(data['operator']):
        output = []
        for c, (low, high) in enumerate(row):
            assert not high
            value = low << data['energies'][r]
            assert not (value & 255)
            output.append(value >> 8)
            while value:
                bit = value & -value
                regular[bit.bit_length()-1] ^= 1 << (N*r+c)
                value ^= bit
        zero_source.append(output)
    changes = [binary_product(generator, value) ^ binary_product(value, generator)
               for value in regular]
    assert all(not value for i, value in enumerate(changes) if i not in (16, 20))
    assert changes[16] == changes[20]
    zero_change = [[(changes[16] >> (N*r+c)) & 1 for c in range(N)] for r in range(N)]
    assert binary_rank([sum(value << c for c, value in enumerate(row)) for row in zero_change]) == 2
    zero_orders, zero_minors = symbolic_minors(zero_source, zero_change)
    assert all(len(record['coefficients']) <= 3 for record in zero_minors)
    # Exactly these nonconstant coefficients need valuation/cancellation control.
    active = {record['size']: record['valuations'][2] for record in zero_minors
              if len(record['valuations']) > 2 and record['valuations'][2] is not None}
    assert active == {9: 26, 10: 35, 11: 44, 12: 54, 13: 64, 14: 75, 15: 86, 16: 98, 17: 116}
    assert all(len(record['valuations']) < 2 or record['valuations'][1] is None for record in zero_minors)
    # At ord(q)=5 the 17-prefix tie is bounded by convexity at 16 and 18.
    prefixes = [0]
    for order in zero_orders:
        prefixes.append(prefixes[-1]+order)
    assert prefixes[16]-16+2*5+prefixes[18] <= 2*prefixes[17]
    infinity_source, infinity_change = regular_matrices(data, generator)
    infinity_orders, infinity_minors = symbolic_minors(infinity_source, infinity_change)
    assert zero_orders == infinity_orders
    controls = minor_controls(zero_source, zero_change, zero_minors)
    controls += minor_controls(infinity_source, infinity_change, infinity_minors)
    assert controls == 800
    cases = []
    for degree in range(1, 12):
        infinity_order = 12-degree
        zero_caps = [[k, value-(16-2*degree if degree <= 7 and 9 <= k <= 16 else 0)]
                     for k, value in enumerate(prefixes)]
        zero_vertices, zero = lower_hull(zero_caps)
        infinity_caps = [[0, 0]]
        for record in infinity_minors:
            valuations = [value+j*infinity_order for j, value in enumerate(record['valuations'])
                          if value is not None]
            minimum = min(valuations)
            if valuations.count(minimum) == 1:
                infinity_caps.append([record['size'], minimum])
        prefix = 0
        for k, order in enumerate(infinity_orders, 1):
            prefix += order
            if order < infinity_order:
                infinity_caps.append([k, prefix])
        infinity_vertices, infinity = lower_hull(infinity_caps)
        cases.append(dict(polynomial_degree=degree, infinity_gauge_order=infinity_order,
                          zero_vertices=zero_vertices, infinity_vertices=infinity_vertices,
                          pressure=pressure_certificate(zero, infinity)))
    assert 2*8**2 > 11**2
    names = ('check_joint_two_place_obstruction.py', 'check_joint_zero_prefix.py',
             'check_joint_endpoint_loops.py', 'check_joint_graded_newton.py',
             'joint_energy_conjugations.py', 'weighted_smith_oracle_f4.py',
             'search_f4_recentered_crossings.py')
    return dict(stride=STRIDE, zero_orders=zero_orders, infinity_orders=infinity_orders,
                zero_minor_witnesses=zero_minors, infinity_minor_witnesses=infinity_minors,
                independent_F16_determinants=controls, cases=cases, rho_floor='8/11',
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names},
                scope='Every characteristic-two field, every nonzero polynomial q of degree at most eleven with q(1)=0, at the inherited physical heights.',
                conclusion='The complete coupled two-place pencil has typed Smith parameter at least 8/11, strictly above 1/sqrt(2). The zero profile is weakly prefix-dominated by the reference for every q.',
                limit='Other heights, actual coefficient decompositions and other generator directions remain open; q=0 is the original reference and is excluded from the 8/11 floor.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-two-place-obstruction.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ('source_sha256', 'zero_minor_witnesses', 'infinity_minor_witnesses', 'cases')}))


if __name__ == '__main__':
    main()
