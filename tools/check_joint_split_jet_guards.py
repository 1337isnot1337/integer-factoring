#!/usr/bin/env python3
"""Exact spectral guards and formal two-jet rigidity for joint 25-to-16 splits."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import comb, gcd
from pathlib import Path

from check_binary_coefficient_entropy import binary_rank
from check_joint_energy_conjugations import polynomial
from check_joint_grid_compressions import multiply, power
from joint_energy_conjugations import kernel, model
from joint_five_letter_splits import HEIGHTS, coordinate_start, crossing, profiles
from physical_split_crossings import linear_combination
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]


def positive_certificate(first, second, p, q):
    f = {0: 1, 1: 1, 3: 1, 4: 1}
    left = multiply({p*k: v for k, v in first.items()}, {p*k: v for k, v in second.items()})
    factor = gcd(q, 4*p)
    difference = Counter(power(left, q//factor))
    for key, value in power({q*k: v for k, v in f.items()}, 4*p//factor).items():
        difference[key] -= value
    difference = {key: value for key, value in difference.items() if value}
    order, degree = min(difference), max(difference)-min(difference)
    terms = {key-order: value for key, value in difference.items()}
    beta = [sum(value*comb(degree-key, i-key) for key, value in terms.items() if key <= i)
            for i in range(degree+1)]
    polynomial_degree = degree
    while min(beta) < 0 and degree < polynomial_degree+4096:
        beta = [beta[0]]+[beta[i]+beta[i-1] for i in range(1, len(beta))]+[beta[-1]]
        degree += 1
    assert min(beta) >= 0 and beta[0] > 0 and beta[-1] > 0
    return dict(rho_floor=[p, q], polynomial=difference, monomial_order=order,
                polynomial_degree=polynomial_degree, Bernstein_degree=degree,
                minimum_integer_coefficient=min(beta),
                coefficient_sha256=sha256(json.dumps(beta, separators=(',', ':')).encode()).hexdigest())


def coefficient(data, side, jet):
    result = []
    for r, row in enumerate(data['operator']):
        values = []
        for c, (value, high) in enumerate(row):
            assert high == 0
            degree = 8+jet-data['energies'][r] if side == 0 else data['energies'][c]+8-jet
            values.append(int(degree >= 0 and (value >> degree) & 1))
        result.append(values)
    return result


def tangent_certificate(data, selected, baseline):
    heights = [a+b for a in HEIGHTS for b in HEIGHTS]
    variables = [(kind, j, r) for kind in (0, 1) for j in range(16) for r in range(25)
                 if data['energies'][r] == heights[j]]
    identity = [[int(r == c) for c in range(25)] for r in range(25)]
    records = []
    for depth in (1, 2):
        matrices = [identity]+[coefficient(data, side, jet) for side in (0, 1) for jet in range(depth)]
        equations = [0]*(len(matrices)*256)
        for v, (kind, j, r) in enumerate(variables):
            for m, matrix in enumerate(matrices):
                for i in range(16):
                    row, col, value = (i, j, matrix[selected[i]][r]) if kind == 0 else (j, i, matrix[r][selected[i]])
                    if value:
                        equations[m*256+row*16+col] |= 1 << v
        nullspace = kernel(equations, len(variables))
        assert binary_rank(equations)+len(nullspace) == len(variables)
        outside = [v for v, (_, _, r) in enumerate(variables) if r not in selected]
        images = [sum(((vector >> v) & 1) << i for i, v in enumerate(outside)) for vector in nullspace]
        record = dict(frozen_jets_per_cut=depth, variables=len(variables),
                      equation_rank=binary_rank(equations), kernel_basis=nullspace,
                      extrinsic_image_dimension=binary_rank(images))
        if depth == 2:
            assert len(variables) == 88 and len(nullspace) == 7 and binary_rank(images) == 0
            generators = []
            for vector in nullspace:
                x = {(j, r): (vector >> v) & 1 for v, (kind, j, r) in enumerate(variables) if kind == 0}
                c = [[x.get((j, selected[i]), 0) for j in range(16)] for i in range(16)]
                for v, (kind, j, r) in enumerate(variables):
                    expected = 0 if r not in selected else c[selected.index(r)][j] if kind == 0 else c[j][selected.index(r)]
                    assert ((vector >> v) & 1) == expected
                assert all(not c[i][j] or heights[i] == heights[j] for i in range(16) for j in range(16))
                for i, j in product(range(16), repeat=2):
                    left = linear_combination((baseline[k][j], c[i][k]) for k in range(16))
                    right = linear_combination((baseline[i][k], c[k][j]) for k in range(16))
                    assert left == right
                generators.append(sum(c[i][j] << (16*i+j) for i in range(16) for j in range(16)))
            assert binary_rank(generators) == 7
            record['whole_polynomial_commutant_generators'] = generators
        records.append(record)
    assert records[0]['extrinsic_image_dimension'] == 2
    return records


def fiber_residual(data, selected):
    matrices = []
    for jet in (0, 1):
        source = coefficient(data, 1, jet)
        result = [[0]*16 for _ in range(16)]
        for x, y, u, v in product(range(4), repeat=4):
            i, j = 4*x+y, 4*v+u
            value = source[selected[i]][selected[j]]
            if i == 10:
                value ^= (source[8][selected[j]] ^ source[16][selected[j]]) << 1
            if j == 10:
                value ^= (source[selected[i]][8] ^ source[selected[i]][16]) << 5
            if i == j == 10:
                value ^= (source[8][8] ^ source[8][16] ^ source[16][8] ^ source[16][16]) << 6
            result[4*x+u][4*v+y] = value
        matrices.append(result)
    leading, first = matrices
    assert all(value in (0, 1) for row in leading for value in row)
    assert leading[0][0] == 1 and binary_rank([sum(value << j for j, value in enumerate(row)) for row in leading]) == 1
    residual = [[first[i][j] ^ (first[0][j] if leading[i][0] else 0)
                 ^ (first[i][0] if leading[0][j] else 0)
                 ^ (first[0][0] if leading[i][0] and leading[0][j] else 0)
                 for j in range(1, 16)] for i in range(1, 16)]
    expected = {(i, j): 1 for i in (1, 4, 7, 13) for j in (1, 4, 7, 13)}
    expected.update({(8, 14): 2, (11, 2): 2, (2, 11): 32, (14, 8): 32})
    assert {(i+1, j+1): value for i, row in enumerate(residual) for j, value in enumerate(row) if value} == expected
    return dict(polynomial_encoding='bit a+5*b encodes alpha^a beta^b; no Boolean reduction',
                nonzero_entries=[[i, j, value] for (i, j), value in sorted(expected.items())],
                rank='1+2*[alpha!=0]+2*[beta!=0] over every characteristic-two field')


def certificate():
    data = model()
    p, l = coordinate_start()
    selected = [column.index(1) for column in p]
    baseline = crossing(p, l)
    assert 2*sum(a+b for a in HEIGHTS for b in HEIGHTS) == 128
    for i, row in enumerate(data['operator']):
        for j, (low, high) in enumerate(row):
            assert high == 0 and not ((low << data['energies'][i]) & 255)
            assert low.bit_length()-1 <= data['energies'][j]+8
    core = [[baseline[4*x+y][4*v+u] for v, u in product(range(2), repeat=2)]
            for x, y in product(range(2), repeat=2)]
    opposite_core = [[baseline[4*x+y][4*v+u] for v, y in product(range(2), repeat=2)]
                     for x, u in product(range(2), repeat=2)]
    core_zero, _ = WeightedSmithOracle([[polynomial(value) for value in row] for row in core]).profile(
        [0, 1, 1, 2], [0]*4)
    core_infinity, _ = WeightedSmithOracle([[polynomial(value) for value in row] for row in opposite_core], infinity=True).profile(
        [0, 1, 0, 1], [0, 0, 1, 1])
    assert [value-8 for value in core_zero] == [value+8 for value in core_infinity] == [0, 1, 3, 4]
    # Both zero jets and the leading infinity jet are independent of the
    # two formal fiber parameters. The next infinity jet is computed below.
    for side, jet in ((0, 0), (0, 1), (1, 0)):
        matrix = coefficient(data, side, jet)
        assert all(not (matrix[8][s] ^ matrix[16][s]) for s in selected)
        assert all(not (matrix[s][8] ^ matrix[s][16]) for s in selected)
        assert not (matrix[8][8] ^ matrix[8][16] ^ matrix[16][8] ^ matrix[16][16])
    tangent = tangent_certificate(data, selected, baseline)
    symbolic = fiber_residual(data, selected)
    q = {0: 1, 1: 1, 3: 1, 4: 1, 10: 12}
    floors = dict(two_zero=positive_certificate(q, {0: 2, 4: 2, 10: 12}, 5, 7),
                  two_one=positive_certificate(q, {0: 1, 1: 2, 6: 1, 10: 12}, 5, 7),
                  nonzero_fiber=positive_certificate(q, {0: 1, 1: 3, 10: 7, 11: 5}, 3, 4))
    assert 2*5*5 > 7*7
    controls = []
    for a, b in product(range(4), repeat=2):
        pp, ll = [row[:] for row in p], [row[:] for row in l]
        for j in (8, 16):
            pp[10][j] ^= b
            ll[10][j] ^= a
        operator = crossing(pp, ll)
        native = profiles(operator)
        opposite = [[None]*16 for _ in range(16)]
        for x, y, u, v in product(range(4), repeat=4):
            opposite[4*x+u][4*v+y] = operator[4*x+y][4*v+u]
        orders, _ = WeightedSmithOracle([[polynomial(value) for value in row] for row in operator]).profile(
            [x+y for x in HEIGHTS for y in HEIGHTS], [0]*16)
        assert [value-8 for value in orders] == native['zero_orders_t']
        orders, _ = WeightedSmithOracle([[polynomial(value) for value in row] for row in opposite], infinity=True).profile(
            [y for _ in HEIGHTS for y in HEIGHTS], [x for x in HEIGHTS for _ in HEIGHTS])
        assert [value+8 for value in orders] == native['infinity_orders_t']
        assert native['infinity_orders_t'].count(0) == 1
        assert native['infinity_orders_t'].count(1) == 1+2*bool(a)+2*bool(b)
        controls.append(dict(alpha=a, beta=b, **native))
    names = (Path(__file__).name, 'joint_five_letter_splits.py', 'joint_energy_conjugations.py',
             'search_f4_recentered_crossings.py', 'physical_split_crossings.py',
             'check_joint_grid_compressions.py', 'weighted_smith_oracle_f4.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    return dict(tangent=tangent, symbolic_residual=symbolic, pressure_floors=floors,
                F4_fiber_controls=controls, independent_fraction_free_orders=512,
                protected_core_both_profiles=[0, 1, 3, 4], protected_core_order_controls=8,
                source_entry_bounds_checked=625,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names},
                scope='Protected-core degree-budget spectral guards; formal rigidity only at the stated coordinate split with both first two matrix jets fixed. No global field-valued classification or factoring improvement.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-split-jet-guards.json')
    args = parser.parse_args()
    result = certificate()
    # Integer polynomial exponent keys have a canonical JSON representation.
    result = json.loads(json.dumps(result))
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps(dict(independent_orders=result['independent_fraction_free_orders'],
                          tangent_dimensions=[len(row['kernel_basis']) for row in result['tangent']],
                          extrinsic_dimensions=[row['extrinsic_image_dimension'] for row in result['tangent']],
                          Bernstein_degrees={key: value['Bernstein_degree'] for key, value in result['pressure_floors'].items()})))


if __name__ == '__main__':
    main()
