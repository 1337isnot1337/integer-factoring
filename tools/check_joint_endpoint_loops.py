#!/usr/bin/env python3
"""Endpoint-singular rational crossings and an all-field rho >= 3/4 floor."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import accumulate, product
import json
from math import comb
from pathlib import Path

from check_joint_energy_conjugations import cut_matrix, weighted_coefficient
from check_joint_energy_riccati import assembly_controls, equations as riccati_equations, supported_rows
from check_joint_evaluation_algebra import polynomial_product
from check_joint_leading_tangent import realign
from joint_energy_conjugations import (N, Space, binary_product, commutator_equations,
                                      kernel, model, scale, unit_characters)
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]


def xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def operator_coefficients(data):
    return [sum(((data['operator'][r][c][0] >> k) & 1) << (N*r+c)
                for r in range(N) for c in range(N)) for k in range(19)]


def hasse(values):
    return [xor(value for j, value in enumerate(values) if j & k == k) for k in range(len(values))]


def endpoint_values(jets, changes, parameters):
    modified = []
    for j, value in enumerate(jets):
        low, high = value, 0
        for pole, parameter in enumerate(parameters, 1):
            a, b = scale((changes[j+pole] if j+pole < len(changes) else 0, 0), parameter)
            low, high = low ^ a, high ^ b
        modified.append((low, high))
    # The Hasse transform at one is its own inverse in characteristic two.
    coefficients = [(xor(pair[0] for j, pair in enumerate(modified) if j & k == k),
                     xor(pair[1] for j, pair in enumerate(modified) if j & k == k))
                    for k in range(19)]
    values = []
    for index in range(625):
        x, y, u, v = ((index//5**i) % 5 for i in range(4))
        position = N*(x+5*y)+v+5*u
        values.append(tuple(sum(((pair[plane] >> position) & 1) << k for k, pair in enumerate(coefficients))
                            for plane in range(2)))
    return values


def regular_matrices(data, generator):
    result = [[[0]*N for _ in range(N)] for _ in range(2)]
    for power, coefficient in data['weighted'].items():
        matrices = (coefficient, binary_product(generator, coefficient) ^ binary_product(coefficient, generator))
        for target, matrix in zip(result, matrices):
            matrix = realign(matrix)
            while matrix:
                bit = matrix & -matrix
                r, c = divmod(bit.bit_length()-1, N)
                target[r][c] ^= 1 << (power+8)
                matrix ^= bit
    return result


def symbolic_minor(source, change, rows):
    size = len(rows)
    determinants = [None]*(1 << size)
    determinants[0] = [1]
    for mask in range(1, 1 << size):
        row, result = rows[mask.bit_count()-1], [0]*(size+1)
        for j, column in enumerate(rows):
            if not mask >> j & 1:
                continue
            for k, value in enumerate(determinants[mask ^ (1 << j)]):
                if value:
                    result[k] ^= polynomial_product(value, source[row][column])
                    result[k+1] ^= polynomial_product(value, change[row][column])
        while result and not result[-1]:
            result.pop()
        determinants[mask] = result
    return determinants[-1]


def integer_product(first, second):
    result = [0]*(len(first)+len(second)-1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i+j] += a*b
    return result


def pressure_certificate(zero, half):
    p = [0]*(3*(max(zero)-8)+1)
    for order in zero:
        p[3*(order-8)] += 1
    lower = [0]*10
    lower[0] = lower[3] = 1
    lower[9] = 7
    first = integer_product(p, lower)
    energy = [0]*(4*max(half)+1)
    for value in half:
        energy[4*value] += 1
    second = integer_product(integer_product(energy, energy), energy)
    for i, value in enumerate(second):
        first[i] -= value
    assert len(first)-1 == 69 and first[:3] == [0]*3
    quotient = first[3:]
    assert len(quotient)-1 == 66
    bernstein = [sum(quotient[i]*comb(66-i, k-i) for i in range(k+1)) for k in range(67)]
    assert min(bernstein) == 2
    return dict(difference_polynomial=first, zero_factor=3, quotient_degree=66,
                scaled_Bernstein_coefficients=bernstein, all_coefficients_positive=True,
                rho_floor='3/4')


def all_rank_reduction(data, jets, radical):
    # Shift D_1,...,D_18 to the generic quadratic builder's C_0,...,C_17.
    weighted = {-8+j: value for j, value in enumerate(jets[1:])}
    initial = [1 << (N*r+c) for r, c in data['positions']]
    records = []
    for pole_order in range(1, 18):
        basis, stages = initial[:], []
        while True:
            size = len(basis)
            blocks, monomials, dimension = riccati_equations(basis, pole_order, weighted)
            rows = [row for block in blocks for row in block if row]
            linear = supported_rows(rows, list(range(size)), dimension)
            square = supported_rows(rows, [monomials[i, i] for i in range(size)], dimension)
            remaining = kernel(linear+square, size)
            stages.append(dict(input_dimension=size, linear_consequences=linear,
                               square_root_consequences=square, remaining_dimension=len(remaining)))
            if len(remaining) == size:
                assert size == 3 and not rows
                break
            space = Space(basis)
            basis = [space.matrix((v, 0))[0] for v in remaining]
        assert len(Space(basis+radical).basis) == 3
        records.append(dict(pole_order=pole_order, stages=stages, terminal_basis=basis))
    return records, assembly_controls(initial, weighted)


def certificate():
    data = model()
    coefficients = operator_coefficients(data)
    jets = hasse(coefficients)
    endpoint_vectors = kernel(commutator_equations(jets[1], data['positions']), 69)
    endpoint = [sum(1 << (N*r+c) for i, (r, c) in enumerate(data['positions']) if v >> i & 1)
                for v in endpoint_vectors]
    assert len(endpoint) == 5
    assert all(binary_product(a, b) == binary_product(b, a) for a, b in product(endpoint, repeat=2))
    squares = [binary_product(a, a) for a in endpoint]
    equations = [sum(((a >> i) & 1) << j for j, a in enumerate(squares)) for i in range(N*N)]
    radical = [Space(endpoint).matrix((v, 0))[0] for v in kernel(equations, 5)]
    assert len(radical) == 3
    assert all(binary_product(a, b) == 0 for a, b in product(radical, repeat=2))
    _, full_characters = unit_characters(data)
    trivial = [data['full'].matrix((v, 0))[0] for v in kernel(full_characters, 4)]
    generator = sum(1 << (N*r+c) for pair in ((4, 20), (7, 11), (8, 16)) for r, c in product(pair, repeat=2))
    assert len(Space(trivial+[generator]+radical).basis) == 3
    assert len(Space(trivial+[generator]).basis) == 3
    assert len(Space(data['full'].basis+[generator]).basis) == 5
    identities = 0
    for coefficient in coefficients:
        for i, a in enumerate(radical):
            assert binary_product(binary_product(a, coefficient), a) == 0
            identities += 1
            for b in radical[i+1:]:
                assert binary_product(binary_product(a, coefficient), b) == binary_product(binary_product(b, coefficient), a)
                identities += 1
    changes = [binary_product(generator, a) ^ binary_product(a, generator) for a in jets]
    assert changes[:4] == [0]*4 and changes[4]
    all_rank_cases, all_rank_controls = all_rank_reduction(data, jets, radical)
    source, change = regular_matrices(data, generator)
    minors = []
    for rows, expected in (([0, 3], [1]), ([0, 2, 3, 7, 8, 11, 12, 13, 17], [52, None, 34, None, 18, None, 4])):
        determinant = symbolic_minor(source, change, rows)
        orders = [(a & -a).bit_length()-1 if a else None for a in determinant]
        assert orders == expected
        if len(rows) == 2:
            assert determinant == [6]
        else:
            for order in (1, 2, 3):
                candidates = [(v+order*i, i) for i, v in enumerate(orders) if v is not None]
                assert min(candidates) == (4+6*order, 6)
                assert sum(v == 4+6*order for v, _ in candidates) == 1
        minors.append(dict(rows=rows, columns=rows, determinant_in_g=determinant,
                           coefficient_valuations=orders))
    half, energies = data['half'], data['energies']
    weights = [half[i//5] for i in range(N)]
    baseline, _ = WeightedSmithOracle(cut_matrix(data['values'], True), infinity=True).profile(weights, weights)
    zero, _ = WeightedSmithOracle(cut_matrix(data['values'], False)).profile(energies, [0]*N)
    assert zero[0] == 8
    pressure = pressure_certificate(zero, half)
    profiles, records = [], []
    for parameters in product(range(4), repeat=3):
        values = endpoint_values(jets, changes, parameters)
        for index, (a, b) in enumerate(values):
            x, y, u, v = ((index//5**i) % 5 for i in range(4))
            assert ((a.bit_count() % 2) | ((b.bit_count() % 2) << 1)) == int(x == v and y == u)
        assert weighted_coefficient(values, energies, -8) == (data['weighted'][-8], 0)
        orders, _ = WeightedSmithOracle(cut_matrix(values, True), infinity=True).profile(weights, weights)
        actual_zero, _ = WeightedSmithOracle(cut_matrix(values, False)).profile(energies, [0]*N)
        assert actual_zero == zero
        assert all(a <= b for a, b in zip(accumulate(orders), accumulate(baseline)))
        if orders not in profiles:
            profiles.append(orders)
        records.append([list(parameters), profiles.index(orders)])
    assert len(profiles) == 4 and len(records) == 64
    return dict(endpoint_commutant_dimension=5, nilpotent_sector_dimension=3,
                trivial_sector_dimension=2, effective_parameters=3, generator=generator,
                quadratic_matrix_identities=identities, commutator_zero_order_at_one=4,
                all_rank_pole_cases=all_rank_cases, all_rank_dense_F4_entry_controls=all_rank_controls,
                symbolic_minors=minors, pressure_certificate=pressure,
                reference_zero=list(zero), F4_profiles=[list(p) for p in profiles], F4_records=records,
                fraction_free_profiles=130, all_F4_profiles_prefix_dominated=True,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('check_joint_endpoint_loops.py', 'joint_energy_conjugations.py',
                                            'check_joint_evaluation_algebra.py', 'check_joint_energy_conjugations.py',
                                            'check_joint_energy_riccati.py',
                                            'check_joint_leading_tangent.py', 'weighted_smith_oracle_f4.py',
                                            'search_f4_recentered_crossings.py')},
                scope='All characteristic-two fields and every rank: U=I+h(t)P with any constant energy-block P satisfying P^2=0, rational h, regular at zero and identity at infinity, yielding a strict polynomial crossing. Also arbitrary rational gauges valued in the endpoint nilpotent sector.',
                conclusion='Every nontrivial gauge reduces to the three endpoint principal parts. Other finite poles force a trivial generator by single-evaluation algebra; all endpoint pole orders reduce by quadratic/Frobenius elimination or the degree bound. Every nontrivial tensor has typed Smith parameter at least 3/4 at the natural heights.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-endpoint-loops.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items() if k not in (
        'source_sha256', 'symbolic_minors', 'pressure_certificate', 'generator', 'F4_profiles', 'F4_records',
        'reference_zero', 'all_rank_pole_cases')}))


if __name__ == '__main__':
    main()
