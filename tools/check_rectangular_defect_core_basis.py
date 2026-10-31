#!/usr/bin/env python3
"""Certify regular coordinate bases for the rectangular defect's core CAR.

All rational coefficients are checked for integrality. The vacuum-projector
row and column choices are certified modulo z, so the resulting complete
Fock bases have unit determinant over the valuation ring.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_rectangular_defect_support import (
    Rat, Z, car_checks, concentrate_infinity, concentrate_zero,
)


def data_for(place, M, full_core=False):
    if full_core:
        from check_rectangular_full_core_support import full_core_data
        return full_core_data(place, M)
    procedure = concentrate_infinity if place == 'infinity' else concentrate_zero
    images, active, orders, trace = procedure(M)
    return {'place': place, 'M': M, 'active_indices': active,
            'active_orders': [orders[i] for i in active],
            'spectator_orders': [orders[i] for i in range(2*M) if i not in active],
            'active_images': {name: {'creation': [a[i].pair() for i in active],
                                    'annihilation': [b[i].pair() for i in active]}
                              for name, (a, b) in images.items()},
            'regular_shears': trace}


def core_operators(data, side):
    """Return integral active CAR pairs; right pairs are pulled through Sigma."""
    M, infinity = data['M'], data['place'] == 'infinity'
    if side == 'left':
        names = ([f'out_{letter}{i}' for letter in ('v', 'y') for i in (0, M-1)]
                 if infinity else [f'left_{letter}{i}' for letter in (('y', 'x') if data.get('full_core') else ('y',)) for i in (0, M-1)])
    elif side == 'right':
        names = ([f'in_{letter}{i}' for letter in (('u', 'x') if data.get('full_core') else ('u',)) for i in (0, M-1)] if infinity
                 else [f'right_{letter}{i}' for letter in ('v', 'u') for i in (0, M-1)])
    else:
        raise ValueError('Side must be left or right.')
    images = {}
    for name in names:
        for suffix in ('create', 'annihilate'):
            image = data['active_images'][name+'_'+suffix]
            a = [Rat(*entry) for entry in image['creation']]
            b = [Rat(*entry) for entry in image['annihilation']]
            if side == 'right':
                a = [entry*Z**(-order) for entry, order in zip(a, data['active_orders'])]
                b = [entry*Z**order for entry, order in zip(b, data['active_orders'])]
            assert all(entry.order() >= 0 for entry in a+b if entry)
            images[name+'_'+suffix] = (a, b)
    car_checks(images, names)
    return names, images


def residue_operator(pair):
    terms = []
    for creation, vector in ((True, pair[0]), (False, pair[1])):
        for mode, coefficient in enumerate(vector):
            if coefficient and coefficient.order() == 0:
                assert coefficient.n & 1 and coefficient.d & 1
                terms.append((mode, creation))
    return terms


def transpose(matrix):
    result = [0]*len(matrix)
    for i, row in enumerate(matrix):
        while row:
            bit = row & -row
            row ^= bit
            result[bit.bit_length()-1] |= 1 << i
    return result


def left_car(matrix, terms):
    result = [0]*len(matrix)
    for mode, creation in terms:
        bit = 1 << mode
        for source, row in enumerate(matrix):
            if bool(source & bit) != creation:
                result[source ^ bit] ^= row
    return result


def right_car(matrix, terms):
    return transpose(left_car(transpose(matrix), [(mode, not creation) for mode, creation in terms]))


def independent_rows(matrix):
    basis, indices = {}, []
    for i, row in enumerate(matrix):
        while row:
            pivot = (row & -row).bit_length()-1
            if pivot not in basis:
                basis[pivot] = row
                indices.append(i)
                break
            row ^= basis[pivot]
    return indices, sorted(basis)


def certificate(data, side):
    names, images = core_operators(data, side)
    k, r = len(data['active_indices']), len(names)
    n, d = 1 << k, 1 << (k-r)
    operators = {name: residue_operator(pair) for name, pair in images.items()}
    projector = [1 << i for i in range(n)]
    for name in names:
        occupied = left_car(left_car(projector, operators[name+'_annihilate']), operators[name+'_create'])
        projector = [a ^ b for a, b in zip(projector, occupied)]
    I, J = independent_rows(projector)
    assert len(I) == len(J) == d
    minor = [sum(((projector[i] >> j)&1) << t for t, j in enumerate(J)) for i in I]
    assert len(independent_rows(minor)[0]) == d
    L, B = [], [0]*n
    for mask in range(1 << r):
        left, right = projector, projector
        for i, name in enumerate(names):
            if mask >> i & 1:
                left = right_car(left, operators[name+'_annihilate'])
                right = left_car(right, operators[name+'_create'])
        L.extend(left[i] for i in I)
        for row, value in enumerate(right):
            B[row] |= sum(((value >> j)&1) << (mask*d+t) for t, j in enumerate(J))
    assert len(independent_rows(L)[0]) == n
    assert len(independent_rows(B)[0]) == n
    for i, name in enumerate(names):
        for suffix, creation in (('create', True), ('annihilate', False)):
            actual = operators[name+'_'+suffix]
            coordinate = [(k-r+i, creation)]
            assert right_car(L, actual) == left_car(L, coordinate)
            assert left_car(B, actual) == right_car(B, coordinate)
    minimum = min(entry.order() for pair in images.values() for vector in pair for entry in vector if entry)
    digest = sha256(json.dumps({'L': [hex(row) for row in L], 'B': [hex(row) for row in B]},
                               sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return {'place': data['place'], 'M': data['M'], 'side': side,
            'full_core': bool(data.get('full_core')),
            'active_modes': k, 'core_modes': r, 'core_names': names,
            'matrix_order': n, 'multiplicity': d,
            'minimum_coefficient_order': minimum,
            'vacuum_rows': I, 'vacuum_columns': J,
            'vacuum_minor_rank_mod_z': d,
            'left_map_rank_mod_z': n, 'right_map_rank_mod_z': n,
            'intertwining_identities_checked': 4*r,
            'coordinate_core_mode_offset': k-r,
            'mod_z_basis_sha256': digest}


def core_basis_data(data, side):
    """Return the exact regular CAR images and certified vacuum pivot choices.

    Coordinate index ``S*multiplicity+j`` means core occupation mask S and
    vacuum multiplicity index j. Left maps use vacuum_rows; right maps use
    vacuum_columns. No inverse of the vacuum minor is needed.
    """
    result = certificate(data, side)
    _, images = core_operators(data, side)
    result['regular_images'] = {
        name: {'creation': [entry.pair() for entry in pair[0]],
               'annihilation': [entry.pair() for entry in pair[1]]}
        for name, pair in images.items()}
    result['core_letters'] = [name.rsplit('_', 1)[-1][0] for name in result['core_names'][::2]]
    result['core_pair_modes'] = {letter: [2*i, 2*i+1] for i, letter in enumerate(result['core_letters'])}
    result['physical_state_to_occupation_mask'] = [1, 0, 2, 3]
    mode_charges = ([1]*len(data['active_indices']) if data['place'] == 'infinity'
                    else [1 if index < data['M'] else -1 for index in data['active_indices']])
    def charge(mask):
        return sum(value for i, value in enumerate(mode_charges) if mask >> i & 1)
    core_charges = []
    for name in result['core_names']:
        values = {mode_charges[i] for i, entry in enumerate(images[name+'_create'][0]) if entry}
        values |= {-mode_charges[i] for i, entry in enumerate(images[name+'_create'][1]) if entry}
        assert len(values) == 1
        value = values.pop()
        assert all(-mode_charges[i] == value for i, entry in enumerate(images[name+'_annihilate'][0]) if entry)
        assert all(mode_charges[i] == value for i, entry in enumerate(images[name+'_annihilate'][1]) if entry)
        core_charges.append(value)
    result['core_charges'] = core_charges
    result['coordinate_row_charges'] = [
        charge(index)+sum(value for i, value in enumerate(core_charges) if mask >> i & 1)
        for mask in range(1 << len(core_charges)) for index in result['vacuum_rows']]
    result['coordinate_column_charges'] = [
        charge(index)+sum(value for i, value in enumerate(core_charges) if mask >> i & 1)
        for mask in range(1 << len(core_charges)) for index in result['vacuum_columns']]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, help='One support artifact; default checks the power-of-two sequence.')
    parser.add_argument('--full-core', action='store_true', help='Check the extended ten-mode support.')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--artifact-directory', type=Path)
    args = parser.parse_args()
    data = ([json.loads(args.data.read_text())] if args.data else
            [data_for(place, M, args.full_core) for M in (2, 4, 8, 16, 32, 64, 128) for place in ('infinity', 'zero')])
    records = []
    for item in data:
        for side in ('left', 'right'):
            full = core_basis_data(item, side)
            if args.artifact_directory:
                args.artifact_directory.mkdir(parents=True, exist_ok=True)
                target = args.artifact_directory/f"core-basis-{item['place']}-M{item['M']}-{side}.json"
                target.write_text(json.dumps(full, separators=(',', ':'))+'\n')
            records.append({key: value for key, value in full.items()
                            if key not in ('regular_images', 'coordinate_row_charges', 'coordinate_column_charges')})
    result = {'status': 'exact rational integrality and unit-determinant core-coordinate witnesses',
              'cases': records}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
