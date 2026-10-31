#!/usr/bin/env python3
"""Exact ten-mode reductions allowing all four physical extreme-core gauges."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_rectangular_defect_support import (
    ZERO, ONE, Z, UNIT, T, Q, vec, pairing, car_checks, apply_shear,
    infinity_images, zero_images, infinity_matrices,
)


def full_infinity_images(M):
    images = infinity_images(M)
    for i in (0, M-1):
        images[f'in_x{i}_create'] = (
            vec(2*M, [(a, UNIT*Z**(1-2*M+2*(a | i))) for a in range(M)]+
                [(M+a, Z**(1+2*a)/UNIT) for a in range(M) if a & i == i]), vec(2*M))
        images[f'in_x{i}_annihilate'] = (
            vec(2*M), vec(2*M, [(M+a, UNIT*Z**(-1-2*a)) for a in range(M) if a & i == a]))
    return images


def full_zero_images(M):
    images = zero_images(M)
    for i in (0, M-1):
        images[f'left_x{i}_create'] = (vec(2*M, [(i, ONE)]), vec(2*M))
        images[f'left_x{i}_annihilate'] = (
            vec(2*M, [(M+a, UNIT*T**(i & a)) for a in range(M)]), vec(2*M, [(i, ONE)]))
    return images


def concentrate_two_sequences(images, offset, M, variable, orders, trace, annihilation):
    """Retain four coordinates for creations {e0,eT,1,s^a} and one covector."""
    high = M-2
    outside = list(range(2, high))
    for j in outside+[high]:
        apply_shear(images, offset+j, offset+1, variable**(j-1), orders, trace)
    for j in outside:
        coefficient = (ONE+variable**(j-1))/(ONE+variable**(high-1))
        apply_shear(images, offset+j, offset+high, coefficient, orders, trace)
    if annihilation == 'constant':
        pivot = ONE
        for j in outside+[high]:
            pivot = pivot+variable**(j-1)
        assert pivot.order() == 0
        for j in outside:
            apply_shear(images, offset+1, offset+j, ONE/pivot, orders, trace)
    elif annihilation == 'inverse':
        pivot = variable**(-high)
        for j in outside:
            coefficient = (ONE+variable**(j-1))/(ONE+variable**(high-1))
            pivot = pivot+coefficient*variable**(-j)
        assert pivot.order() == -2*high
        for j in outside:
            apply_shear(images, offset+high, offset+j, variable**(-j)/pivot, orders, trace)
    else:
        raise ValueError('The one spreading annihilation covector must be constant or inverse.')
    return {offset+j for j in (0, 1, high, M-1)}


def concentrate_four_sequences(images, offset, M, variable, orders, trace):
    """Retain six coordinates for {e0,eT,1,s^a} and {e0,eT,1,s^-a}."""
    high, high_b, low_b = M-2, M-3, 2
    kept = {0, 1, low_b, high_b, high, M-1}
    outside = [j for j in range(M) if j not in kept]
    for j in outside+[high]:
        apply_shear(images, offset+j, offset+1, variable**(j-1), orders, trace)
    for j in outside:
        apply_shear(images, offset+j, offset+high,
                    (ONE+variable**(j-1))/(ONE+variable**(high-1)), orders, trace)
    for j in outside+[low_b]:
        apply_shear(images, offset+high_b, offset+j, variable**(high_b-j), orders, trace)
    for j in outside:
        apply_shear(images, offset+low_b, offset+j,
                    (ONE+variable**(high_b-j))/(ONE+variable**(high_b-low_b)), orders, trace)
    return {offset+j for j in kept}


def concentrate_full_infinity(M):
    images = full_infinity_images(M)
    orders = [1-2*M+2*a for a in range(2*M)]
    trace, diagonal = [], []
    if M < 8:
        return images, list(range(2*M)), orders, trace, diagonal
    # This regular diagonal unit keeps Sigma fixed with its inverse companion.
    for a in range(M):
        factor = Q**(-a)
        for creation, annihilation in images.values():
            creation[M+a] = creation[M+a]*factor
            annihilation[M+a] = annihilation[M+a]/factor
        diagonal.append({'index': M+a, 'coefficient': factor.pair()})
    first = concentrate_two_sequences(images, 0, M, T, orders, trace, 'constant')
    second = concentrate_four_sequences(images, M, M, T/Q, orders, trace)
    active = sorted(first | second)
    assert all(not a[i] and not b[i] for a, b in images.values() for i in range(2*M) if i not in active)
    return images, active, orders, trace, diagonal


def concentrate_full_zero(M):
    images = full_zero_images(M)
    orders = [2*a+1 for _ in range(2) for a in range(M)]
    trace = []
    if M < 8:
        return images, list(range(2*M)), orders, trace, []
    first = concentrate_four_sequences(images, 0, M, T, orders, trace)
    second = concentrate_two_sequences(images, M, M, T, orders, trace, 'inverse')
    active = sorted(first | second)
    assert all(not a[i] and not b[i] for a, b in images.values() for i in range(2*M) if i not in active)
    return images, active, orders, trace, []


def full_core_data(place, M):
    procedure = concentrate_full_infinity if place == 'infinity' else concentrate_full_zero
    images, active, orders, trace, diagonal = procedure(M)
    return {'place': place, 'M': M, 'full_core': True,
            'active_indices': active, 'active_orders': [orders[i] for i in active],
            'spectator_orders': [orders[i] for i in range(2*M) if i not in active],
            'active_images': {name: {'creation': [a[i].pair() for i in active],
                                    'annihilation': [b[i].pair() for i in active]}
                              for name, (a, b) in images.items()},
            'regular_shears': trace, 'regular_diagonal': diagonal}


def check_x_transport(M):
    _, _, _, V, Vi, sigma = infinity_matrices(M)
    images = full_infinity_images(M)
    for i in (0, M-1):
        create, _ = images[f'in_x{i}_create']
        _, annihilate = images[f'in_x{i}_annihilate']
        assert create == [sigma[j]*Vi[j][M+i] for j in range(2*M)]
        assert annihilate == [V[M+i][j]/sigma[j] for j in range(2*M)]
    return {'M': M, 'independent_transport_components_checked': 8*M}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--artifact-directory', type=Path)
    args = parser.parse_args()
    records = []
    for M in (2, 4, 8, 16, 32, 64, 128):
        for place, initial, procedure, groups in (
            ('infinity', full_infinity_images, concentrate_full_infinity,
             [[f'out_{letter}{i}' for letter in ('v', 'y') for i in (0, M-1)],
              [f'in_{letter}{i}' for letter in ('u', 'x') for i in (0, M-1)]]),
            ('zero', full_zero_images, concentrate_full_zero,
             [[f'left_{letter}{i}' for letter in ('y', 'x') for i in (0, M-1)],
              [f'right_{letter}{i}' for letter in ('v', 'u') for i in (0, M-1)]])):
            before = initial(M)
            images, active, orders, trace, diagonal = procedure(M)
            checks = sum(car_checks(before, group)+car_checks(images, group) for group in groups)
            for name, first in before.items():
                for other, second in before.items():
                    assert pairing(first, second) == pairing(images[name], images[other])
            data = full_core_data(place, M)
            if args.artifact_directory:
                args.artifact_directory.mkdir(parents=True, exist_ok=True)
                target = args.artifact_directory/f'full-core-{place}-M{M}.json'
                target.write_text(json.dumps(data, separators=(',', ':'))+'\n')
            records.append({'place': place, 'M': M, 'active_modes': len(active),
                            'spectator_modes': 2*M-len(active), 'regular_shears': len(trace),
                            'regular_diagonal_units': len(diagonal), 'canonical_pairings_checked': checks,
                            'all_cross_pairings_checked': len(images)**2,
                            'coefficient_sha256': sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()})
    result = {'status': 'exact full-core support in at most ten active modes; no pressure claim',
              'original_x_transport_checks': [check_x_transport(M) for M in (2, 4, 8, 16)],
              'cases': records}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
