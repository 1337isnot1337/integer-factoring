#!/usr/bin/env python3
"""Full physical CAR support for arbitrary hot reflection-core state gauges.

Both physical creation and annihilation operators are concentrated by
regular Gaussian-preserving changes. The uniform bound is min(M,11k).
"""
from __future__ import annotations

from functools import lru_cache

from check_rectangular_defect_support import ZERO, ONE, Z, apply_shear, car_checks, pairing
from reflection_hot_block_profiles import validate_core


def full_images(interval, ramification=1, core=1, internal_regular=False):
    validate_core(interval, core)
    if not isinstance(ramification, int) or ramification < 1:
        raise ValueError('Use a positive integer substitution degree.')
    m, p = interval // 2, ramification
    c, s, u = ONE + Z ** (2 * p), Z ** (4 * p), Z ** (2 * p)
    f, b = Z ** (2 * interval * p), ONE + Z ** (2 * interval * p)
    answer = {}
    def vector():
        return [ZERO] * interval
    for j in range(core):
        # All indices below are reversed diagonal indices r=m-1-i.
        cre, ann = vector(), vector();cre[m + j] = ONE
        answer[f'left_first{j}_create'] = (cre, ann)
        cre, ann = vector(), vector()
        for r in range(m):
            cre[r] = c * s ** (j & r)
            ann[r] = c * s ** (j - (j & r))
        ann[m + j] = ONE
        answer[f'left_first{j}_annihilate'] = (cre, ann)
        cre, ann = vector(), vector();cre[j] = ONE
        for r in range(m):
            cre[m + r] = c * s ** (r - (j & r))
        answer[f'left_second{j}_create'] = (cre, ann)
        cre, ann = vector(), vector();cre[m - 1 - j] = b;ann[j] = ONE
        for r in range(m):
            cre[m + r] = c * f * s ** (j & r)
        answer[f'left_second{j}_annihilate'] = (cre, ann)
        cre, ann = vector(), vector();cre[j] = u * s ** j
        ann[m - 1 - j] = b / (u * s ** (m - 1 - j))
        for r in range(m):
            ann[m + r] = c * s ** ((j & r) - r) / u
        answer[f'right_first{j}_create'] = (cre, ann)
        cre, ann = vector(), vector();ann[j] = ONE / (u * s ** j)
        for r in range(m):
            ann[m + r] = c * s ** (-(j & r)) / (f * u)
        answer[f'right_first{j}_annihilate'] = (cre, ann)
        cre, ann = vector(), vector();cre[m + j] = f * u * s ** j
        for r in range(m):
            cre[r] = c * u * s ** (r + j - (j & r))
            ann[r] = c * s ** ((j & r) - r) / u
        answer[f'right_second{j}_create'] = (cre, ann)
        cre, ann = vector(), vector();ann[m + j] = ONE / (f * u * s ** j)
        answer[f'right_second{j}_annihilate'] = (cre, ann)
    if internal_regular:
        # Remove the regular internal Gaussian factors before transporting
        # the state algebra. This is usable only for weights under which
        # those factors and their inverses remain regular.
        for name, (cre, ann) in list(answer.items()):
            if name.startswith('left') and name.endswith('annihilate'):
                answer[name] = (vector(), ann)
            if name.startswith('right') and name.endswith('create'):
                answer[name] = (cre, vector())
    return answer


def concentrate_residue(images, offset, residue, step, count, variable, orders, trace, both_ends):
    """Four exact sweeps for two spreading vectors and their duals."""
    coordinate = lambda j: offset + residue + step * j
    high = count - 2 if both_ends else count - 1
    high_dual = count - 3 if both_ends else count - 2
    retained = {0, 1, 2, high_dual, high}
    if both_ends:
        retained.add(count - 1)
    # First concentrate the geometric and constant creation sequences.
    moving = [j for j in range(2, count) if not (both_ends and j == count - 1)]
    for j in moving:
        apply_shear(images, coordinate(j), coordinate(1), variable ** (j - 1), orders, trace)
    for j in moving:
        if j != high:
            apply_shear(images, coordinate(j), coordinate(high),
                        (ONE + variable ** (j - 1)) / (ONE + variable ** (high - 1)), orders, trace)
    # Then concentrate the inverse-geometric and constant annihilators.
    outside = [j for j in range(count) if j not in retained]
    for j in outside + [2]:
        apply_shear(images, coordinate(high_dual), coordinate(j),
                    variable ** (high_dual - j), orders, trace)
    for j in outside:
        apply_shear(images, coordinate(2), coordinate(j),
                    (ONE + variable ** (high_dual - j)) / (ONE + variable ** (high_dual - 2)), orders, trace)
    return {coordinate(j) for j in retained}


def concentrate(interval, ramification=1, core=1, internal_regular=False):
    before = full_images(interval, ramification, core, internal_regular)
    images = {name: (list(first), list(second)) for name, (first, second) in before.items()}
    m, k = interval // 2, core
    orders = [ramification * (2 + 4 * r) for r in range(interval)]
    trace = []
    count = m // k
    if internal_regular:
        variable = Z ** (4 * ramification * k)
        active = set()
        if count < 2:
            return before, images, list(range(interval)), orders, trace
        geometric = sum((variable ** j for j in range(count-1)), ZERO)
        for offset in (0, m):
            for residue in range(k):
                pivot = offset + residue + k
                active |= {offset + residue, pivot}
                for j in range(2, count):
                    apply_shear(images, offset + residue + k*j, pivot,
                                variable ** (j-1), orders, trace)
                for j in range(2, count):
                    apply_shear(images, pivot, offset + residue + k*j,
                                ONE / geometric, orders, trace)
        active = sorted(active)
        assert all(not cre[i] and not ann[i] for cre, ann in images.values()
                   for i in range(interval) if i not in active)
        return before, images, active, orders, trace
    if count < 8:
        return before, images, list(range(interval)), orders, trace
    variable = Z ** (4 * ramification * k)
    active = set()
    for residue in range(k):
        active |= concentrate_residue(images, 0, residue, k, count, variable, orders, trace, True)
        active |= concentrate_residue(images, m, residue, k, count, variable, orders, trace, False)
    active = sorted(active)
    assert len(active) == 11 * core
    assert all(not cre[i] and not ann[i] for cre, ann in images.values()
               for i in range(interval) if i not in active)
    return before, images, active, orders, trace


@lru_cache(maxsize=32)
def core_data(interval, place='zero', ramification=1, core=1, internal_regular=False):
    if place not in ('zero', 'infinity'):
        raise ValueError('Choose zero or infinity.')
    _, images, active, orders, trace = concentrate(interval, ramification, core, internal_regular)
    left = ('x', 'y') if place == 'zero' else ('u', 'x')
    right = ('u', 'v') if place == 'zero' else ('v', 'y')
    names = {}
    renamed = {}
    for side, letters in (('left', left), ('right', right)):
        names[side] = []
        for which, letter in zip(('first', 'second'), letters):
            for j in range(core):
                prefix = f'{side}_{letter}{j}'
                names[side].append(prefix)
                for kind in ('create', 'annihilate'):
                    creation, annihilation = images[f'{side}_{which}{j}_{kind}']
                    renamed[prefix + '_' + kind] = {
                        'creation': [creation[i].pair() for i in active],
                        'annihilation': [annihilation[i].pair() for i in active]}
    return {'M': interval, 'place': place, 'ramification': ramification, 'core': core,
            'internal_gaussians_removed': internal_regular,
            'active_indices': active, 'active_orders': [orders[i] for i in active],
            'spectator_orders': [order for i, order in enumerate(orders) if i not in active],
            'active_images': renamed, 'core_names': names, 'regular_shears': trace,
            'logical_core_indices': [interval // 2 - 1 - j for j in range(core)],
            'scalar_order': (1 if place == 'zero' else -1) * ramification * interval * interval // 2}


def check_support(interval, ramification=1, core=1, internal_regular=False):
    before, after, active, orders, trace = concentrate(interval, ramification, core, internal_regular)
    groups = [[f'{side}_{which}{j}' for which in ('first', 'second') for j in range(core)]
              for side in ('left', 'right')]
    canonical = sum(car_checks(before, group) + car_checks(after, group) for group in groups)
    for first, a in before.items():
        for second, b in before.items():
            assert pairing(a, b) == pairing(after[first], after[second])
    for name, (creation, annihilation) in after.items():
        if name.startswith('right'):
            creation = [value * Z ** (-order) for value, order in zip(creation, orders)]
            annihilation = [value * Z ** order for value, order in zip(annihilation, orders)]
        assert all(value.order() >= 0 for value in creation + annihilation if value)
    return {'interval': interval, 'ramification': ramification, 'core': core,
            'internal_gaussians_removed': internal_regular,
            'active_modes': len(active), 'spectator_modes': interval - len(active),
            'regular_shears': len(trace), 'canonical_pairings': canonical,
            'all_cross_pairings': len(before) ** 2}
