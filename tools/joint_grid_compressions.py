#!/usr/bin/env python3
"""Joint coordinate pair compression with a new additive physical grid."""
from collections import Counter, defaultdict
from itertools import combinations_with_replacement, permutations, product

from check_joint_energy_conjugations import polynomial
from check_triangular_auxiliary_row_fiber import HEIGHT
from joint_height_conjugations import source_matrix
from weighted_smith_oracle_f4 import WeightedSmithOracle


def energy_grids(qx, qy):
    """Sorted nonnegative grids with one vacuum and available total heights."""
    available = Counter(a+b for a in HEIGHT for b in HEIGHT)
    for first in combinations_with_replacement(range(1, 9), qx-1):
        for second in combinations_with_replacement(range(1, 9), qy-1):
            hx, hy = (0,)+first, (0,)+second
            needed = Counter(a+b for a in hx for b in hy)
            if all(count <= available[h] for h, count in needed.items()):
                yield hx, hy


def coordinate_maps(hx, hy):
    fibers, target = defaultdict(list), defaultdict(list)
    for a, b in product(range(4), repeat=2):
        fibers[HEIGHT[a]+HEIGHT[b]].append(4*a+b)
    for x, y in product(range(len(hx)), range(len(hy))):
        target[hx[x]+hy[y]].append(len(hy)*x+y)
    levels = sorted(target)
    choices = [tuple(permutations(fibers[h], len(target[h]))) for h in levels]
    for selected in product(*choices):
        mapping = [None]*(len(hx)*len(hy))
        for h, images in zip(levels, selected):
            for position, image in zip(target[h], images):
                mapping[position] = image
        yield tuple(mapping)


def crossing(hx, hy, mapping):
    qx, qy = len(hx), len(hy)
    if len(mapping) != qx*qy or len(set(mapping)) != len(mapping):
        raise ValueError('Use an injective joint coordinate map of the full new grid.')
    if any(not 0 <= image < 16 or HEIGHT[image//4]+HEIGHT[image%4] != hx[x]+hy[y]
           for (x, y), image in zip(product(range(qx), range(qy)), mapping)):
        raise ValueError('The joint map must preserve the sum of physical heights.')
    source, values = source_matrix(), []
    zero = [[0]*(qx*qy) for _ in range(qx*qy)]
    infinity = [[0]*(qx*qy) for _ in range(qx*qy)]
    for x, y, u, v in product(range(qx), range(qy), range(qy), range(qx)):
        row, column = mapping[x*qy+y], mapping[v*qy+u]
        value = source[row][4*(column%4)+column//4]
        assert value.bit_count() % 2 == int(x == v and y == u)
        values.append(value)
        normalized = value << (hx[x]+hy[y])
        assert normalized & 15 == 0
        zero[x*qy+y][v*qy+u] = polynomial((normalized >> 4, 0))
        bound, reversed_value = hx[v]+hy[u]+4, 0
        for degree in range(value.bit_length()):
            if value >> degree & 1:
                assert degree <= bound
                reversed_value ^= 1 << (bound-degree)
        infinity[x*qy+u][v*qy+y] = polynomial((reversed_value, 0))
    return dict(q_x=qx, q_y=qy, heights_x=list(hx), heights_y=list(hy),
                mapping=list(mapping), values=values, zero=zero, infinity=infinity)


def profiles(code):
    size = code['q_x']*code['q_y']
    return {side+'_orders': list(WeightedSmithOracle(code[side]).profile([0]*size, [0]*size)[0])
            for side in ('zero', 'infinity')}
