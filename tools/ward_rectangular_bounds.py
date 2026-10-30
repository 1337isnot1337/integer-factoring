#!/usr/bin/env python3
"""Complete rectangular Pfaffian bounds and the Ward row-grade descent."""

from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import comb, isqrt, log2

from check_modular_crossing import rank_two
from stable_pfaffian_bounds import stable_multiplicities, weighted_masks, centered_row_bounds


GRADES = ((1, 0), (0, 1), (2, 1), (1, 2))
ANCILLAS = ((1, 0), (1, 1), (1, 2), (1, 3), (0, 0), (0, 2))


@lru_cache(None)
def rectangular_block(m, n, ra, ry, cv, cu, pairs, ar, ac):
    a, b = GRADES[ar]
    c, d = GRADES[ac]
    numerators = (ry + b - ra - a + m, cu + d - cv - c + m,
                  ra + ry + a + b - 2 * pairs - 2 + n,
                  cv + cu + c + d - 2 * pairs - 2 + n)
    if any(value % 2 for value in numerators):
        return (), 0
    ur, uc, vr, vc = (value // 2 for value in numerators)
    fr, fs = weighted_masks(m, ur)
    fc, _ = weighted_masks(m, uc)
    sr, ss = weighted_masks(n, vr)
    sc, _ = weighted_masks(n, vc)
    if not fr or not fc or not sr or not sc:
        return (), 0
    first = [[fs[x & y] for y in fc] for x in fr]
    second = [[ss[x & y] for y in sc] for x in sr]
    stride = m + n + 5
    values = {}
    for w in {v for row in first for v in row}:
        for z in {v for row in second for v in row}:
            assert min(m - ur - uc + w, uc - w, vr - z, z) >= 0
            value = 1
            for which, exponent in enumerate((uc - w, vr - z, z)):
                while exponent:
                    bit = exponent & -exponent
                    if which == 0:
                        value ^= value << bit
                    elif which == 1:
                        value ^= value << (stride * bit)
                    else:
                        value ^= (value << bit) ^ (value << (stride * bit))
                    exponent ^= bit
            values[w, z] = value
    rdim, cdim = len(fr) * len(sr), len(fc) * len(sc)
    rows = {}
    for i, k in product(range(len(fr)), range(len(sr))):
        ri = i * len(sr) + k
        for j, h in product(range(len(fc)), range(len(sc))):
            cj = j * len(sc) + h
            value = values[first[i][j], second[k][h]]
            while value:
                bit = value & -value
                degree = bit.bit_length() - 1
                key = degree % stride, degree // stride
                rows.setdefault(key, [0] * rdim)[ri] ^= 1 << cj
                value ^= bit
    return tuple((a, b, rank_two(matrix)) for (a, b), matrix in sorted(rows.items())), max(rdim, cdim)


@lru_cache(None)
def rectangular_bounds(pairs, ra, ry, cv, cu):
    if min(ra, ry, cv, cu) < 0 or max(ra, ry, cv, cu) > 2 * pairs:
        return {}, 0
    cap = min(comb(2 * pairs, ra) * comb(2 * pairs, ry),
              comb(2 * pairs, cv) * comb(2 * pairs, cu))
    candidates, largest = [], 0
    for ar, ac in ANCILLAS:
        rows = defaultdict(int)
        for m, n, twist, count in stable_multiplicities(pairs):
            profile, size = rectangular_block(m, n, ra, ry, cv, cu, pairs, ar, ac)
            largest = max(largest, size)
            for a, b, rank in profile:
                if a + twist + b <= ra:
                    rows[a + twist, b] += count * rank
        candidates.append(rows)
    return {(a, b): min(cap, *(row.get((a, b), 0) for row in candidates))
            for a in range(ra + 1) for b in range(ra + 1 - a)}, largest


def ward_row_bounds(pairs, length, forms):
    answer = defaultdict(int)
    largest = 0
    for j in range(forms + 1):
        for h in range(forms - j + 1):
            ra, ry = length - forms - j, length - 3 * forms + j + 2 * h
            profile, size = rectangular_bounds(pairs, ra, ry, length, length)
            largest = max(largest, size)
            for (a, b), rank in profile.items():
                answer[a, b + forms + j] += rank
    return dict(answer), largest


def ward_column_bounds(pairs, length, forms, channel_counts=None):
    answer = defaultdict(int)
    largest = 0
    for j in range(forms + 1):
        ell = forms - j
        channels = (2 ** ell * comb(pairs, ell) if channel_counts is None
                    else channel_counts[ell])
        for h in range(ell + 1):
            ra, ry = length - forms - j, length - 3 * forms + j + 2 * h
            profile, size = rectangular_bounds(pairs, ra, length, length, ry)
            largest = max(largest, size)
            for (b, a), rank in profile.items():
                answer[a, b + forms + j] += channels * rank
    return dict(answer), largest


def ward_score(pairs, length, forms, labels, channel_counts=None):
    ward, largest = ward_row_bounds(pairs, length, forms)
    column, column_largest = ward_column_bounds(pairs, length, forms, channel_counts)
    original, _ = rectangular_bounds(pairs, length, length, length, length)
    total, profile = 0, []
    for a in range(length + 1):
        for b in range(length + 1 - a):
            if b < forms:
                row = col = 0
            else:
                row = min(labels * labels, ward.get((a, b), 0), original.get((a, b), 0))
                col = min(labels * labels, original.get((b, a), 0), column.get((a, b), 0))
            if row or col:
                value = row * col
                root = isqrt(value)
                total += root + int(root * root < value)
                profile.append([a, b, row, col])
    return {'pairs': pairs, 'length': length, 'forms': forms, 'labels': labels,
            'integer_rank_sum_upper_bound': total, 'profile': profile,
            'largest_explicit_block': max(largest, column_largest),
            'rho_of_upper_bound': log2(total) / (2 * log2(labels))}
