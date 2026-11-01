#!/usr/bin/env python3
"""Homogeneous coherent physical codes from Boolean derivative inverse systems.

The derivative closure of homogeneous top forms has a homogeneous reduced
basis. Its leading occupation masks form a bit-deletion downset. Coordinate
projection onto those masks is a compatible deletion-closed retraction.
"""
from __future__ import annotations

from coherent_hot_compounds import height, support, validate_maps


def derivative(vector, deleted):
    result = 0
    for word in support(vector):
        if word & deleted == deleted:
            result ^= 1 << (word ^ deleted)
    return result


def inverse_system(top_forms, modes):
    top_forms = tuple(top_forms)
    if not 1 <= modes <= 16 or not top_forms:
        raise ValueError('Use nonempty top forms on one through sixteen modes.')
    for value in top_forms:
        if not isinstance(value, int) or value <= 0 or value.bit_length() > 1 << modes:
            raise ValueError('Every top form must be a nonzero binary occupation vector.')
        words = support(value)
        if len({height(w) for w in words}) != 1 or len({w.bit_count() for w in words}) != 1:
            raise ValueError('Each top form must have one height and one occupation.')
    pivots = {}
    derivative_count = 0
    for top in top_forms:
        used = 0
        for word in support(top):
            used |= word
        deleted = used
        while True:
            value = derivative(top, deleted)
            derivative_count += 1
            while value:
                pivot = value.bit_length() - 1
                if pivot in pivots:
                    value ^= pivots[pivot]
                else:
                    pivots[pivot] = value
                    break
            if deleted == 0:
                break
            deleted = (deleted - 1) & used
    for pivot in sorted(pivots):
        for old in pivots:
            if old != pivot and pivots[old] >> pivot & 1:
                pivots[old] ^= pivots[pivot]
    words = sorted(pivots, key=lambda w: (height(w), w))
    embedding, retraction = tuple(pivots[w] for w in words), tuple(1 << w for w in words)
    _, _, heights, charges, deletion_checks = validate_maps(embedding, retraction, modes)
    return {'modes': modes, 'top_forms': list(top_forms), 'pivot_words': words,
            'embedding': list(embedding), 'retraction': list(retraction),
            'heights': list(heights), 'charges': list(charges),
            'derivative_generators': derivative_count, 'deletion_checks': deletion_checks}
