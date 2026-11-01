#!/usr/bin/env python3
"""Homogeneous deletion-closed codes defined by differential equations.

All bit derivatives commute. Therefore the common kernel of any collection
of homogeneous Boolean differential operators, inside a height downset, is
deletion closed. A reduced leading-word basis gives a coordinate retract.
"""
from __future__ import annotations

from collections import defaultdict

from coherent_hot_compounds import height, support, validate_maps
from coherent_inverse_systems import derivative


def differential_kernel(operators, modes, height_cap):
    operators = tuple(operators)
    if not 1 <= modes <= 16 or not operators or not isinstance(height_cap, int) or height_cap < 0:
        raise ValueError('Use homogeneous nonzero operators, up to sixteen modes, and a nonnegative cap.')
    for operator in operators:
        if not isinstance(operator, int) or operator <= 0 or operator.bit_length() > 1 << modes:
            raise ValueError('Each differential operator must be a nonzero binary occupation vector.')
        words = support(operator)
        if (len({height(word) for word in words}) != 1
                or len({word.bit_count() for word in words}) != 1):
            raise ValueError('Each operator must have one height and one differential degree.')
    cells = defaultdict(list)
    for word in range(1 << modes):
        if height(word) <= height_cap:
            cells[height(word), word.bit_count()].append(word)
    ambient_width = 1 << modes
    kernel_generators, ambient_count = [], 0
    for words in cells.values():
        image_pivots = {}
        for word in words:
            ambient_count += 1
            preimage, image = 1 << word, 0
            for index, operator in enumerate(operators):
                value = 0
                for removed in support(operator):
                    value ^= derivative(preimage, removed)
                image ^= value << (index * ambient_width)
            while image:
                pivot = image.bit_length() - 1
                if pivot in image_pivots:
                    old_image, old_preimage = image_pivots[pivot]
                    image ^= old_image
                    preimage ^= old_preimage
                else:
                    image_pivots[pivot] = image, preimage
                    break
            if image == 0:
                kernel_generators.append(preimage)
    pivots = {}
    for value in kernel_generators:
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    for pivot in sorted(pivots):
        for other in pivots:
            if other != pivot and pivots[other] >> pivot & 1:
                pivots[other] ^= pivots[pivot]
    words = sorted(pivots, key=lambda word: (height(word), word))
    embedding = [pivots[word] for word in words]
    retraction = [1 << word for word in words]
    _, _, heights, charges, checks = validate_maps(embedding, retraction, modes)
    for vector in embedding:
        for operator in operators:
            value = 0
            for removed in support(operator):
                value ^= derivative(vector, removed)
            assert value == 0
    return {'modes': modes, 'height_cap': height_cap, 'operators': list(operators),
            'ambient_dimension': ambient_count, 'q': len(embedding),
            'embedding': embedding, 'retraction': retraction, 'pivot_words': words,
            'heights': list(heights), 'charges': list(charges), 'deletion_checks': checks}
