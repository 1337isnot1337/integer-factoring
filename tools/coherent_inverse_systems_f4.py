#!/usr/bin/env python3
"""Weighted homogeneous derivative inverse systems over F4.

This extends the repository's binary inverse-system constructor. Sparse
occupation vectors use the existing F4 encoding and ordinary bilinear
pairing; no conjugation is introduced. Leading occupation words form a
deletion downset, giving a compatible coordinate retraction.
"""
from __future__ import annotations

from coherent_hot_compounds import height
from coherent_hot_compounds_f4 import axpy, mul, validate_maps


def derivative(vector, deleted):
    return {word ^ deleted: scalar for word, scalar in vector.items()
            if word & deleted == deleted}


def inverse_system(top_forms, modes):
    top_forms = [dict(vector) for vector in top_forms]
    if not 1 <= modes <= 16 or not top_forms:
        raise ValueError('Use nonempty top forms on one through sixteen modes.')
    for vector in top_forms:
        if not vector or any(not isinstance(word, int) or not 0 <= word < 1 << modes
                             or scalar not in (1, 2, 3) for word, scalar in vector.items()):
            raise ValueError('Every top form must be a nonzero sparse F4 occupation vector.')
        if len({height(word) for word in vector}) != 1 or len({word.bit_count() for word in vector}) != 1:
            raise ValueError('Each top form must have one height and one occupation.')
    pivots, queue = {}, []

    def insert(vector):
        vector = dict(vector)
        while vector:
            pivot = max(vector)
            scalar = vector[pivot]
            if pivot in pivots:
                vector = axpy(vector, scalar, pivots[pivot])
            else:
                inverse = mul(scalar, scalar)
                value = {word: mul(inverse, coefficient) for word, coefficient in vector.items()}
                pivots[pivot] = value
                queue.append(value)
                return

    for vector in top_forms:
        insert(vector)
    # Derivatives of a basis generate the entire closure, so only independent
    # vectors enter the queue. There are at most 2^modes such vectors.
    for vector in queue:
        for bit in range(modes):
            insert(derivative(vector, 1 << bit))
    for pivot in sorted(pivots):
        for old in sorted(pivots):
            if old != pivot and pivot in pivots[old]:
                pivots[old] = axpy(pivots[old], pivots[old][pivot], pivots[pivot])
    words = sorted(pivots, key=lambda word: (height(word), word))
    assert all((word ^ (1 << bit)) in pivots for word in words
               for bit in range(modes) if word >> bit & 1)
    embedding = [pivots[word] for word in words]
    retraction = [{word: 1} for word in words]
    heights, charges = validate_maps(embedding, retraction, modes)
    return dict(modes=modes, field='F4: omega^2=omega+1', top_forms=top_forms,
                pivot_words=words, embedding=embedding, retraction=retraction,
                heights=heights, charges=charges, derivative_images=len(queue)*modes)
