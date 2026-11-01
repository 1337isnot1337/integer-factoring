#!/usr/bin/env python3
"""Nondegenerate bilinear Gram pairings for homogeneous F4 physical codes."""
from __future__ import annotations

import time

from coherent_hot_compounds_f4 import axpy, compound_block, dot, mul, tensor_vector, validate_maps


def paired_retraction(embedding, dual_space, modes):
    size = len(embedding)
    if len(dual_space) != size:
        raise ValueError('The embedding and dual-space bases must have equal size.')
    gram = [[dot(a, b) for b in embedding] for a in dual_space]
    augmented = [row+[int(i == j) for j in range(size)] for i, row in enumerate(gram)]
    for c in range(size):
        pivot = next((r for r in range(c, size) if augmented[r][c]), None)
        if pivot is None:
            raise ValueError('The ordinary bilinear pairing matrix is singular.')
        augmented[c], augmented[pivot] = augmented[pivot], augmented[c]
        inverse = mul(augmented[c][c], augmented[c][c])
        augmented[c] = [mul(inverse, value) for value in augmented[c]]
        for r in range(size):
            if r != c and augmented[r][c]:
                scalar = augmented[r][c]
                augmented[r] = [a ^ mul(scalar, b) for a, b in zip(augmented[r], augmented[c])]
    inverse = [row[size:] for row in augmented]
    retraction = []
    for row in inverse:
        vector = {}
        for scalar, source in zip(row, dual_space):
            vector = axpy(vector, scalar, source)
        retraction.append(vector)
    heights, charges = validate_maps(embedding, retraction, modes)
    return dict(embedding=embedding, retraction=retraction, gram=gram, inverse_gram=inverse,
                heights=heights, charges=charges)


def orthogonal_retraction(embedding, modes):
    result = paired_retraction(embedding, embedding, modes)
    assert all(not result['gram'][i][j]
               or (result['heights'][i], result['charges'][i]) == (result['heights'][j], result['charges'][j])
               for i in range(len(embedding)) for j in range(len(embedding)))
    return result


def orthogonal_profiles(embedding, modes, keep_blocks=False):
    """Compute one cut; both regular cuts are equivalent by constant Gram units."""
    code = orthogonal_retraction(embedding, modes)
    retraction, charges = code['retraction'], code['charges']
    size = len(embedding)
    started, orders, blocks = time.monotonic(), [], []
    for charge in sorted({a+b for a in charges for b in charges}):
        pairs = [(x, y) for x in range(size) for y in range(size) if charges[x]+charges[y] == charge]
        rows = [tensor_vector(retraction[x], retraction[y], modes) for x, y in pairs]
        columns = [tensor_vector(embedding[y], embedding[x], modes) for x, y in pairs]
        block = compound_block(modes, rows, columns, keep_blocks)
        blocks.append(block)
        orders += block['orders_t']
    assert len(orders) == size*size
    orders.sort()
    return dict(q=size, modes=modes, heights=code['heights'], charges=charges,
                zero_orders_t=orders, infinity_orders_t=orders,
                cuts=[dict(cut='zero', orders_t=orders, blocks=blocks)],
                infinity_derived_by_constant_Gram_equivalence=True,
                seconds=time.monotonic()-started)
