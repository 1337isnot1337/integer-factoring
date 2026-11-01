#!/usr/bin/env python3
"""Height-fifteen graph kernels with an F16 cubic differential.

The physical space is ker(D_7 + a D_0D_1D_5 + b D_0D_2D_4
                         + c D_1D_2D_3), truncated by inherited height.
The graph coordinates omit bit7. Occupation number may change by two.
"""
from collections import defaultdict

import tensor_hot_compounds_f16 as field


def validate_maps(embedding, retraction, modes):
    if not 1 <= modes <= 16 or not embedding or len(embedding) != len(retraction):
        raise ValueError('Invalid mixed-map dimensions.')
    heights, parities = [], []
    for p, left in zip(embedding, retraction):
        if not p or not left or any(type(w) is not int or not 0 <= w < 1 << modes or not 0 < c < 16
                                   for vector in (p, left) for w, c in vector.items()):
            raise ValueError('Invalid sparse F16 vector.')
        words = set(p) | set(left)
        hh = {field.height(w) for w in words}
        pp = {w.bit_count() % 2 for w in words}
        if len(hh) != 1 or len(pp) != 1:
            raise ValueError('Maps must preserve height and occupation parity.')
        heights.append(hh.pop())
        parities.append(pp.pop())
    if any(field.dot(left, p) != int(i == j)
           for i, left in enumerate(retraction) for j, p in enumerate(embedding)):
        raise ValueError('LP=I failed.')
    for vectors in (embedding, retraction):
        pivots = field.basis(vectors)
        if len(pivots) != len(embedding):
            raise ValueError('Dependent map.')
        for vector in vectors:
            for bit in range(modes):
                derivative = {w ^ (1 << bit): c for w, c in vector.items() if w >> bit & 1}
                if not field.in_span(derivative, pivots):
                    raise ValueError('Deletion closure failed.')
    return heights, parities


def graph_code(coefficients, height_cap, modes=16):
    if len(coefficients) != 3 or any(not 0 <= x < 16 for x in coefficients):
        raise ValueError('Expected three F16 cubic coefficients.')
    if not 8 <= modes <= 16 or type(height_cap) is not int or height_cap < 0:
        raise ValueError('Invalid graph dimensions.')
    cubic = {mask: coefficient for mask, coefficient in zip((35, 21, 14), coefficients) if coefficient}
    words = sorted((w for w in range(1 << modes) if not (w >> 7) & 1 and field.height(w) <= height_cap),
                   key=lambda w: (field.height(w), w))
    embedding = []
    for word in words:
        vector = {word: 1}
        vector.update({target | 128: c for target, c in field.apply_operator({word: 1}, cubic).items()})
        embedding.append(vector)
    left = [{w: 1} for w in words]
    heights, parities = validate_maps(embedding, left, modes)
    operator = {128: 1, **cubic}
    assert all(not field.apply_operator(vector, operator) for vector in embedding)
    return {'modes': modes, 'height_cap': height_cap, 'q': len(words), 'coefficients': coefficients,
            'operator': operator, 'embedding': embedding, 'retraction': left,
            'pivot_words': words, 'heights': heights, 'parities': parities}


def orthogonal_retraction(code):
    embedding = code['embedding']
    cells = defaultdict(list)
    for i, height in enumerate(code['heights']):
        cells[height].append(i)
    left = [{} for _ in embedding]
    for indices in cells.values():
        gram = [[field.dot(embedding[i], embedding[j]) for j in indices] for i in indices]
        inverse = field.inverse(gram)
        if inverse is None:
            return None
        for i, row in zip(indices, inverse):
            for coefficient, j in zip(row, indices):
                left[i] = field.axpy(left[i], coefficient, embedding[j])
    validate_maps(embedding, left, code['modes'])
    return left
