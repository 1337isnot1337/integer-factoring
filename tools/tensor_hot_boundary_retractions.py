#!/usr/bin/env python3
"""Homogeneous low-cell retractions with every other height cell retained.

These maps are not required to be deletion closed. They must be evaluated
with the literal regular tensor, not the compressed-compound shortcut.
"""
from collections import defaultdict

import tensor_hot_compounds_f16 as field


def validate_height_maps(embedding, retraction, modes):
    if not 1 <= modes <= 16 or not embedding or len(embedding) != len(retraction):
        raise ValueError('Invalid physical-map dimensions.')
    heights = []
    for p, left in zip(embedding, retraction):
        if not p or not left or any(type(w) is not int or not 0 <= w < 1 << modes or not 0 < c < 16
                                   for vector in (p, left) for w, c in vector.items()):
            raise ValueError('Invalid sparse F16 vector.')
        hh = {field.height(w) for w in set(p) | set(left)}
        if len(hh) != 1:
            raise ValueError('The physical maps mix inherited heights.')
        heights.append(hh.pop())
    if any(field.dot(left, p) != int(i == j)
           for i, left in enumerate(retraction) for j, p in enumerate(embedding)):
        raise ValueError('LP=I failed.')
    diagnostics = {}
    for name, vectors in (('embedding', embedding), ('retraction', retraction)):
        pivots = field.basis(vectors)
        violations = []
        for i, vector in enumerate(vectors):
            for bit in range(modes):
                derivative = {w ^ (1 << bit): c for w, c in vector.items() if w >> bit & 1}
                if not field.in_span(derivative, pivots):
                    violations.append({'label': i, 'deleted_mode': bit, 'derivative': derivative})
        diagnostics[name] = {'deletion_closed': not violations, 'violations': violations}
    return heights, diagnostics


def boundary_code(height_cap, boundary_height, coefficients, modes=16, orthogonal=True):
    if not 0 <= boundary_height <= height_cap:
        raise ValueError('The altered height cell must be retained.')
    words = sorted((w for w in range(1 << modes) if field.height(w) <= height_cap),
                   key=lambda w: (field.height(w), w))
    boundary = [w for w in words if field.height(w) == boundary_height]
    if len(boundary) < 2 or len(coefficients) != len(boundary) - 1 or any(not 0 <= a < 16 for a in coefficients):
        raise ValueError('Supply one coefficient for each retained boundary word.')
    removed = boundary[-1]
    coefficients_by_word = dict(zip(boundary[:-1], coefficients))
    labels = [w for w in words if w != removed]
    embedding, left = [], []
    for word in labels:
        vector = {word: 1}
        coefficient = coefficients_by_word.get(word, 0)
        if coefficient:
            vector[removed] = coefficient
        embedding.append(vector)
        left.append({word: 1})
    if orthogonal:
        indices = [i for i, w in enumerate(labels) if field.height(w) == boundary_height]
        gram = [[field.dot(embedding[i], embedding[j]) for j in indices] for i in indices]
        inverse = field.inverse(gram)
        if inverse is None:
            raise ValueError('The boundary Gram matrix is singular.')
        for i, row in zip(indices, inverse):
            left[i] = {}
            for coefficient, j in zip(row, indices):
                left[i] = field.axpy(left[i], coefficient, embedding[j])
    heights, diagnostics = validate_height_maps(embedding, left, modes)
    return {'modes': modes, 'height_cap': height_cap, 'boundary_height': boundary_height,
            'boundary_words': boundary, 'removed_pivot_word': removed,
            'coefficients': coefficients, 'q': len(labels), 'ambient_dimension': len(words),
            'pivot_words': labels, 'embedding': embedding, 'retraction': left,
            'retraction_kind': 'orthogonal' if orthogonal else 'coordinate',
            'heights': heights, 'deletion_diagnostics': diagnostics}
