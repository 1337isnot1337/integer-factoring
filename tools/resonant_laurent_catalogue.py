#!/usr/bin/env python3
"""Exhaust three-term F4 Laurent resonance blocks by exact product classes."""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import product

import numpy as np

from f4_repaired_identity_filters import MUL


def polynomial(code, lower):
    return {lower+i: (code >> (2*i)) & 3 for i in range(3) if (code >> (2*i)) & 3}


def encoded_product(first, second):
    value = 0
    for i in range(3):
        for j in range(3):
            value ^= MUL[(first >> (2*i)) & 3][(second >> (2*j)) & 3] << (2*(i+j))
    return value


def first_index(code, digits):
    return next((i for i in range(digits) if (code >> (2*i)) & 3), 1000000)


def last_index(code, digits):
    return next((i for i in range(digits-1, -1, -1) if (code >> (2*i)) & 3), -1000000)


def endpoint(code, digits):
    value = 0
    for i in range(digits):
        value ^= (code >> (2*i)) & 3
    return value


@lru_cache(maxsize=2)
def offdiagonal_classes(lower):
    result = {}
    for first, second in product(range(64), repeat=2):
        key = (encoded_product(first, second),
               lower+min(first_index(first, 3), first_index(second, 3)),
               lower+max(last_index(first, 3), last_index(second, 3)))
        if key not in result:
            result[key] = [0, first, second]
        result[key][0] += 1
    assert len(result) == 607 and sum(value[0] for value in result.values()) == 4096
    return tuple((key, tuple(value)) for key, value in sorted(result.items()))


@lru_cache(maxsize=2)
def tables(lower):
    rows = offdiagonal_classes(lower)
    return {
        'rows': rows,
        'products': np.array([key[0] for key, value in rows], dtype=np.int64),
        'minimum': np.array([key[1] for key, value in rows], dtype=np.int64),
        'maximum': np.array([key[2] for key, value in rows], dtype=np.int64),
        'counts': np.array([value[0] for key, value in rows], dtype=np.int64),
        'multiply': np.array([[encoded_product(a, b) for b in range(64)] for a in range(64)], dtype=np.int64),
        'val3': np.array([first_index(value, 3) for value in range(64)], dtype=np.int64),
        'deg3': np.array([last_index(value, 3) for value in range(64)], dtype=np.int64),
        'val5': np.array([first_index(value, 5) for value in range(1024)], dtype=np.int64),
        'deg5': np.array([last_index(value, 5) for value in range(1024)], dtype=np.int64),
        'endpoint5': np.array([endpoint(value, 5) for value in range(1024)], dtype=np.int64),
    }


def diagonal_signature_groups(n, P, lower, first, second):
    """Return signatures, counts, and representative offdiagonals for two fixed diagonals."""
    if lower not in (-2, -1) or n*P < 4:
        raise ValueError('This exact band formula covers windows[-2,0] and[-1,1], with nP≥4.')
    m = n*P
    data = tables(lower)
    rows, products = data['rows'], data['products']
    min_off, max_off, counts = data['minimum'], data['maximum'], data['counts']
    constant_code = 1 << (-2*lower)
    delta_a, delta_b, trace = first ^ constant_code, second ^ constant_code, first ^ second
    determinant = data['multiply'][first, second] ^ products
    det_beta = data['multiply'][delta_a, delta_b] ^ products
    endpoint_values = data['endpoint5'][determinant]
    mask = endpoint_values != 0
    if not np.any(mask):
        return ()
    val_delta_a = lower+int(data['val3'][delta_a]) if delta_a else 2*m
    val_delta_b = lower+int(data['val3'][delta_b]) if delta_b else 2*m
    first_zero = np.minimum(min(val_delta_a, val_delta_b), min_off)
    val_trace = lower+int(data['val3'][trace]) if trace else 4*m
    second_band = 2*m+val_trace if trace else 4*m
    total_zero = np.where(det_beta != 0, 2*lower+data['val5'][det_beta], second_band)
    first_infinity = -2*np.maximum(max(lower+data['deg3'][first], lower+data['deg3'][second]), max_off)
    total_infinity = -2*(2*lower+data['deg5'][determinant])
    single_a = lower+int(data['val3'][delta_a]) if delta_a else 3*m
    single_b = lower+int(data['val3'][delta_b]) if delta_b else 3*m
    signatures = np.stack((np.full(len(rows), val_trace), np.full(len(rows), single_a),
                           np.full(len(rows), single_b), first_zero, total_zero-first_zero,
                           first_infinity, total_infinity-first_infinity), axis=1)[mask]
    valid_indices = np.flatnonzero(mask)
    unique, first_occurrence, inverse = np.unique(signatures, axis=0, return_index=True, return_inverse=True)
    weights = np.bincount(inverse, weights=counts[mask]).astype(np.int64)
    return tuple((tuple(map(int, value)), int(count), rows[int(valid_indices[index])][1][1:])
                 for value, count, index in zip(unique, weights, first_occurrence))


def complete_catalogue(n, P, lower):
    classes = {}
    for first, second in product(range(64), repeat=2):
        for values, count, offdiagonal in diagonal_signature_groups(n, P, lower, first, second):
            if values not in classes:
                classes[values] = {'signature': values, 'count': 0, 'window_lower': lower,
                                   'entry_codes': (first, offdiagonal[0], offdiagonal[1], second)}
            classes[values]['count'] += count
    total = sum(row['count'] for row in classes.values())
    assert total == 180*4**8 == 11796480
    return tuple(classes.values())
