#!/usr/bin/env python3
"""Exact four-copy Pauli crossing compressed by a ternary CSS split.

All polynomial coefficients survive. The code is the coset incidence of
<1111> in the kernel of 1122 over F3, with transpose left inverse over F4.
"""
from __future__ import annotations

from hashlib import sha256
from itertools import product

import numpy as np

from cyclic_pauli_crossings import construct, literal_seed


def code():
    a, b = np.ones(4, dtype=np.int64), np.array([1, 1, 2, 2])
    representatives = np.array([s for s in product(range(3), repeat=4)
                                if s[0] == 0 and np.dot(s, b) % 3 == 0])
    words = (representatives[:, None, :] + np.arange(3)[None, :, None]*a) % 3
    indices = words @ (3**np.arange(3, -1, -1))
    embedding = np.zeros((81, 9), dtype=np.uint8)
    for j, row in enumerate(indices):
        embedding[row, j] = 1
    assert np.array_equal((embedding.T.astype(int) @ embedding) % 2, np.eye(9))
    return words, embedding


def binary_multiply(a, b, degree):
    result = np.zeros_like(a, dtype=np.uint64)
    for j in range(degree+1):
        result ^= (a << np.uint64(j))*((b >> np.uint64(j)) & np.uint64(1))
    return result


def polynomial_multiply(left, right, right_degree=3):
    a, b = left
    c, d = right
    bd = binary_multiply(b, d, right_degree)
    return [binary_multiply(a, c, right_degree) ^ bd,
            binary_multiply(a, d, right_degree) ^ binary_multiply(b, c, right_degree) ^ bd]


def compress():
    """Return all 9^4 entries in (x,y)|(v,u) order, degree at most 12."""
    words, _ = code()
    cell, cell_metadata = construct(1)
    # Each logical quadruple has 3^4 literal summands from independent legs.
    logical = np.array(list(product(range(9), repeat=4)), dtype=np.int64)
    shifts = np.array(list(product(range(3), repeat=4)), dtype=np.int64)
    physical = [words[logical[:, leg, None], shifts[None, :, leg], :]
                for leg in range(4)]
    result = [np.ones((9**4, 3**4), dtype=np.uint64),
              np.zeros((9**4, 3**4), dtype=np.uint64)]
    for i in range(4):
        x, y, v, u = (s[:, :, i] for s in physical)
        local = [p[3*x+y, 3*v+u] for p in cell]
        result = polynomial_multiply(result, local, cell_metadata['degree'])
    planes = [np.bitwise_xor.reduce(p, axis=1).reshape(81, 81) for p in result]
    from cyclic_pauli_crossings import parity
    assert np.array_equal(parity(planes[0]), np.eye(81))
    assert not np.any(parity(planes[1]))
    assert all(not np.any(p >> np.uint64(13)) for p in planes)
    metadata = {'q': 9, 'degree': 12, 'physical_length': 4, 'logical_length': 2,
                'a': [1, 1, 1, 1], 'b': [1, 1, 2, 2],
                'code_words': words.tolist(), 'cell': cell_metadata,
                'planes_sha256': sha256(b''.join(p.astype('<u8').tobytes() for p in planes)).hexdigest()}
    return planes, metadata


def seed():
    planes, metadata = compress()
    return literal_seed(planes, metadata), planes, metadata
