#!/usr/bin/env python3
"""Exact fixed-point-row jets for canonical Boolean Cauchy packets.

Unsigned wraparound preserves residues modulo 2^precision for precision <=63.
No floating-point arithmetic is used.
"""
from functools import lru_cache
import numpy as np
from cauchy_packet_polynomials import cauchy_polynomials, shift

@lru_cache(maxsize=8)
def _local_kernels(length, precision):
    modulus = 1 << precision
    kernels = cauchy_polynomials(length)
    return tuple({q: np.array([v.numerator*pow(v.denominator,-1,modulus)%modulus
                              for v in shift(poly,center,2*length)],dtype=np.uint64)
                  for q,poly in kernels.items()}
                 for center in range(-length+1,length,2))

def fixed_row_jets(length, row_state, precision):
    modulus, mask = 1 << precision, np.uint64((1 << precision)-1)
    if not isinstance(precision,int) or not 1<=precision<=63:
        raise ValueError("Require precision between 1 and 63")
    if not isinstance(row_state,int) or not 0<=row_state<(1<<length):
        raise ValueError("Require a Boolean row state")
    degree, count = length*(2*length-1), 1 << length
    blocks = []
    for local in _local_kernels(length,precision):
        block = np.zeros((degree+1,count),dtype=np.uint64)
        for state in range(count):
            polynomial = np.array([1],dtype=np.uint64)
            for x in range(length):
                image = sum((local[y-x] for y in range(length) if state>>y&1),
                            np.zeros(2*length,dtype=np.uint64)) & mask
                if not row_state>>x&1:
                    image = -image
                    image[0] += np.uint64(1)
                polynomial = np.convolve(polynomial,image) & mask
            block[:,state] = polynomial
        blocks.append(block)
    return np.concatenate(blocks)

