#!/usr/bin/env python3
"""Exact sufficient all-scale comparison of products of profile sums.

If every ascending candidate convolution prefix sum is at most its
reference counterpart, then the candidate sum of exp(-s*x) is at least
the reference sum for every s>0. Thus the candidate cannot improve the
profile-product cost at unchanged physical energies.
"""
from __future__ import annotations


def convolution(profiles):
    if len(profiles) != 2 or any(not row for row in profiles):
        raise ValueError('Use two nonempty complete order lists.')
    if any(not isinstance(value, int) for row in profiles for value in row):
        raise ValueError('This exact comparison requires integer orders.')
    return tuple(sorted(a+b for a in profiles[0] for b in profiles[1]))


def dominated(profiles, reference_convolution):
    candidate = convolution(profiles)
    if len(candidate) != len(reference_convolution):
        raise ValueError('Compare convolution spectra of the same cardinality.')
    difference = 0
    for a, b in zip(candidate, reference_convolution):
        difference += a-b
        if difference > 0:
            return False
    return True

