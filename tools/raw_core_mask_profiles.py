#!/usr/bin/env python3
"""Arbitrary endpoint-preserving constant masks on the raw extreme core.

Physical labels are ordered (x,y,u,v). A mask has four entries per
label and satisfies m(x,y,y,x)=1. The original Gaussian endpoint is
triangular on typed words, so this condition preserves its crossing.
"""
from __future__ import annotations

from itertools import product

from blocked_core_charge_pressure import BlockMatrixOracle
from strict_core_charge_gauges import CoreChargePressure, gaussian_coordinates, projection_data


def index(labels, alphabet=4):
    value = 0
    for label in labels:
        value = alphabet * value + label
    return value


def normalize_mask(mask, alphabet=4):
    if callable(mask):
        mask = tuple(mask(*labels) for labels in product(range(alphabet), repeat=4))
    else:
        mask = tuple(mask)
    if len(mask) != alphabet ** 4 or any(value not in (0, 1, 2, 3) for value in mask):
        raise ValueError('Supply one F4 coefficient per physical core tuple.')
    if any(mask[index((x, y, y, x), alphabet)] != 1
           for x, y in product(range(alphabet), repeat=2)):
        raise ValueError('The mask must equal one on every retained diagonal tuple.')
    return mask


def pair_partition_mask(partition, alphabet=4):
    """Keep a transition when its two unordered core pairs share a class."""
    pairs = [(i, j) for i in range(alphabet) for j in range(i, alphabet)]
    if len(partition) != len(pairs):
        raise ValueError('Give one class per unordered pair, including repetitions.')
    classes = dict(zip(pairs, partition))
    return normalize_mask(lambda x, y, u, v:
        int(classes[tuple(sorted((x, u)))] == classes[tuple(sorted((y, v)))]), alphabet)


def additive_pair_partitions(primes=(3, 5, 7, 11, 13)):
    """Finite representatives only; no completeness claim over all groups."""
    found = {}
    for prime in primes:
        for rest in product(range(prime), repeat=3):
            charges = (0,) + rest
            first = next((value for value in charges if value), None)
            if first not in (None, 1):
                continue
            sums = [(charges[i] + charges[j]) % prime for i in range(4) for j in range(i, 4)]
            labels = {}
            signature = tuple(labels.setdefault(value, len(labels)) for value in sums)
            found.setdefault(signature, {'modulus': prime, 'charges': charges})
    return found


class CoreMaskPressure(CoreChargePressure):
    """The full physical objective with exact active constant-mask matrices."""

    def __init__(self, interval, mask, family='raw_extremes'):
        self.interval, self.core, self.q = interval, 2, 4
        self.mask = normalize_mask(mask)
        self.family = family
        self.records, self.oracles = [], []
        self.queries = 0
        for place in ('zero', 'infinity'):
            data = projection_data(place, interval, 2, family)
            gaussian, shifts, metadata, left, right = gaussian_coordinates(place, interval, 2, family)
            left_letters, right_letters = set(left['labels']), set(right['labels'])
            assert left_letters | right_letters == set('xyuv') and not left_letters & right_letters
            matrix = []
            for i, source in enumerate(gaussian):
                row = []
                for j, pair in enumerate(source):
                    if pair == (0, 0):
                        row.append(pair)
                        continue
                    labels = tuple((left['labels'][letter][i] if letter in left_letters
                                    else right['labels'][letter][j]) for letter in 'xyuv')
                    scalar = self.mask[index(labels)]
                    assert pair[1] == 0
                    row.append((pair[0] if scalar & 1 else 0, pair[0] if scalar & 2 else 0))
                matrix.append(row)
            self.records.append((data, metadata, left, right))
            self.oracles.append(BlockMatrixOracle(matrix, shifts, left['row_charges'],
                                                   right['column_charges'], metadata['profile_precision']))
        self.core_indices = self.records[0][0]['physical_core_indices']
        self.physical_spectators = [i for i in range(interval) if i not in self.core_indices]
