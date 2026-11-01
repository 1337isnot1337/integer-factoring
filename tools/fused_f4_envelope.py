#!/usr/bin/env python3
"""Exact F4 scalar extension of the literal ordered two-copy envelope.

Coordinates are two binary bit planes in the retained 33-vector basis.
No rank-sum improvement is asserted by this representation.
"""
from __future__ import annotations

from check_fused_scalar_projection import (TARGET, coefficient_matrices, coefficient_space,
                                         combine, partitions, projected_tensor)
from f4_lattice_linear import INV, add, coefficient, eliminate, scale
from joint_energy_conjugations import Space


def normalize(coordinates):
    if not coordinates[0] | coordinates[1]:
        return 0, 0
    pivot = (coordinates[0] | coordinates[1]).bit_length()-1
    return scale(coordinates, INV[coefficient(coordinates, pivot)])


class Envelope:
    def __init__(self):
        self.cuts, self.labels = coefficient_space(ordered=True)
        assert len(self.labels) == 33
        self.groups = partitions()
        for basis, groups in zip(self.cuts, self.groups):
            masks = {r: sum(1 << c for c in group) for group in groups for r in group}
            for matrix in basis:
                assert all(((matrix >> (81*r)) & ((1 << 81)-1)) & ~masks[r] == 0
                           for r in range(81))
        self.coordinate_space = Space(self.cuts[0])
        self.binary_cache = {}
        self.rank_cache = {}

    def matrices(self, coordinates):
        low, high = coordinates
        if low < 0 or high < 0 or (low | high) >> 33:
            raise ValueError('Use two nonnegative 33-bit coordinate planes.')
        for value in (low, high):
            if value not in self.binary_cache:
                self.binary_cache[value] = tuple(combine(basis, value) for basis in self.cuts)
        return tuple((self.binary_cache[low][side], self.binary_cache[high][side])
                     for side in range(2))

    def ranks(self, coordinates):
        coordinates = normalize(coordinates)
        if coordinates not in self.rank_cache:
            ranks = []
            for (low, high), groups in zip(self.matrices(coordinates), self.groups):
                rank = 0
                for group in groups:
                    rows = [(sum(((low >> (81*r+c)) & 1) << j for j, c in enumerate(group)),
                             sum(((high >> (81*r+c)) & 1) << j for j, c in enumerate(group)))
                            for r in group]
                    rank += len(eliminate(rows))
                ranks.append(rank)
            self.rank_cache[coordinates] = tuple(ranks)
        return self.rank_cache[coordinates]

    def projected_piece(self, i, j, scalar=1):
        # In F16=F4[u], the constant F4 coefficient has binary digits
        # chi(z),chi(omega*z). This is a projection, not the field trace.
        from f4_lattice_linear import f16_multiply
        row, column = coefficient_matrices()
        result = []
        for multiplier in (scalar, f16_multiply(2, scalar)):
            matrix = projected_tensor(row[i], row[j], multiplier)
            coordinates = self.coordinate_space.binary_coordinates(matrix)
            assert coordinates is not None
            assert combine(self.cuts[1], coordinates) == projected_tensor(column[i], column[j], multiplier)
            result.append(coordinates)
        return tuple(result)

    def decomposition(self, pieces):
        total = (0, 0)
        ranks = []
        for piece in pieces:
            total = add(total, piece)
            ranks.append(self.ranks(piece))
        if total != (TARGET, 0):
            raise ValueError('The proposed pieces do not sum to the crossing.')
        return ranks
