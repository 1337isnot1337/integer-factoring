#!/usr/bin/env python3
"""Choose equal-order Smith channel bases by exact matroid optimization.

This extends adaptive_smith_lattices: a finite state set is closed under
every nonzero binary combination inside each equal-order channel block.
A positive rational potential selects a minimum-cost basis in each block.
Every selected physical factorization is reconstructed entry by entry.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import product

import adaptive_smith_lattices as adaptive
import check_cancellation_module_automata as arithmetic
import check_exact_module_orbits as lattices
import check_minimal_boundary_modules as linear


ZERO, ONE = arithmetic.ZERO, arithmetic.ONE


def capped_lattice(generators, dimension, precision):
    residues = arithmetic.module_basis(
        [tuple(arithmetic.jet(value, precision) for value in row)
         for row in generators], dimension, precision)
    lifted = arithmetic.lift_generators(residues, dimension, precision)
    return lattices.hermite_lattice(
        tuple(tuple(arithmetic.rat(value) for value in row) for row in lifted),
        dimension)


def binary_plane_lattices(precision):
    """Every normalized full-tail lattice in O^2, in a fixed order."""
    assert precision >= 1
    result = [((ONE, ZERO), (ZERO, ONE))]
    for depth in range(1, precision + 1):
        for residue in range(1 << depth):
            result.append(lattices.hermite_lattice(
                ((ONE, arithmetic.rat(residue)),
                 (ZERO, arithmetic.rat(1 << depth))), 2))
        for residue in range(0, 1 << depth, 2):
            result.append(lattices.hermite_lattice(
                ((arithmetic.rat(1 << depth), ZERO),
                 (arithmetic.rat(residue), ONE)), 2))
    assert len(result) == len(set(result)) == 3 * (1 << precision) - 2
    return tuple(result)


def combine(mask, vectors):
    result = [ZERO] * len(vectors[0])
    for index, vector in enumerate(vectors):
        if mask >> index & 1:
            result = [arithmetic.radd(left, right)
                      for left, right in zip(result, vector)]
    return tuple(result)


def binary_inverse(rows):
    size = len(rows)
    assert all(0 < row < (1 << size) for row in rows)
    work = [row | (1 << (size + index)) for index, row in enumerate(rows)]
    for column in range(size):
        pivot = next(index for index in range(column, size)
                     if work[index] >> column & 1)
        work[column], work[pivot] = work[pivot], work[column]
        for index in range(size):
            if index != column and work[index] >> column & 1:
                work[index] ^= work[column]
    assert all(row & ((1 << size) - 1) == (1 << index)
               for index, row in enumerate(work))
    inverse = tuple(row >> size for row in work)
    for i, row in enumerate(rows):
        value = 0
        for j in range(size):
            if row >> j & 1:
                value ^= inverse[j]
        assert value == 1 << i
    return inverse


def greedy_basis(candidates, rank, potential):
    """Minimum-weight basis of the binary vector matroid, with exact costs."""
    basis, chosen = [], []
    for mask, target in sorted(candidates, key=lambda item: (potential[item[1]], item[0])):
        enlarged = arithmetic.binary_basis(basis + [mask])
        if len(enlarged) > len(basis):
            basis = list(enlarged)
            chosen.append((mask, target))
            if len(basis) == rank:
                break
    assert len(basis) == rank
    return tuple(chosen)


class MatroidSmithLattices:
    """Full-tail states closed under equal-order channel combinations."""

    def __init__(self, local, degree, *, alphabet=2, cut=0, side=0,
                 infinity=False, initial=None, terminal=None, precision=2,
                 seed_states=None, maximum_states=1000, maximum_block_rank=16):
        dimension = len(next(iter(local.values())))
        self.dimension, self.alphabet, self.degree = dimension, alphabet, degree
        self.cut, self.side, self.infinity = cut, side, infinity
        self.precision = precision
        assert alphabet >= 2 and degree >= 0 and precision >= 1
        assert cut in (0, 1) and side in (0, 1)
        assert set(local) == set(product(range(alphabet), repeat=4))
        assert all(len(matrix) == dimension and all(len(row) == dimension for row in matrix)
                   for matrix in local.values())
        assert all(isinstance(value, int) and 0 <= value < (1 << (degree + 1))
                   for matrix in local.values() for row in matrix for value in row)
        self.initial = tuple(initial or (ONE,) + (ZERO,) * (dimension - 1))
        self.terminal = tuple(terminal or (ONE,) + (ZERO,) * (dimension - 1))
        assert len(self.initial) == len(self.terminal) == dimension
        assert all(value in (ZERO, ONE) for value in self.initial + self.terminal)
        assert any(value != ZERO for value in self.initial)
        first = capped_lattice((self.initial,), dimension, precision)
        self.states = list(seed_states) if seed_states is not None else [first]
        assert len(self.states) == len(set(self.states))
        assert all(state == capped_lattice(state, dimension, precision)
                   for state in self.states)
        assert all(min(arithmetic.rorder(value) for row in state for value in row
                       if value != ZERO) == 0 for state in self.states)
        assert first in self.states
        indices = {state: index for index, state in enumerate(self.states)}
        self.start = indices[first]
        self.initial_coordinates = lattices.membership(self.initial, first)
        self.blocks, self.flattenings = [], []
        self.candidate_membership_checks = 0
        cursor, half_size = 0, alphabet * alphabet
        while cursor < len(self.states):
            basis = self.states[cursor]
            matrix = adaptive.flatten_lattice_action(
                basis, local, alphabet, degree, cut, side, infinity)
            factors = adaptive.rational_schur_factors(matrix)
            grouped = defaultdict(list)
            for order, left, right in factors:
                grouped[order].append((left, right))
            blocks = []
            for order, pairs in sorted(grouped.items()):
                rank = len(pairs)
                if rank > maximum_block_rank:
                    raise adaptive.InconclusiveSearch("equal-order block rank budget reached")
                lefts, rights = tuple(zip(*pairs))
                candidates = []
                for mask in range(1, 1 << rank):
                    left = combine(mask, lefts)
                    assert min(arithmetic.rorder(value) for value in left if value != ZERO) == 0
                    rows = tuple(tuple(left[(i * half_size + physical) * dimension + j]
                                       for j in range(dimension))
                                 for i in range(len(basis)) for physical in range(half_size))
                    target = capped_lattice(rows, dimension, precision)
                    if target not in indices:
                        if len(self.states) >= maximum_states:
                            raise adaptive.InconclusiveSearch("matroid state budget reached before closure")
                        indices[target] = len(self.states)
                        self.states.append(target)
                    for row in rows:
                        lattices.membership(row, target)
                    self.candidate_membership_checks += len(rows)
                    candidates.append((mask, indices[target]))
                blocks.append({"order": order, "lefts": lefts, "rights": rights,
                               "candidates": tuple(candidates)})
            self.blocks.append(tuple(blocks))
            self.flattenings.append(matrix)
            cursor += 1
        self.terminal_orders = []
        for basis in self.states:
            values = tuple(linear.scalar_product(row, self.terminal) for row in basis)
            orders = [arithmetic.rorder(value) for value in values if value != ZERO]
            self.terminal_orders.append(min(orders) if orders else None)
        self.closed = True

    def select(self, potential, point):
        """Return and verify the exact factorization selected by a potential."""
        assert len(potential) == len(self.states)
        assert all(isinstance(value, (int, Fraction)) and value > 0 for value in potential)
        assert isinstance(point, Fraction) and point > 1
        transfer, policy, certificates = [], [], []
        entry_checks = 0
        half_size, dimension = self.alphabet * self.alphabet, self.dimension
        for state_index, blocks in enumerate(self.blocks):
            row, choices, factors = {}, [], []
            for block in blocks:
                lefts, rights = block["lefts"], block["rights"]
                selected = greedy_basis(block["candidates"], len(lefts), potential)
                choices.append(selected)
                inverse = binary_inverse(tuple(mask for mask, _ in selected))
                for channel, (mask, target) in enumerate(selected):
                    left = combine(mask, lefts)
                    # L' = L S^T, hence R' = (S^-1)^T R.
                    right_mask = sum(((inverse[i] >> channel) & 1) << i
                                     for i in range(len(lefts)))
                    right = combine(right_mask, rights)
                    rows = tuple(tuple(left[(i * half_size + physical) * dimension + j]
                                       for j in range(dimension))
                                 for i in range(len(self.states[state_index]))
                                 for physical in range(half_size))
                    coordinates = tuple(lattices.membership(values, self.states[target])
                                        for values in rows)
                    order = block["order"]
                    power = (self.degree if self.infinity else 0) - order
                    row[target] = row.get(target, Fraction(0)) + point ** power
                    factors.append({"order": order, "target": target, "mask": mask,
                                    "coordinates": coordinates, "right": right,
                                    "left": left})
            original = self.flattenings[state_index]
            for i, j in product(range(len(original)), range(half_size)):
                reconstructed = ZERO
                for factor in factors:
                    term = arithmetic.rshift(arithmetic.rmul(
                        factor["left"][i], factor["right"][j]), factor["order"])
                    reconstructed = arithmetic.radd(reconstructed, term)
                assert reconstructed == original[i][j]
                entry_checks += 1
            transfer.append(row)
            policy.append(tuple(choices))
            certificates.append(tuple(factors))
        terminal = tuple(point ** (-order) if order is not None else Fraction(0)
                         for order in self.terminal_orders)
        return {"transfer": transfer, "terminal": terminal, "policy": policy,
                "factors": certificates, "entry_reconstruction_checks": entry_checks}

    def summary(self):
        return {"closed": self.closed, "dimension": self.dimension,
                "alphabet": self.alphabet, "degree": self.degree,
                "cut": self.cut, "side": self.side, "infinity": self.infinity,
                "precision": self.precision, "states": self.states, "start": self.start,
                "terminal_orders": self.terminal_orders,
                "initial_coordinates": self.initial_coordinates,
                "candidate_membership_checks": self.candidate_membership_checks,
                "blocks": [[{"order": block["order"], "rank": len(block["lefts"]),
                             "candidates": block["candidates"]} for block in blocks]
                           for blocks in self.blocks]}
