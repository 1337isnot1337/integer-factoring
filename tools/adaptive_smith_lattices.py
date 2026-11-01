#!/usr/bin/env python3
"""State-dependent Smith factorizations over F_2(tau).

The input is an ordered polynomial crossing.  At each exact auxiliary row
lattice L, this code factors the complete physical flattening of B_L K.
Unlike a fixed-channel module bound, it may change the channel basis at
every state.  Closure and every local factorization are checked exactly.
Exhausting a state or degree budget is an inconclusive search outcome.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import check_cancellation_module_automata as arithmetic
import check_exact_module_orbits as lattices
import check_minimal_boundary_modules as linear


ZERO, ONE = arithmetic.ZERO, arithmetic.ONE


class InconclusiveSearch(RuntimeError):
    """A finite resource limit was reached without a closed certificate."""


def rational_schur_factors(matrix):
    """Exact regular factors of any rank, with a checked reconstruction.

    Each output is (order, left, right), representing tau**order left right.
    Both vectors are regular, and each contains a unit.  The entry of least
    valuation is used at every step; identically zero residuals stop exactly.
    """
    original = tuple(tuple(row) for row in matrix)
    width = len(original[0]) if original else 0
    assert all(len(row) == width for row in original)
    assert all(lattices.regular(value) for row in original for value in row)
    work = [list(row) for row in original]
    factors = []
    while True:
        choices = [(arithmetic.rorder(value), i, j)
                   for i, row in enumerate(work)
                   for j, value in enumerate(row) if value != ZERO]
        if not choices:
            break
        order, row, column = min(choices)
        pivot = work[row][column]
        left = tuple(arithmetic.rshift(values[column], -order) for values in work)
        right = tuple(arithmetic.rdiv(value, pivot) for value in work[row])
        assert all(lattices.regular(value) for value in left + right)
        assert min(arithmetic.rorder(value) for value in left if value != ZERO) == 0
        assert min(arithmetic.rorder(value) for value in right if value != ZERO) == 0
        factors.append((order, left, right))
        for i, j in product(range(len(work)), range(width)):
            correction = arithmetic.rshift(arithmetic.rmul(left[i], right[j]), order)
            work[i][j] = arithmetic.radd(work[i][j], correction)
    for i, j in product(range(len(original)), range(width)):
        reconstructed = ZERO
        for order, left, right in factors:
            reconstructed = arithmetic.radd(
                reconstructed,
                arithmetic.rshift(arithmetic.rmul(left[i], right[j]), order),
            )
        assert reconstructed == original[i][j]
    return tuple(factors)


def flatten_lattice_action(basis, local, alphabet, degree, cut, side, infinity):
    """Rows (lattice coordinate, physical half, outgoing auxiliary index)."""
    dimension = len(next(iter(local.values())))
    half_size = alphabet * alphabet
    matrix = [[ZERO] * half_size
              for _ in range(len(basis) * half_size * dimension)]
    for (x, y, u, v), value in local.items():
        left, right = ((alphabet * x + y, alphabet * u + v) if cut == 0
                       else (alphabet * x + u, alphabet * y + v))
        physical, other = (left, right) if side == 0 else (right, left)
        if infinity:
            transformed = tuple(tuple(arithmetic.rat(sum(
                1 << (degree - k)
                for k in range(entry.bit_length()) if entry >> k & 1
            )) for entry in row) for row in value)
        else:
            transformed = tuple(tuple(arithmetic.rat(entry) for entry in row)
                                for row in value)
        result = linear.matrix_product(basis, transformed)
        for i, row in enumerate(result):
            for j, entry in enumerate(row):
                matrix[(i * half_size + physical) * dimension + j][other] = entry
    return tuple(tuple(row) for row in matrix)


class AdaptiveSmithLattices:
    """A closed state-dependent factorization graph, exact or full-tail capped.

    The default states are exact DVR lattices and may be rank deficient.
    With precision P, every image is enlarged by tau**P O**D, giving a
    finite-state sufficient bound.  No tail is silently added in exact mode.
    """

    def __init__(self, local, degree, *, alphabet=2, cut=0, side=0,
                 infinity=False, initial=None, terminal=None, precision=None,
                 maximum_states=1000, maximum_degree=1024):
        self.dimension = len(next(iter(local.values())))
        self.alphabet, self.degree = alphabet, degree
        self.cut, self.side, self.infinity = cut, side, infinity
        self.precision = precision
        dimension, half_size = self.dimension, alphabet * alphabet
        assert alphabet >= 2 and degree >= 0
        assert cut in (0, 1) and side in (0, 1)
        assert precision is None or precision >= 1
        assert set(local) == set(product(range(alphabet), repeat=4))
        assert all(len(matrix) == dimension and all(len(row) == dimension for row in matrix)
                   for matrix in local.values())
        assert all(isinstance(value, int) and 0 <= value < (1 << (degree + 1))
                   for matrix in local.values() for row in matrix for value in row)
        self.initial = tuple(initial or (ONE,) + (ZERO,) * (dimension - 1))
        self.terminal = tuple(terminal or (ONE,) + (ZERO,) * (dimension - 1))
        assert len(self.initial) == len(self.terminal) == dimension
        assert all(value in (ZERO, ONE) for value in self.initial + self.terminal), \
            "this polynomial-crossing interface uses constant binary boundaries"
        assert any(value != ZERO for value in self.initial)
        assert min(arithmetic.rorder(value) for value in self.initial if value != ZERO) == 0

        def state_basis(generators):
            if precision is None:
                return lattices.hermite_lattice(generators, dimension)
            residues = arithmetic.module_basis(
                [tuple(arithmetic.jet(value, precision) for value in row)
                 for row in generators], dimension, precision)
            lifted = arithmetic.lift_generators(residues, dimension, precision)
            return lattices.hermite_lattice(
                tuple(tuple(arithmetic.rat(value) for value in row) for row in lifted),
                dimension)

        first = state_basis((self.initial,))
        self.initial_coordinates = lattices.membership(self.initial, first)
        self.states = [first]
        indices = {first: 0}
        self.transitions = []
        self.factors = []
        self.local_entry_checks = 0
        self.regular_coordinate_checks = 0
        cursor = 0
        while cursor < len(self.states):
            basis = self.states[cursor]
            matrix = flatten_lattice_action(
                basis, local, alphabet, degree, cut, side, infinity)
            pieces = rational_schur_factors(matrix)
            self.local_entry_checks += len(matrix) * half_size
            transitions, certificates = [], []
            for order, left, right in pieces:
                rows = tuple(tuple(left[(i * half_size + physical) * dimension + j]
                                   for j in range(dimension))
                             for i in range(len(basis)) for physical in range(half_size))
                target = state_basis(rows)
                assert target
                assert min(arithmetic.rorder(value) for row in target for value in row
                           if value != ZERO) == 0
                if max(max(value[0].bit_length(), value[1].bit_length())
                       for row in target for value in row) > maximum_degree:
                    raise InconclusiveSearch("rational coefficient degree budget reached")
                coordinates = tuple(lattices.membership(row, target) for row in rows)
                assert all(lattices.regular(value) for row in coordinates for value in row)
                self.regular_coordinate_checks += len(rows)
                if target not in indices:
                    if len(self.states) >= maximum_states:
                        raise InconclusiveSearch("state budget reached before closure")
                    indices[target] = len(self.states)
                    self.states.append(target)
                target_index = indices[target]
                transitions.append((order, target_index))
                certificates.append({
                    "order": order, "target": target_index,
                    "coordinates": coordinates, "right": right,
                })
            self.transitions.append(tuple(transitions))
            self.factors.append(tuple(certificates))
            cursor += 1
        self.terminal_orders = []
        for basis in self.states:
            values = tuple(linear.scalar_product(row, self.terminal) for row in basis)
            orders = [arithmetic.rorder(value) for value in values if value != ZERO]
            self.terminal_orders.append(min(orders) if orders else None)
        self.closed = True

    def weighted_matrix(self, point):
        assert isinstance(point, Fraction) and point > 1
        matrix = []
        for transitions in self.transitions:
            row = {}
            for order, target in transitions:
                power = (self.degree if self.infinity else 0) - order
                row[target] = row.get(target, Fraction(0)) + point ** power
            matrix.append(row)
        terminal = tuple(point ** (-order) if order is not None else Fraction(0)
                         for order in self.terminal_orders)
        return matrix, terminal

    def summary(self, include_factors=False):
        result = {
            "closed": self.closed, "dimension": self.dimension,
            "alphabet": self.alphabet, "degree": self.degree,
            "cut": self.cut, "side": self.side, "infinity": self.infinity,
            "precision": self.precision,
            "states": self.states, "transitions": self.transitions,
            "terminal_orders": self.terminal_orders,
            "initial_coordinates": self.initial_coordinates,
            "local_entry_checks": self.local_entry_checks,
            "regular_coordinate_checks": self.regular_coordinate_checks,
        }
        if include_factors:
            result["factors"] = self.factors
        return result
