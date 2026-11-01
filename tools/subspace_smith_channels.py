#!/usr/bin/env python3
"""Polynomial-size residue-subspace oracle for regular channel bases.

For each charged order and declared full-tail target lattice, solve a
binary linear system and retain polynomial lifts of its residue space.
Minimum-weight basis selection then needs only these subspace bases,
not an enumeration of primitive channel jets or residue vectors.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import adaptive_smith_lattices as adaptive
import check_cancellation_module_automata as arithmetic
import check_exact_module_orbits as lattices
import check_minimal_boundary_modules as linear
from matroid_smith_channels import capped_lattice


ZERO, ONE = arithmetic.ZERO, arithmetic.ONE


def binary_remainder(value, basis):
    for row in basis:
        if value >> (row.bit_length() - 1) & 1:
            value ^= row
    return value


def kernel_relations(columns):
    """A basis of the binary kernel, retaining original input coordinates."""
    image_basis, relations = {}, []
    for index, column in enumerate(columns):
        value, coefficients = column, 1 << index
        for pivot in sorted(image_basis, reverse=True):
            if value >> pivot & 1:
                old_value, old_coefficients = image_basis[pivot]
                value ^= old_value
                coefficients ^= old_coefficients
        if value:
            image_basis[value.bit_length() - 1] = value, coefficients
        else:
            relations.append(coefficients)
    assert len(relations) + len(image_basis) == len(columns)
    return tuple(relations)


def projected_basis(relations, residue_columns):
    """Canonical residue-space basis with preimages in variable coordinates."""
    pivots = {}
    for coefficients in relations:
        residue, remaining = 0, coefficients
        while remaining:
            bit = remaining & -remaining
            residue ^= residue_columns[bit.bit_length() - 1]
            remaining ^= bit
        for pivot in sorted(pivots, reverse=True):
            if residue >> pivot & 1:
                old_residue, old_coefficients = pivots[pivot]
                residue ^= old_residue
                coefficients ^= old_coefficients
        if not residue:
            continue
        pivot = residue.bit_length() - 1
        for old_pivot in tuple(pivots):
            old_residue, old_coefficients = pivots[old_pivot]
            if old_residue >> pivot & 1:
                pivots[old_pivot] = old_residue ^ residue, old_coefficients ^ coefficients
        pivots[pivot] = residue, coefficients
    return tuple(pivots[pivot] for pivot in sorted(pivots, reverse=True))


def regular_combination(coefficients, factors, order):
    """Form tau^-order Uv; extra common valuation is permitted."""
    result = [ZERO] * len(factors[0][1])
    for coefficient, (valuation, left, _) in zip(coefficients, factors):
        if coefficient:
            scalar = arithmetic.rshift(arithmetic.rat(coefficient), valuation - order)
            result = [arithmetic.radd(old, arithmetic.rmul(scalar, value))
                      for old, value in zip(result, left)]
    assert all(lattices.regular(value) for value in result)
    assert any(value != ZERO for value in result)
    return tuple(result)


def residue_subspace(factors, target, dimension, precision, order):
    """Compute residues of all v with Uv in tau^order times target.

    Flattened left vectors are grouped into consecutive auxiliary rows
    of length dimension. The normalized Schur columns must be independent
    modulo tau. Return a residue basis and polynomial coefficient lifts.
    """
    rank = len(factors)
    length = len(factors[0][1])
    assert length % dimension == 0 and precision >= 1 and order >= 0
    assert all(len(left) == length for _, left, _ in factors)
    row_count = length // dimension
    mask = (1 << precision) - 1
    target_basis = arithmetic.module_basis(
        [tuple(arithmetic.jet(value, precision) for value in row) for row in target],
        dimension, precision,
    )
    jets = [tuple(arithmetic.jet(value, precision) for value in left)
            for _, left, _ in factors]
    variable_terms, columns, residue_columns = [], [], []
    for index, ((valuation, _, _), jet) in enumerate(zip(factors, jets)):
        delta = valuation - order
        if delta >= precision:
            # Its constant coefficient is invisible modulo the target tail.
            variable_terms.append((index, 0))
            columns.append(0)
            residue_columns.append(1 << index)
            continue
        for degree in range(max(0, delta), precision):
            obstruction = 0
            for row in range(row_count):
                packed = arithmetic.pack(
                    tuple((jet[row * dimension + coordinate] << degree) & mask
                          for coordinate in range(dimension)), precision,
                )
                remainder = binary_remainder(packed, target_basis)
                obstruction |= remainder << (row * dimension * precision)
            exponent = degree - delta
            variable_terms.append((index, exponent))
            columns.append(obstruction)
            residue_columns.append((1 << index) if exponent == 0 else 0)
    relations = kernel_relations(columns)
    projected = projected_basis(relations, residue_columns)
    answer = []
    for residue, variables in projected:
        coefficients = [0] * rank
        remaining = variables
        while remaining:
            bit = remaining & -remaining
            index, exponent = variable_terms[bit.bit_length() - 1]
            coefficients[index] ^= 1 << exponent
            remaining ^= bit
        assert sum((value & 1) << index for index, value in enumerate(coefficients)) == residue
        answer.append((residue, tuple(coefficients)))
    return {
        "basis": tuple(answer),
        "variables": len(columns),
        "constraint_bits": row_count * dimension * precision,
        "kernel_dimension": len(relations),
        "residue_dimension": len(answer),
    }


class SubspaceSmithLattices:
    """Regular channel policies using a fixed declared family of targets.

    Membership in a target is sufficient; the exact generated lattice
    may be smaller. A charged order may be below the actual order.
    The target family is fixed, so every selected transition is closed.
    """

    def __init__(self, local, degree, *, seed_states, precision=2, alphabet=2,
                 cut=0, side=0, infinity=False, initial=None, terminal=None,
                 include_lower_orders=True):
        self.dimension = len(next(iter(local.values())))
        dimension, half_size = self.dimension, alphabet * alphabet
        self.degree, self.alphabet, self.precision = degree, alphabet, precision
        self.cut, self.side, self.infinity = cut, side, infinity
        self.include_lower_orders = include_lower_orders
        assert precision >= 1 and degree >= 0 and alphabet >= 2
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
        self.states = list(seed_states)
        assert len(self.states) == len(set(self.states))
        assert all(state == capped_lattice(state, dimension, precision) for state in self.states)
        assert all(min(arithmetic.rorder(value) for row in state for value in row
                       if value != ZERO) == 0 for state in self.states)
        first = capped_lattice((self.initial,), dimension, precision)
        assert first in self.states
        self.start = self.states.index(first)
        self.initial_coordinates = lattices.membership(self.initial, first)
        self.local_factors, self.flattenings, self.actions = [], [], []
        self.linear_systems = self.total_variables = self.total_residue_dimensions = 0
        self.maximum_variables = 0
        for basis in self.states:
            matrix = adaptive.flatten_lattice_action(
                basis, local, alphabet, degree, cut, side, infinity)
            factors = adaptive.rational_schur_factors(matrix)
            assert factors, "the declared state has zero local action"
            rank = len(factors)
            orders = tuple(factor[0] for factor in factors)
            residues = [sum((arithmetic.jet(value, 1) & 1) << i
                            for i, value in enumerate(factor[1])) for factor in factors]
            assert len(arithmetic.binary_basis(residues)) == rank
            actions = []
            first_order = 0 if include_lower_orders else min(orders)
            for order in range(first_order, max(orders) + 1):
                for target_index, target in enumerate(self.states):
                    subspace = residue_subspace(factors, target, dimension, precision, order)
                    self.linear_systems += 1
                    self.total_variables += subspace["variables"]
                    self.maximum_variables = max(self.maximum_variables, subspace["variables"])
                    self.total_residue_dimensions += subspace["residue_dimension"]
                    if subspace["basis"]:
                        actions.append({"target": target_index, "order": order, **subspace})
            assert len(arithmetic.binary_basis(
                [residue for action in actions for residue, _ in action["basis"]])) == rank
            self.actions.append(tuple(actions))
            self.local_factors.append(factors)
            self.flattenings.append(matrix)
        self.terminal_orders = []
        for state in self.states:
            values = [linear.scalar_product(row, self.terminal) for row in state]
            orders = [arithmetic.rorder(value) for value in values if value != ZERO]
            self.terminal_orders.append(min(orders) if orders else None)
        self.closed = True

    def choose(self, potential, point):
        assert len(potential) == len(self.states)
        assert all(isinstance(value, (int, Fraction)) and value > 0 for value in potential)
        assert isinstance(point, Fraction) and point > 1
        choices, transfer = [], []
        for index, actions in enumerate(self.actions):
            rank = len(self.local_factors[index])
            weighted = sorted(actions, key=lambda action: (
                point ** ((self.degree if self.infinity else 0) - action["order"])
                * potential[action["target"]], action["order"], action["target"],
            ))
            basis, selected, row = [], [], {}
            for action in weighted:
                target, order = action["target"], action["order"]
                for residue, coefficients in action["basis"]:
                    enlarged = arithmetic.binary_basis(basis + [residue])
                    if len(enlarged) > len(basis):
                        basis = list(enlarged)
                        selected.append(((residue, target, order), coefficients))
                        exponent = (self.degree if self.infinity else 0) - order
                        row[target] = row.get(target, Fraction(0)) + point ** exponent
                        if len(basis) == rank:
                            break
                if len(basis) == rank:
                    break
            assert len(basis) == rank
            choices.append(tuple(selected))
            transfer.append(row)
        return transfer, choices

    def verify(self, choices):
        certificates, entry_checks, membership_checks = [], 0, 0
        dimension, half_size = self.dimension, self.alphabet * self.alphabet
        for state_index, selected in enumerate(choices):
            original = self.local_factors[state_index]
            rank = len(original)
            change = tuple(tuple(arithmetic.rat(value) for value in coefficients)
                           for _, coefficients in selected)
            inverse = linear.inverse(change)
            assert all(lattices.regular(value) for row in inverse for value in row)
            factors = []
            for channel, ((residue, target, order), coefficients) in enumerate(selected):
                assert sum((value & 1) << i for i, value in enumerate(coefficients)) == residue
                left = regular_combination(coefficients, original, order)
                right = (ZERO,) * half_size
                for i in range(rank):
                    right = tuple(arithmetic.radd(old, arithmetic.rmul(inverse[i][channel], value))
                                  for old, value in zip(right, original[i][2]))
                assert all(lattices.regular(value) for value in right)
                rows = tuple(tuple(left[(i * half_size + physical) * dimension + j]
                                   for j in range(dimension))
                             for i in range(len(self.states[state_index]))
                             for physical in range(half_size))
                coordinates = tuple(lattices.membership(row, self.states[target]) for row in rows)
                membership_checks += len(rows)
                factors.append({"order": order, "target": target, "left": left,
                                "right": right, "coordinates": coordinates,
                                "residue": residue, "basis_row_polynomials": coefficients})
            matrix = self.flattenings[state_index]
            for i, j in product(range(len(matrix)), range(half_size)):
                reconstructed = ZERO
                for factor in factors:
                    term = arithmetic.rshift(arithmetic.rmul(factor["left"][i], factor["right"][j]),
                                             factor["order"])
                    reconstructed = arithmetic.radd(reconstructed, term)
                assert reconstructed == matrix[i][j]
                entry_checks += 1
            certificates.append({"change": change, "inverse": inverse, "factors": factors})
        return {"states": certificates, "local_entry_checks": entry_checks,
                "regular_membership_checks": membership_checks}

    def terminal_weights(self, point):
        return tuple(point ** (-order) if order is not None else Fraction(0)
                     for order in self.terminal_orders)

    def summary(self):
        return {"closed": self.closed, "dimension": self.dimension,
                "alphabet": self.alphabet, "degree": self.degree,
                "cut": self.cut, "side": self.side, "infinity": self.infinity,
                "precision": self.precision, "states": self.states, "start": self.start,
                "terminal_orders": self.terminal_orders,
                "linear_systems": self.linear_systems,
                "total_variables": self.total_variables,
                "maximum_variables": self.maximum_variables,
                "total_residue_dimensions": self.total_residue_dimensions,
                "action_counts": list(map(len, self.actions)),
                "contains_enlarged_targets": True,
                "contains_undercharged_orders": True,
                "includes_orders_below_minimum_smith_order": self.include_lower_orders}
