#!/usr/bin/env python3
"""Complete finite-precision channel-basis choices over F_2[[tau]].

At every state this enumerates the valuation, full-tail target, and
residue vector of every primitive regular combination of Smith columns.
The enumeration is finite and supplies polynomial representatives.
Unimodular basis selection is a vector-matroid optimization on residues.
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


def weighted_combination(coefficients, factors, order):
    result = [ZERO] * len(factors[0][1])
    for coefficient, (valuation, left, _) in zip(coefficients, factors):
        if coefficient:
            scalar = arithmetic.rshift(arithmetic.rat(coefficient), valuation - order)
            result = [arithmetic.radd(value, arithmetic.rmul(scalar, entry))
                      for value, entry in zip(result, left)]
    assert all(lattices.regular(value) for value in result)
    assert min(arithmetic.rorder(value) for value in result if value != ZERO) == 0
    return tuple(result)


class DVRSmithLattices:
    """Enumerate all regular basis actions on a declared full-tail state set."""

    def __init__(self, local, degree, *, seed_states, precision=2, alphabet=2,
                 cut=0, side=0, infinity=False, initial=None, terminal=None,
                 maximum_states=1000, maximum_jet_vectors=1000000):
        self.dimension = len(next(iter(local.values())))
        dimension, half_size = self.dimension, alphabet * alphabet
        self.degree, self.alphabet, self.precision = degree, alphabet, precision
        self.cut, self.side, self.infinity = cut, side, infinity
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
        indices = {state: index for index, state in enumerate(self.states)}
        self.start = indices[first]
        self.initial_coordinates = lattices.membership(self.initial, first)
        self.local_factors, self.flattenings = [], []
        self.candidates, self.constant_candidates = [], []
        self.jet_vectors, self.jet_bounds = 0, []
        cursor = 0
        while cursor < len(self.states):
            basis = self.states[cursor]
            matrix = adaptive.flatten_lattice_action(
                basis, local, alphabet, degree, cut, side, infinity)
            factors = adaptive.rational_schur_factors(matrix)
            assert factors, "the declared state has zero local action"
            rank = len(factors)
            orders = tuple(factor[0] for factor in factors)
            jets = tuple(tuple(arithmetic.jet(value, precision) for value in factor[1])
                         for factor in factors)
            residue_columns = [sum((value & 1) << i for i, value in enumerate(vector))
                               for vector in jets]
            assert len(arithmetic.binary_basis(residue_columns)) == rank
            self.jet_bounds.append(precision + max(orders) - min(orders))
            choices, constant_choices, target_cache = {}, {}, {}

            def target_of_jet(vector):
                if vector in target_cache:
                    return target_cache[vector]
                rows = tuple(tuple(arithmetic.rat(vector[(i * half_size + physical) * dimension + j])
                                   for j in range(dimension))
                             for i in range(len(basis)) for physical in range(half_size))
                target = capped_lattice(rows, dimension, precision)
                if target not in indices:
                    if len(self.states) >= maximum_states:
                        raise adaptive.InconclusiveSearch("DVR basis state budget reached before closure")
                    indices[target] = len(self.states)
                    self.states.append(target)
                target_cache[vector] = indices[target]
                return indices[target]

            # Write w_i=tau^(nu_i-a) v_i.  The w_i are regular and at
            # least one has nonzero constant coefficient.  Enumerating
            # their P-jets avoids enumerating irrelevant high v_i terms.
            for order in range(min(orders), max(orders) + 1):
                ranges = [range(0, 1 << precision,
                                1 << max(0, min(precision, valuation - order)))
                          for valuation in orders]
                for weights in product(*ranges):
                    self.jet_vectors += 1
                    if self.jet_vectors > maximum_jet_vectors:
                        raise adaptive.InconclusiveSearch("DVR basis jet budget reached before closure")
                    if not any(weight & 1 for weight in weights):
                        continue
                    residue, free, representative = 0, [], []
                    for i, (valuation, weight) in enumerate(zip(orders, weights)):
                        delta = valuation - order
                        if delta <= 0:
                            representative.append(weight << (-delta))
                            if delta == 0 and weight & 1:
                                residue |= 1 << i
                        elif delta < precision:
                            representative.append(weight >> delta)
                            if weight >> delta & 1:
                                residue |= 1 << i
                        else:
                            representative.append(0)
                            free.append(i)
                    if not residue and not free:
                        continue
                    vector = [0] * len(jets[0])
                    for weight, left in zip(weights, jets):
                        if weight:
                            vector = [entry ^ (arithmetic.pmul(weight, value) & ((1 << precision) - 1))
                                      for entry, value in zip(vector, left)]
                    target = target_of_jet(tuple(vector))
                    for mask in range(1 << len(free)):
                        full_residue, coefficients = residue, list(representative)
                        for bit, i in enumerate(free):
                            if mask >> bit & 1:
                                full_residue |= 1 << i
                                coefficients[i] = 1
                        if full_residue:
                            choices.setdefault((full_residue, target, order), tuple(coefficients))

            # Independently construct the complete constant-basis subfamily.
            for residue in range(1, 1 << rank):
                coefficients = tuple((residue >> i) & 1 for i in range(rank))
                order = min(valuation for valuation, value in zip(orders, coefficients) if value)
                left = weighted_combination(coefficients, factors, order)
                target = target_of_jet(tuple(arithmetic.jet(value, precision) for value in left))
                key = residue, target, order
                assert key in choices
                constant_choices[key] = coefficients
            self.candidates.append(choices)
            self.constant_candidates.append(constant_choices)
            self.local_factors.append(factors)
            self.flattenings.append(matrix)
            cursor += 1
        self.terminal_orders = []
        for state in self.states:
            values = [linear.scalar_product(row, self.terminal) for row in state]
            orders = [arithmetic.rorder(value) for value in values if value != ZERO]
            self.terminal_orders.append(min(orders) if orders else None)
        self.closed = True

    def choose(self, potential, point, *, constant_only=False):
        assert len(potential) == len(self.states)
        assert all(isinstance(value, (int, Fraction)) and value > 0 for value in potential)
        assert isinstance(point, Fraction) and point > 1
        families = self.constant_candidates if constant_only else self.candidates
        choices, transfer = [], []
        for index, family in enumerate(families):
            rank = len(self.local_factors[index])
            best_by_residue = {}
            for key, coefficients in family.items():
                residue, target, order = key
                weight = point ** ((self.degree if self.infinity else 0) - order) * potential[target]
                candidate = weight, key, coefficients
                if residue not in best_by_residue or candidate < best_by_residue[residue]:
                    best_by_residue[residue] = candidate
            basis, selected, row = [], [], {}
            for _, key, coefficients in sorted(best_by_residue.values()):
                residue, target, order = key
                enlarged = arithmetic.binary_basis(basis + [residue])
                if len(enlarged) > len(basis):
                    basis = list(enlarged)
                    selected.append((key, coefficients))
                    power = (self.degree if self.infinity else 0) - order
                    row[target] = row.get(target, Fraction(0)) + point ** power
                    if len(basis) == rank:
                        break
            assert len(basis) == rank
            choices.append(tuple(selected))
            transfer.append(row)
        return transfer, choices

    def verify(self, choices):
        """Reconstruct actual rational factors for a selected unimodular policy."""
        certificates, entry_checks, membership_checks = [], 0, 0
        dimension, half_size = self.dimension, self.alphabet * self.alphabet
        for state_index, selected in enumerate(choices):
            original_factors = self.local_factors[state_index]
            rank = len(original_factors)
            change = tuple(tuple(arithmetic.rat(value) for value in coefficients)
                           for _, coefficients in selected)
            inverse = linear.inverse(change)
            assert all(lattices.regular(value) for row in inverse for value in row)
            factors = []
            for channel, ((residue, target, order), coefficients) in enumerate(selected):
                assert sum((value & 1) << i for i, value in enumerate(coefficients)) == residue
                left = weighted_combination(coefficients, original_factors, order)
                right = tuple(ZERO for _ in range(half_size))
                for i in range(rank):
                    right = tuple(arithmetic.radd(old, arithmetic.rmul(inverse[i][channel], value))
                                  for old, value in zip(right, original_factors[i][2]))
                assert all(lattices.regular(value) for value in right)
                rows = tuple(tuple(left[(i * half_size + physical) * dimension + j]
                                   for j in range(dimension))
                             for i in range(len(self.states[state_index]))
                             for physical in range(half_size))
                assert capped_lattice(rows, dimension, self.precision) == self.states[target]
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
                "terminal_orders": self.terminal_orders, "jet_vectors": self.jet_vectors,
                "complete_jet_bounds": self.jet_bounds,
                "candidate_counts": list(map(len, self.candidates)),
                "constant_candidate_counts": list(map(len, self.constant_candidates))}
