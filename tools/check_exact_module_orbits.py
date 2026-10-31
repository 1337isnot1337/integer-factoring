#!/usr/bin/env python3
"""Exact row-lattice orbits, including rank-deficient states, over F_2[[tau]].

All matrix entries are rational functions in tau.  A completed orbit is an
exact certificate; exhausting a resource budget is an inconclusive result.
The checker has no numerical or third-party dependencies.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

import check_cancellation_module_automata as module
import check_minimal_boundary_modules as minimal


ZERO, ONE = module.ZERO, module.ONE


def regular(value):
    return value == ZERO or module.rorder(value) >= 0


def combination(coefficients, vectors, width):
    answer = [ZERO] * width
    for coefficient, vector in zip(coefficients, vectors):
        for j, value in enumerate(vector):
            answer[j] = module.radd(answer[j], module.rmul(coefficient, value))
    return tuple(answer)


def membership(vector, basis):
    """Return regular coefficients expressing vector in a Hermite basis."""
    remainder = list(vector)
    coefficients = []
    for row in basis:
        pivot = next(j for j, value in enumerate(row) if value != ZERO)
        coefficient = module.rdiv(remainder[pivot], row[pivot])
        assert regular(coefficient), "vector is not in this row lattice"
        coefficients.append(coefficient)
        remainder = [module.radd(value, module.rmul(coefficient, old))
                     for value, old in zip(remainder, row)]
    assert all(value == ZERO for value in remainder), "vector is outside the row span"
    assert combination(coefficients, basis, len(vector)) == tuple(vector)
    return tuple(coefficients)


def hermite_lattice(generators, dimension):
    """Canonical DVR row Hermite basis, with exact two-sided span witnesses.

    The leftmost nonzero entry of each row is a pure power of tau.  Entries
    above that pivot are polynomials of smaller degree.  Missing pivots are
    retained as missing pivots; no finite-precision tail is inserted.
    """
    original = tuple(tuple(row) for row in generators)
    assert all(len(row) == dimension for row in original)
    assert all(regular(value) for row in original for value in row)
    rows, transforms = [], []
    for i, row in enumerate(original):
        if any(value != ZERO for value in row):
            rows.append(list(row))
            transforms.append([ONE if i == j else ZERO for j in range(len(original))])
    used = 0
    for column in range(dimension):
        choices = [(module.rorder(rows[i][column]), i)
                   for i in range(used, len(rows)) if rows[i][column] != ZERO]
        if not choices:
            continue
        order, index = min(choices)
        rows[used], rows[index] = rows[index], rows[used]
        transforms[used], transforms[index] = transforms[index], transforms[used]
        unit = module.rshift(rows[used][column], -order)
        assert module.rorder(unit) == 0
        rows[used] = [module.rdiv(value, unit) for value in rows[used]]
        transforms[used] = [module.rdiv(value, unit) for value in transforms[used]]
        assert rows[used][column] == (1 << order, 1)
        for i in range(used + 1, len(rows)):
            factor = module.rshift(rows[i][column], -order)
            assert regular(factor)
            if factor != ZERO:
                rows[i] = [module.radd(value, module.rmul(factor, pivot))
                           for value, pivot in zip(rows[i], rows[used])]
                transforms[i] = [module.radd(value, module.rmul(factor, pivot))
                                 for value, pivot in zip(transforms[i], transforms[used])]
        for i in range(used):
            remainder = (module.jet(rows[i][column], order), 1) if order else ZERO
            factor = module.rshift(module.radd(rows[i][column], remainder), -order)
            assert regular(factor)
            rows[i] = [module.radd(value, module.rmul(factor, pivot))
                       for value, pivot in zip(rows[i], rows[used])]
            transforms[i] = [module.radd(value, module.rmul(factor, pivot))
                             for value, pivot in zip(transforms[i], transforms[used])]
            assert rows[i][column] == remainder
        used += 1
        if used == len(rows):
            break
    assert all(all(value == ZERO for value in row) for row in rows[used:])
    basis = tuple(tuple(row) for row in rows[:used])
    for row, coefficients in zip(basis, transforms[:used]):
        assert all(regular(value) for value in coefficients)
        assert combination(coefficients, original, dimension) == row
    for row in original:
        membership(row, basis)
    return basis


class ExactModuleOrbit:
    """Enumerate normalized exact row lattices from the free initial lattice."""

    def __init__(self, channels, maximum_states=1000, maximum_degree=1024):
        self.channels = tuple(channels)
        self.dimension = len(channels[0][0])
        dimension = self.dimension
        assert all(len(matrix) == dimension and all(len(row) == dimension for row in matrix)
                   for choices in channels for matrix in choices)
        assert all(regular(value) for choices in channels for matrix in choices
                   for row in matrix for value in row)
        initial = tuple(tuple(ONE if i == j else ZERO for j in range(dimension))
                        for i in range(dimension))
        self.states = [initial]
        state_index = {initial: 0}
        self.transitions = {}
        position = 0
        while position < len(self.states):
            state = self.states[position]
            for symbol, choices in enumerate(channels):
                generators = [module.rational_vector_product(row, matrix)
                              for row, matrix in product(state, choices)]
                image = hermite_lattice(generators, dimension)
                if not image:
                    self.transitions[position, symbol] = None
                    continue
                order = min(module.rorder(value) for row in image for value in row
                            if value != ZERO)
                target = tuple(tuple(module.rshift(value, -order) for value in row)
                               for row in image)
                assert target == hermite_lattice(target, dimension)
                assert min(module.rorder(value) for row in target for value in row
                           if value != ZERO) == 0
                assert max(max(value[0].bit_length(), value[1].bit_length())
                           for row in target for value in row) <= maximum_degree, "degree budget"
                if target not in state_index:
                    assert len(self.states) < maximum_states, "state budget"
                    state_index[target] = len(self.states)
                    self.states.append(target)
                self.transitions[position, symbol] = order, state_index[target]
            position += 1
        self.closed = True

    def word_order(self, word):
        state, order = 0, 0
        for symbol in word:
            entry = self.transitions[state, symbol]
            if entry is None:
                return None
            shift, state = entry
            order += shift
        return order

    def weighted_matrix(self, exponents, point):
        matrix = [dict() for _ in self.states]
        for (row, symbol), entry in self.transitions.items():
            if entry is not None:
                order, column = entry
                matrix[row][column] = matrix[row].get(column, Fraction(0)) + point ** (
                    exponents[symbol] - order)
        return matrix


def companion_two():
    """Conjugated version of C A(z) + z(E+C) e f, with A^2=z I."""
    local = {}
    change = ((ONE, ZERO), (ONE, ONE))
    for x, y, u, v in product(range(2), repeat=4):
        identity = int(x == u and y == v)
        crossing = int(x == v and y == u)
        nilpotent = identity ^ crossing
        original = (((nilpotent << 1, 1), (identity, 1)),
                    ((crossing << 1, 1), ZERO))
        transformed = minimal.matrix_product(minimal.matrix_product(change, original), change)
        assert all(value[1] == 1 for row in transformed for value in row)
        polynomial = [[value[0] for value in row] for row in transformed]
        assert all(polynomial[i][0].bit_count() % 2 == int(i == 0 and crossing)
                   for i in range(2))
        local[x, y, u, v] = polynomial
    return local


def rational_potential(matrix, rate):
    """Solve (rate I - T)h=1 exactly, then verify h>0 and Th<rate h."""
    size = len(matrix)
    work = [[(rate if i == j else Fraction(0)) - matrix[i].get(j, Fraction(0))
             for j in range(size)] + [Fraction(1)] for i in range(size)]
    for column in range(size):
        pivot = next(i for i in range(column, size) if work[i][column])
        work[column], work[pivot] = work[pivot], work[column]
        divisor = work[column][column]
        work[column] = [value / divisor for value in work[column]]
        for i in range(size):
            if i != column:
                factor = work[i][column]
                work[i] = [value - factor * old for value, old in zip(work[i], work[column])]
    vector = tuple(row[-1] for row in work)
    assert all(value > 0 for value in vector)
    for i, row in enumerate(matrix):
        assert sum(value * vector[j] for j, value in row.items()) == rate * vector[i] - 1
    return vector


def lump_matrix(matrix, classes):
    size = max(classes) + 1
    result = [None] * size
    for i, row in enumerate(matrix):
        aggregate = {}
        for j, value in row.items():
            target = classes[j]
            aggregate[target] = aggregate.get(target, Fraction(0)) + value
        if result[classes[i]] is None:
            result[classes[i]] = aggregate
        else:
            assert result[classes[i]] == aggregate
    assert all(row is not None for row in result)
    return result


def symbolic_quotient(orbit, exponents, classes):
    """Check lumping as a Laurent-polynomial identity, for every positive s."""
    rows = [dict() for _ in orbit.states]
    for (source, symbol), entry in orbit.transitions.items():
        if entry is not None:
            order, target = entry
            key = classes[target], exponents[symbol] - order
            rows[source][key] = rows[source].get(key, 0) + 1
    result = [None] * (max(classes) + 1)
    for source, row in enumerate(rows):
        if result[classes[source]] is None:
            result[classes[source]] = row
        else:
            assert result[classes[source]] == row
    return result


def check_random_lattices():
    rng = Random(0xD7A)
    checked = 0
    for dimension in range(1, 4):
        for _ in range(20):
            generators = [tuple(module.rat(rng.randrange(32), 2 * rng.randrange(16) + 1)
                                for _ in range(dimension)) for _ in range(rng.randrange(1, 7))]
            basis = hermite_lattice(generators, dimension)
            for precision in range(1, 6):
                left = module.module_basis(
                    [tuple(module.jet(value, precision) for value in row) for row in generators],
                    dimension, precision)
                right = module.module_basis(
                    [tuple(module.jet(value, precision) for value in row) for row in basis],
                    dimension, precision)
                assert left == right
                checked += 1
            rng.shuffle(generators)
            assert hermite_lattice(generators, dimension) == basis
    return checked


def check():
    local = companion_two()
    records, word_checks = [], 0
    point = Fraction(5)
    rates = (Fraction(266, 25), Fraction(607, 500))
    for cut, infinity in ((0, True), (1, False)):
        factors = module.local_factors(local, cut, 0, (0, 0), infinity)
        representation = minimal.regularize(
            minimal.minimize(factors["channels"], (ONE, ZERO), (ONE, ZERO)), factors["exponents"])
        assert representation["dimension"] == 2
        orbit = ExactModuleOrbit(representation["channels"])
        assert len(orbit.states) == (2 if infinity else 7)
        assert sorted(set(map(len, orbit.states))) == ([2] if infinity else [1, 2])
        for length in range(4):
            for word in product(range(len(orbit.channels)), repeat=length):
                literal = minimal.word_order(orbit.channels, representation["initial"],
                                             representation["terminal"], word, free=True)
                assert literal == orbit.word_order(word)
                word_checks += 1
        matrix = orbit.weighted_matrix(representation["exponents"], point)
        potential = rational_potential(matrix, rates[cut])
        if infinity:
            assert matrix == [{1: 4 * point}, {0: Fraction(3), 1: point}]
            classes = (0, 1)
            assert symbolic_quotient(orbit, representation["exponents"], classes) == [
                {(1, 1): 4}, {(0, 0): 3, (1, 1): 1}]
        else:
            classes = (0, 1, 2, 2, 3, 1, 2)
            quotient = lump_matrix(matrix, classes)
            assert quotient == [
                {0: 1 / point, 1: Fraction(1), 2: 2 / point},
                {0: 1 / point, 2: 3 / point},
                {1: Fraction(1), 2: 2 / point, 3: 1 / point},
                {1: 1 + 1 / point, 2: 2 / point ** 2},
            ]
            assert symbolic_quotient(orbit, representation["exponents"], classes) == [
                {(0, -1): 1, (1, 0): 1, (2, -1): 2},
                {(0, -1): 1, (2, -1): 3},
                {(1, 0): 1, (2, -1): 2, (3, -1): 1},
                {(1, 0): 1, (1, -1): 1, (2, -2): 2},
            ]
        records.append({"cut": cut, "infinity": infinity, "states": orbit.states,
                        "exponents": representation["exponents"],
                        "transitions": [{"state": state, "symbol": symbol,
                                         "image": entry}
                                        for (state, symbol), entry in orbit.transitions.items()],
                        "classes": classes, "point": str(point), "upper_rate": str(rates[cut]),
                        "potential": [str(value) for value in potential]})
    assert 44 ** 2 - 4 * 3 * 240 == -944
    assert Fraction(10) > Fraction(256, 27)
    result = {"schema": "exact-rational-module-orbits-v1", "field": "F2",
              "scope": "fixed companion factorization; exact grouped-channel word pressure",
              "random_finite_jet_span_checks": check_random_lattices(),
              "literal_channel_word_checks": word_checks,
              "orbits": records, "uniform_lower_bound_B": "10",
              "upper_bound_point": str(point), "upper_bound_B": str(rates[0] * rates[1]),
              "factoring_exponent_improvement": False}
    directory = Path(__file__).resolve().parent
    result["sources"] = {name: sha256((directory / name).read_bytes()).hexdigest()
                         for name in ("check_cancellation_module_automata.py",
                                      "check_minimal_boundary_modules.py")}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = check()
    payload = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
