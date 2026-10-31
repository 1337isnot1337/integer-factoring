#!/usr/bin/env python3
"""Exact controls for one-sided charge pruning and both-cut Smith monotonicity."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random

from weighted_smith_oracle import WeightedSmithOracle, polynomial_multiply


def rank_binary(matrix):
    basis = {}
    for row in matrix:
        value = sum((entry & 1) << j for j, entry in enumerate(row))
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def flatten(tensor, q, cut):
    matrix = [[0] * (q * q) for _ in range(q * q)]
    for (x, y, u, v), value in tensor.items():
        row, column = ((q*x+y, q*u+v) if cut == 0 else (q*x+u, q*y+v))
        matrix[row][column] = value
    return matrix


def weights(gauges, q, cut):
    row_legs, column_legs = (((0, 1), (2, 3)) if cut == 0 else ((0, 2), (1, 3)))
    return [[gauges[legs[0]][i] + gauges[legs[1]][j]
             for i in range(q) for j in range(q)] for legs in (row_legs, column_legs)]


def scalar_tensor(q, a, b, rng):
    base, changed = {}, {}
    for labels in product(range(q), repeat=4):
        x, y, u, v = labels
        charge = a*x + b*y - b*u - a*v
        value = int(x == v and y == u)
        if charge == 0:
            value ^= polynomial_multiply(3, rng.randrange(4))
        base[labels] = value
        changed[labels] = value ^ (polynomial_multiply(3, rng.randrange(4)) if charge > 0 else 0)
        assert base[labels].bit_count() % 2 == int(x == v and y == u)
        assert changed[labels].bit_count() % 2 == int(x == v and y == u)
    return base, changed


def word_entry(tensor, words):
    value = 1
    for labels in zip(*words):
        value = polynomial_multiply(value, tensor[labels])
    return value


def matrix_multiply(a, b):
    return tuple(tuple(polynomial_multiply(a[i][0], b[0][j]) ^
                       polynomial_multiply(a[i][1], b[1][j]) for j in range(2)) for i in range(2))


def auxiliary_controls(rng):
    q, a, b = 3, 1, -1
    zero = ((0, 0), (0, 0))
    identity = ((1, 0), (0, 1))
    base, changed = {}, {}
    for labels in product(range(q), repeat=4):
        x, y, u, v = labels
        charge = a*x + b*y - b*u - a*v
        original = identity if x == v and y == u else zero
        base[labels] = tuple(tuple(original[i][j] ^ (polynomial_multiply(3, rng.randrange(4))
                                  if charge == 0 else 0) for j in range(2)) for i in range(2))
        changed[labels] = tuple(tuple(base[labels][i][j] ^ (polynomial_multiply(3, rng.randrange(4))
                                     if charge > 0 else 0) for j in range(2)) for i in range(2))
    neutral = [base[key] for key in base if key[0]-key[1]+key[2]-key[3] == 0]
    assert any(matrix_multiply(x, y) != matrix_multiply(y, x) for x in neutral for y in neutral)
    words = sorted(set(permutations(range(q))))
    checks = 0
    for boundary in product(words, repeat=4):
        first = second = identity
        for labels in zip(*boundary):
            first = matrix_multiply(first, base[labels])
            second = matrix_multiply(second, changed[labels])
        assert first == second
        checks += 1
    return checks


def quartic_counts():
    signs = []
    copies = []
    for i in range(16):
        leg, species = i // 4, i % 2
        signs.append(1 if ((leg in (1, 3)) == (species == 0)) else -1)
        copies.append((i % 4) // 2)
    neutral = 0
    positive, negative = Counter(), Counter()
    for subset in combinations(range(16), 4):
        vector = tuple(sum(signs[i] for i in subset if copies[i] == copy) for copy in range(2))
        neutral += vector == (0, 0)
        if sum(vector) == 4:
            positive[vector] += 1
        if sum(vector) == -4:
            negative[vector] += 1
    pairs = sum(n * negative[tuple(-v for v in vector)] for vector, n in positive.items())
    assert neutral == 328 and pairs == 1810
    assert sorted(positive.values()) == [1, 1, 16, 16, 36]
    return {'neutral_per_copy': neutral, 'opposite_vector_pairs': pairs,
            'positive_copy_occupation_counts': [positive[(k, 4-k)] for k in range(5)]}


def check():
    rng = random.Random(73100021)
    prefixes = pressure_checks = rank_checks = word_checks = 0
    for q in (2, 3):
        for a, b in product((-1, 1), repeat=2):
            for sample in range(4):
                base, changed = scalar_tensor(q, a, b, rng)
                for cut in (0, 1):
                    original = flatten(base, q, cut)
                    perturbed = flatten(changed, q, cut)
                    for infinity in (False, True):
                        gauges = [[rng.randrange(-4, 5) for _ in range(q)] for _ in range(4)]
                        rows, columns = weights(gauges, q, cut)
                        nu, _ = WeightedSmithOracle(original, infinity).profile(rows, columns)
                        mu, _ = WeightedSmithOracle(perturbed, infinity).profile(rows, columns)
                        assert sum(nu) == sum(mu)
                        for k in range(1, q*q+1):
                            assert sum(mu[:k]) <= sum(nu[:k])
                            prefixes += 1
                        for s in (Fraction(2), Fraction(3, 2), Fraction(2, 3)):
                            assert sum(s**(-x) for x in mu) >= sum(s**(-x) for x in nu)
                            pressure_checks += 1
                    degree = max(value.bit_length() for row in perturbed for value in row)
                    for functional in range(1 << degree):
                        r0 = rank_binary([[(value & functional).bit_count() % 2 for value in row] for row in original])
                        r1 = rank_binary([[(value & functional).bit_count() % 2 for value in row] for row in perturbed])
                        assert r1 >= r0
                        rank_checks += 1
                if sample == 0:
                    words = sorted(set(permutations(range(q))))
                    for boundary in product(words, repeat=4):
                        assert word_entry(base, boundary) == word_entry(changed, boundary)
                        word_checks += 1
    original = WeightedSmithOracle([[1, 0], [0, 1]]).profile([0, 0], [0, 0])[0]
    bidirectional = WeightedSmithOracle([[1, 3], [3, 1]]).profile([0, 0], [0, 0])[0]
    assert original == (0, 0) and bidirectional == (0, 2)
    assert sum(Fraction(2)**(-x) for x in bidirectional) < sum(Fraction(2)**(-x) for x in original)
    return {'status': 'pass', 'scalar_tensor_samples': 32,
            'both_cut_both_place_prefix_comparisons': prefixes,
            'exact_convex_pressure_comparisons': pressure_checks,
            'coefficient_functional_rank_comparisons': rank_checks,
            'common_type_polynomial_entries': word_checks,
            'noncommuting_auxiliary_word_entries': auxiliary_controls(rng),
            'quartic_enumeration': quartic_counts(),
            'opposite_direction_control_profiles': [list(original), list(bidirectional)],
            'source_sha256': {p.name: sha256(p.read_bytes()).hexdigest() for p in
                              (Path(__file__), Path(__file__).with_name('weighted_smith_oracle.py'))},
            'scope': 'One-sided perturbations only; no improved factoring exponent.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rendered = json.dumps(check(), indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
