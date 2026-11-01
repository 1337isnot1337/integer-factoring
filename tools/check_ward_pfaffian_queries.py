#!/usr/bin/env python3
"""Independent controls for streamed Ward Pfaffian entry queries."""

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
from itertools import product
import json
from pathlib import Path
import random

from ward_pfaffian_queries import BinaryField, pfaffian, squarefree_coefficient_stream


def direct_coefficient(forms, sites, field):
    values = [0] * (1 << forms)
    values[0] = 1
    for column, ordinary, inserted in sites:
        following = [field.mul(ordinary, value) for value in values]
        for mask, value in enumerate(values):
            if not value:
                continue
            for a, weight in enumerate(column):
                if not (mask >> a) & 1:
                    following[mask | (1 << a)] ^= field.mul(value, field.mul(inserted, weight))
        values = following
    return values[-1]


def recursive_pfaffian(matrix, field):
    @lru_cache(None)
    def visit(indices):
        if not indices:
            return 1
        value = 0
        first = indices[0]
        for at in range(1, len(indices)):
            remaining = indices[1:at] + indices[at + 1:]
            value ^= field.mul(matrix[first][indices[at]], visit(remaining))
        return value
    return visit(tuple(range(len(matrix))))


def scalar_controls():
    rng = random.Random(882031)
    counts = Counter()

    def check(forms, sites, field, kind):
        streamed = squarefree_coefficient_stream(forms, iter(sites), field)
        expected = direct_coefficient(forms, sites, field)
        assert streamed.value == expected, (field.size, forms, sites, streamed, expected)
        counts[kind] += 1
        counts['nonzero_queries'] += bool(expected)
        counts['queries_with_forced_sites'] += bool(streamed.forced_sites)
        counts['odd_residual_queries'] += streamed.residual_forms % 2
        if streamed.zero_reason:
            counts[streamed.zero_reason] += 1

    binary = BinaryField(3)
    for forms in range(4):
        choices = [(column, m, b) for column in product(range(2), repeat=forms)
                   for m, b in product(range(2), repeat=2)]
        for pairs in range(3):
            for sites in product(choices, repeat=pairs):
                check(forms, sites, binary, 'exhaustive_binary_queries')
    choices = [(column, m, b) for column in product(range(2), repeat=2)
               for m, b in product(range(2), repeat=2)]
    for sites in product(choices, repeat=3):
        check(2, sites, binary, 'exhaustive_binary_queries')

    for modulus in (3, 7, 11, 19):
        field = BinaryField(modulus)
        for value in range(1, field.size):
            assert field.mul(value, field.inverse(value)) == 1
        for forms in range(9):
            for sample in range(60):
                pairs = rng.randrange(13) if sample % 2 else forms + rng.randrange(5)
                forced = rng.randrange(min(forms, pairs) + 1)
                sites = []
                for i in range(pairs):
                    column = [rng.randrange(field.size) for _ in range(forms)]
                    if i < forced and sample % 3 == 0:
                        column = [int(a == i) for a in range(forms)]
                    ordinary = 0 if i < forced else rng.randrange(1, field.size)
                    inserted = rng.randrange(field.size)
                    sites.append((column, ordinary, inserted))
                rng.shuffle(sites)
                check(forms, sites, field, 'random_extension_field_queries')
        for size in range(0, 12, 2):
            for _ in range(40):
                matrix = [[0] * size for _ in range(size)]
                for i in range(size):
                    for j in range(i + 1, size):
                        matrix[i][j] = matrix[j][i] = rng.randrange(field.size)
                assert pfaffian(matrix, field) == recursive_pfaffian(matrix, field)
                counts['paired_elimination_vs_matching_controls'] += 1
    assert counts['dependent forced columns'] and counts['too many forced sites']
    assert counts['zero local factor'] and counts['nonzero_queries']
    return dict(counts)


def matrix_inverse(matrix, field):
    size = len(matrix)
    values = [list(row) + [int(i == j) for j in range(size)]
              for i, row in enumerate(matrix)]
    for col in range(size):
        pivot = next(i for i in range(col, size) if values[i][col])
        values[col], values[pivot] = values[pivot], values[col]
        reciprocal = field.inverse(values[col][col])
        values[col] = [field.mul(reciprocal, x) for x in values[col]]
        for i in range(size):
            if i != col:
                value = values[i][col]
                values[i] = [x ^ field.mul(value, y) for x, y in zip(values[i], values[col])]
    return [row[size:] for row in values]


def interpolation_controls():
    rng = random.Random(240820)
    field = BinaryField(19)
    cases = points = coefficient_checks = zero_ordinary_points = 0

    def evaluate(poly, x, y):
        value = 0
        for (a, b), scalar in poly.items():
            value ^= field.mul(scalar, field.mul(field.power(x, a), field.power(y, b)))
        return value

    def polynomial_add_product(target, first, second, weight=1):
        for (a, b), value in first.items():
            for (c, d), other in second.items():
                key = a + c, b + d
                target[key] = target.get(key, 0) ^ field.mul(weight, field.mul(value, other))

    for pairs in range(1, 4):
        for forms in range(pairs + 1):
            for _ in range(3):
                sites = []
                for i in range(pairs):
                    column = [rng.randrange(2) for _ in range(forms)]
                    ordinary = {key: rng.randrange(2) for key in ((0, 0), (1, 0), (0, 1))}
                    inserted = {key: rng.randrange(2) for key in ((0, 0), (1, 0), (0, 1), (1, 1))}
                    sites.append((column, ordinary, inserted))
                direct = [{(0, 0): 1}] + [{} for _ in range((1 << forms) - 1)]
                for column, ordinary, inserted in sites:
                    following = [{} for _ in direct]
                    for mask, poly in enumerate(direct):
                        polynomial_add_product(following[mask], poly, ordinary)
                        for a, weight in enumerate(column):
                            if not (mask >> a) & 1:
                                polynomial_add_product(following[mask | (1 << a)], poly, inserted, weight)
                    direct = following
                degree = pairs + forms
                nodes = list(range(degree + 1))
                vandermonde = [[field.power(node, exponent) for exponent in range(degree + 1)]
                               for node in nodes]
                inverse = matrix_inverse(vandermonde, field)
                grid = []
                for x in nodes:
                    row = []
                    for y in nodes:
                        evaluated = [(a, evaluate(m, x, y), evaluate(b, x, y)) for a, m, b in sites]
                        row.append(squarefree_coefficient_stream(forms, iter(evaluated), field).value)
                        points += 1
                        zero_ordinary_points += any(not m for _, m, _ in evaluated)
                    grid.append(row)
                first = [[0] * len(nodes) for _ in nodes]
                for a in nodes:
                    for y in nodes:
                        for x in nodes:
                            first[a][y] ^= field.mul(inverse[a][x], grid[x][y])
                for a in nodes:
                    for b in nodes:
                        value = 0
                        for y in nodes:
                            value ^= field.mul(first[a][y], inverse[b][y])
                        assert value == direct[-1].get((a, b), 0)
                        assert value in (0, 1)
                        coefficient_checks += 1
                cases += 1
    return {'polynomial_query_cases': cases, 'interpolation_points': points,
            'coefficient_comparisons': coefficient_checks,
            'grid_points_with_zero_ordinary_entries': zero_ordinary_points}


def physical_kernel_controls():
    from check_ward_stable_operators import bits, local_operators
    from explore_complement_edge_pencils import MONOMIALS
    from pfaffian_particle_hole import particle_hole_mask

    rng = random.Random(670881)
    field = BinaryField(7)
    incidence, creation, *_ = local_operators()
    max_degree = max(sum(term) for terms in MONOMIALS for term in terms)
    assert max_degree == 4
    checked = nonzero = 0
    for pairs in (1, 2, 3):
        for sample in range(50 if pairs < 3 else 8):
            forms = rng.randrange(pairs + 2)
            weights = [[rng.randrange(field.size) for _ in range(forms)] for _ in range(pairs)]
            eta = rng.randrange(field.size)
            edges = [rng.randrange(field.size) for _ in range(6)]
            ordinary = [[0] * 16 for _ in range(16)]
            for x, y, u, v in product(range(4), repeat=4):
                value = 0
                for term in MONOMIALS[((x * 4 + y) * 4 + u) * 4 + v]:
                    monomial = 1
                    for edge, exponent in zip(edges, term):
                        monomial = field.mul(monomial, field.power(edge, exponent))
                    value ^= monomial
                ordinary[particle_hole_mask(x, y)][particle_hole_mask(v, u)] = value
            acted = [[0] * 16 for _ in range(16)]
            for source in range(16):
                for target in bits(incidence[source]):
                    for col in range(16):
                        acted[target][col] ^= ordinary[source][col]
                for target in bits(creation[source]):
                    for col in range(16):
                        acted[target][col] ^= field.mul(eta, ordinary[source][col])
            source = rng.randrange(16 ** pairs)
            target = rng.randrange(16 ** pairs)
            vector = [1]
            for i in range(pairs):
                digit = source >> (4 * i) & 15
                vector = [field.mul(old, ordinary[row][digit])
                          for row in range(16) for old in vector]
            for a in range(forms):
                following = [0] * len(vector)
                for at, value in enumerate(vector):
                    if not value:
                        continue
                    for i, column in enumerate(weights):
                        if not column[a]:
                            continue
                        digit = at >> (4 * i) & 15
                        base = at ^ (digit << (4 * i))
                        scalar = field.mul(value, column[a])
                        for other in bits(incidence[digit]):
                            following[base | (other << (4 * i))] ^= scalar
                        for other in bits(creation[digit]):
                            following[base | (other << (4 * i))] ^= field.mul(eta, scalar)
                vector = following
            if sample % 2 == 0:
                support = [i for i, value in enumerate(vector) if value]
                if support:
                    target = rng.choice(support)
            sites = []
            for i, column in enumerate(weights):
                row, col = target >> (4 * i) & 15, source >> (4 * i) & 15
                sites.append((column, ordinary[row][col], acted[row][col]))
            query = squarefree_coefficient_stream(forms, iter(sites), field).value
            assert query == vector[target]
            checked += 1
            nonzero += bool(query)
    return {'literal_global_operator_entries': checked, 'nonzero_literal_entries': nonzero,
            'maximum_local_pfaffian_total_degree': max_degree}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    names = (Path(__file__).name, 'ward_pfaffian_queries.py',
             'check_ward_stable_operators.py', 'explore_complement_edge_pencils.py',
             'pfaffian_particle_hole.py')
    result = {'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in names},
              'scalar_controls': scalar_controls(),
              'interpolation_controls': interpolation_controls(),
              'physical_kernel_controls': physical_kernel_controls(),
              'scope': 'Exact physical-entry and formal coefficient queries in characteristic two; no coefficient-rank or factoring-score claim.'}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
