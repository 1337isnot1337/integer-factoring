#!/usr/bin/env python3
"""Exact full rank-one-image cubic span and column rank reduction.

The physical control is a literal four-state Pfaffian contraction. Its local
edge convention and cubic input are from PR #29, pinned in the report. The
normal form and rank formula are independently implemented here.
"""

import argparse
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import random


def poly_mul(a, b):
    out = 0
    while b:
        bit = b & -b
        out ^= a << (bit.bit_length() - 1)
        b ^= bit
    return out


@lru_cache(None)
def local_pfaffian(x, y, u, v):
    vertices = [(leg, bit) for leg, mask in enumerate((u, x, v, y))
                for bit in range(2) if mask >> bit & 1]
    if len(vertices) % 2:
        return 0
    edges = [[0] * len(vertices) for _ in vertices]
    for i, (first, a) in enumerate(vertices):
        for j in range(i + 1, len(vertices)):
            second, b = vertices[j]
            if (first, second) == (0, 1) and a == b:
                edges[i][j] = 2  # p, with r=1
            elif (first, second) == (1, 3) and a != b:
                edges[i][j] = 16  # q, stored at bit degree four
            elif ((first, second) in ((1, 2), (2, 3), (0, 3)) and a == b
                  or (first, second) == (0, 2) and a != b):
                edges[i][j] = 1

    @lru_cache(None)
    def recurse(mask):
        if not mask:
            return 1
        first_bit = mask & -mask
        i = first_bit.bit_length() - 1
        rest, out = mask ^ first_bit, 0
        while rest:
            bit = rest & -rest
            j = bit.bit_length() - 1
            out ^= poly_mul(edges[i][j], recurse(mask ^ first_bit ^ bit))
            rest ^= bit
        return out

    return recurse((1 << len(vertices)) - 1)


class PhysicalCubic:
    def __init__(self, d=7):
        self.d = d
        self.degrees = (0, 1, 2, 4, 12)

        def occupation(i, j, high=False):
            return tuple((2 if high else 1) if k == i else 3 if k == j else 0
                         for k in range(d))

        self.inputs = [[occupation(i, j) for j in range(d) if j != i]
                       for i in range(d)]
        self.outputs = [occupation(i, (i + 1) % d, high)
                        for high in (False, True) for i in range(d)]
        self.local = tuple(local_pfaffian(*q) for q in product(range(4), repeat=4))
        self.parent = [[[[0] * (2 * d) for _ in range(d)] for _ in range(d)]
                       for _ in self.degrees]
        for x, u, y, v in product(range(2 * d), range(2 * d), range(d), range(d)):
            value = self.kernel(self.outputs[x], self.outputs[u], y, v)
            scalars = ((1 << 0) | (1 << 8), (1 << 1) | (1 << 9),
                       (1 << 2) | (1 << 6), (1 << 4) | (1 << 8),
                       (1 << 8) | (1 << 12))
            reconstructed = 0
            for degree, scalar in zip(self.degrees, scalars):
                if value >> degree & 1:
                    reconstructed ^= scalar
            assert value == reconstructed
            for k, degree in enumerate(self.degrees):
                if value >> degree & 1:
                    self.parent[k][y][v][x] ^= 1 << u

    @lru_cache(None)
    def pattern(self, indices):
        out = 1
        for i in indices:
            out = poly_mul(out, self.local[i])
            if not out:
                break
        return out

    def kernel(self, first, second, y, v):
        parity = sum(((a.bit_count() + b.bit_count()) % 2) << i
                     for i, (a, b) in enumerate(zip(first, second)))
        if parity ^ (1 << y) ^ (1 << v):
            return 0
        value = 0
        for yy, vv in product(self.inputs[y], self.inputs[v]):
            indices = tuple(sorted(((a * 4 + b) * 4 + c) * 4 + e
                                   for a, b, c, e in zip(first, yy, second, vv)))
            value ^= self.pattern(indices)
        return value

    def matrices(self, complement):
        d = self.d
        supports = [[i] + [d + j for j in range(d) if complement[i] >> j & 1]
                    for i in range(d)]
        output = [(1 << i) | (complement[i] << d) for i in range(d)]
        rows = [[0] * (d * d) for _ in self.degrees]
        columns = [[0] * (d * d) for _ in self.degrees]
        for k, y, v, x in product(range(5), range(d), range(d), range(d)):
            value = 0
            for i in supports[x]:
                value ^= self.parent[k][y][v][i]
            for u in range(d):
                if (value & output[u]).bit_count() % 2:
                    rows[k][d * x + y] ^= 1 << (d * v + u)
                    columns[k][d * x + u] ^= 1 << (d * v + y)
        return rows, columns


def add_basis(basis, value):
    while value:
        pivot = value.bit_length() - 1
        if pivot not in basis:
            basis[pivot] = value
            return True
        value ^= basis[pivot]
    return False


def flat(rows):
    return sum(row << (len(rows) * i) for i, row in enumerate(rows))


def gf_mul(a, b, field):
    if field == 2:
        return a & b
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return out


def rank_dense(matrix, field, reverse=False):
    if not matrix:
        return 0
    rows = [row[:] for row in matrix]
    if reverse:
        rows.reverse()
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        value = rows[rank][column]
        inv = next(i for i in range(1, field) if gf_mul(value, i, field) == 1)
        rows[rank] = [gf_mul(x, inv, field) for x in rows[rank]]
        for i in range(rank + 1, len(rows)):
            scalar = rows[i][column]
            if scalar:
                rows[i] = [x ^ gf_mul(scalar, y, field)
                           for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def pairs(d):
    return [(i, i) for i in range(d)] + list(combinations(range(d), 2))


def cycle(i, j, d):
    return int(i != j and (i - j) % d in (1, d - 1))


def build_cuts(d, e, a, q):
    labels = pairs(d)
    lookup = {pair: i for i, pair in enumerate(labels)}
    rows = [[0] * (d * d) for _ in range(d * d)]
    columns = [[0] * (d * d) for _ in range(d * d)]
    for x, y, v, u in product(range(d), repeat=4):
        value = q[lookup[tuple(sorted((x, u)))]][lookup[tuple(sorted((y, v)))]]
        if x == v and y == u:
            value ^= e ^ (a if cycle(x, y, d) else 0)
        rows[d * x + y][d * v + u] = value
        columns[d * x + u][d * v + y] = value
    return rows, columns


def reduced_column_rank(d, e, a, q, field):
    labels = pairs(d)
    zeros = [k for k, (i, j) in enumerate(labels) if i != j
             and (e ^ (a if cycle(i, j, d) else 0)) == 0]
    keep = list(range(d)) + zeros
    matrix = [[q[i][j] ^ (e if i == j and i < d else 0) for j in keep]
              for i in keep]
    nonzero_edges = d * (d - 1) // 2 - len(zeros)
    return 2 * nonzero_edges + rank_dense(matrix, field), len(matrix)


def coordinates(binary_columns, d):
    """Decode e, a and Q, checking every entry of the proposed normal form."""
    edge = 1
    nonedge = 2
    e = ((binary_columns[nonedge] >> nonedge) & 1) ^ (
        (binary_columns[nonedge] >> (nonedge * d)) & 1)
    edge_value = ((binary_columns[edge] >> edge) & 1) ^ (
        (binary_columns[edge] >> (edge * d)) & 1)
    a = edge_value ^ e
    labels = pairs(d)
    q = [[((binary_columns[d * x + u] >> (d * v + y)) & 1)
          ^ ((e ^ (a if cycle(x, y, d) else 0)) if x == v and y == u else 0)
          for y, v in labels] for x, u in labels]
    dense = build_cuts(d, e, a, q)
    rebuilt = [[sum(value << j for j, value in enumerate(row)) for row in cut]
               for cut in dense]
    assert rebuilt[1] == binary_columns
    coordinate = e | (a << 1)
    for i, row in enumerate(q):
        coordinate |= sum(value << (2 + len(labels) * i + j) for j, value in enumerate(row))
    return coordinate, rebuilt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    d = 7
    physical = PhysicalCubic(d)
    zero, symmetry, endpoint = 0, 0, 0
    for k, y, v, i, j in product(range(5), range(d), range(d), range(d), range(d)):
        parent = physical.parent[k][y][v]
        assert not (parent[i] >> (d + j)) & 1
        assert not (parent[d + i] >> j) & 1
        zero += 2
        assert ((parent[d + i] >> (d + j)) & 1) == (
            (physical.parent[k][v][y][d + j] >> (d + i)) & 1)
        symmetry += 1

    major_symmetry, q3_parent_entries = 0, 0
    for k, y, v, x, u in product(range(5), range(d), range(d), range(2 * d), range(2 * d)):
        value = (physical.parent[k][y][v][x] >> u) & 1
        assert value == ((physical.parent[k][v][y][u] >> x) & 1)
        major_symmetry += 1
        if k == 4:
            assert value == int(x == d + y and u == d + v)
            q3_parent_entries += 1

    baseline = physical.matrices([0] * d)
    for x, y, v, u in product(range(d), repeat=4):
        crossing = int(x == v and y == u)
        noncrossing = int(x == u and y == v)
        diagonal = int(x == y == v == u)
        graph = crossing & cycle(x, y, d)
        expected = (crossing, graph ^ noncrossing, noncrossing ^ diagonal, graph, 0)
        for k, value in enumerate(expected):
            assert ((baseline[0][k][d * x + y] >> (d * v + u)) & 1) == value

    square_vectors = [1 << i for i in range(d)] + [(1 << i) | (1 << j)
                                                           for i, j in combinations(range(d), 2)]
    full_basis, q3_basis = {}, {}
    matrix_entries = 0
    maps = [None] + list(product(square_vectors, repeat=2))
    for map_index, uv in enumerate(maps):
        complement = [0] * d if uv is None else [uv[1] if uv[0] >> i & 1 else 0
                                               for i in range(d)]
        rows, columns = physical.matrices(complement)
        assert rows[0] == columns[0] == [1 << i for i in range(d * d)]
        endpoint += d ** 4
        for k, (row, column) in enumerate(zip(rows, columns)):
            coordinate, rebuilt = coordinates(column, d)
            assert rebuilt == [row, column]
            add_basis(full_basis, coordinate)
            matrix_entries += 2 * d ** 4
            if k == 4 and uv is not None:
                assert coordinate & 3 == 0
                add_basis(q3_basis, coordinate)
                u, v = uv
                for x, y, z, w in product(range(d), repeat=4):
                    expected = ((u >> x) & 1) & ((u >> w) & 1) & ((v >> y) & 1) & ((v >> z) & 1)
                    assert ((row[x * d + y] >> (z * d + w)) & 1) == expected

    assert len(q3_basis) == 784
    assert len(full_basis) == 786

    # The q^3 formula checked above is quadratic in all entries of C. One
    # matrix unit or the sum of two matrix units has rank at most two.
    # Their vectorized outer squares are a basis of the full symmetric space.
    major_basis = {}
    width = d * d
    rank_two_maps = [1 << i for i in range(width)] + [(1 << i) | (1 << j)
                                                     for i, j in combinations(range(width), 2)]
    for vector in rank_two_maps:
        outer = [vector if vector >> i & 1 else 0 for i in range(width)]
        assert add_basis(major_basis, flat(outer))
    assert len(major_basis) == width * (width + 1) // 2 == 1225
    rng = random.Random(309786)
    rank_controls, profiles = 0, []
    for field in (2, 4):
        for e, a in product(range(field), repeat=2):
            for trial in range(6):
                q = [[rng.randrange(field) for _ in range(28)] for _ in range(28)]
                if trial < 2:
                    left = [[rng.randrange(field) for _ in range(trial + 1)] for _ in range(28)]
                    right = [[rng.randrange(field) for _ in range(28)] for _ in range(trial + 1)]
                    q = [[0] * 28 for _ in range(28)]
                    for i, j, k in product(range(28), range(28), range(trial + 1)):
                        q[i][j] ^= gf_mul(left[i][k], right[k][j], field)
                _, column = build_cuts(d, e, a, q)
                actual = rank_dense(column, field)
                assert actual == rank_dense(column, field, reverse=True)
                reduced, size = reduced_column_rank(d, e, a, q, field)
                assert actual == reduced
                rank_controls += 1
                if trial == 0:
                    profiles.append({'field': field, 'e': e, 'a': a,
                                     'residual_matrix_dimension': size,
                                     'rank': actual})

    result = {
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'prior_source': '1337isnot1337/integer-factoring@7154cec09b982b54a0fcc0bb7a303318fc1f90da',
        'prior_packages': ['cubic-varying-image-affordable-escape.md',
                           'cubic-varying-image-normal-form-and-peels.md',
                           'cubic-fixed-image-rank-one-barrier.md'],
        'labels': d, 'field_scope': 'Every field of characteristic two',
        'rank_one_maps': 784, 'zero_map_included': True,
        'square_vector_masks': square_vectors,
        'symmetric_update_dimension': len(q3_basis),
        'complete_raw_span_dimension': len(full_basis),
        'physical_raw_matrices_regenerated': 5 * len(maps),
        'normal_form_entries_checked_both_cuts': matrix_entries,
        'endpoint_entries_checked': endpoint,
        'parent_mixed_zero_entries_checked': zero,
        'parent_swap_symmetry_entries_checked': symmetry,
        'parent_major_symmetry_entries_checked': major_symmetry,
        'q3_parent_formula_entries_checked': q3_parent_entries,
        'complete_five_polynomial_parent_entries_checked': (2 * d) ** 2 * d ** 2,
        'all_complement_span_dimension': len(major_basis),
        'rank_at_most_two_maps_suffice': len(rank_two_maps),
        'exact_column_rank_controls': rank_controls,
        'column_rank_control_profiles': profiles,
        'rank_formula': '2 times nonzero off-diagonal graph weights plus rank of Q restricted to diagonals and zero-weight pairs, with e added to its diagonal-label diagonal',
        'strict_seed_obtained': False,
        'optimality_claim': False,
        'scope': 'All rank-one complementary maps give the 786-dimensional normal form; all complementary maps give the 1225-dimensional major-symmetric space and rank-at-most-two maps already span it. The fixed seven-label cubic input/output selector is retained. Both descriptions hold over every characteristic-two field. No factoring exponent improvement is claimed.',
    }
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ('column_rank_control_profiles', 'square_vector_masks')}))


if __name__ == '__main__':
    main()
