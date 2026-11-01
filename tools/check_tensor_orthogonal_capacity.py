#!/usr/bin/env python3
"""Exact controls for the general all-subset Gram capacity theorems.

The proofs are in research/tensor-orthogonal-capacity.md. Finite checks
confirm witnesses, binary classification, and the subspace identity;
they do not replace the infinite-family proofs.
"""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path


class Field:
    def __init__(self, degree, modulus):
        self.degree, self.q, self.modulus = degree, 1 << degree, modulus
        self.table = [[self._multiply(a, b) for b in range(self.q)]
                      for a in range(self.q)]
        self.inverses = [0] + [next(b for b in range(1, self.q)
                                  if self.table[a][b] == 1)
                              for a in range(1, self.q)]

    def _multiply(self, a, b):
        out = 0
        while b:
            if b & 1:
                out ^= a
            a <<= 1
            if a & self.q:
                a ^= self.modulus
            b >>= 1
        return out

    def rank(self, rows):
        if not rows:
            return 0
        a = [list(row) for row in rows]
        rank = 0
        for j in range(len(a[0])):
            pivot = next((i for i in range(rank, len(a)) if a[i][j]), None)
            if pivot is None:
                continue
            a[rank], a[pivot] = a[pivot], a[rank]
            inv = self.inverses[a[rank][j]]
            a[rank] = [self.table[inv][x] for x in a[rank]]
            for i in range(len(a)):
                if i != rank and a[i][j]:
                    scale = a[i][j]
                    a[i] = [x ^ self.table[scale][y]
                            for x, y in zip(a[i], a[rank])]
            rank += 1
            if rank == len(a):
                break
        return rank


def capacity_upper(degree, rank):
    q = 1 << degree
    value = Fraction(degree * rank * (q ** rank - 1),
                     q ** rank - q ** (rank - 1))
    return value.numerator // value.denominator


def subspace_upper(degree, rank, k):
    q = 1 << degree
    value = Fraction(degree * (rank * (rank + 1) - k * (k + 1))
                     * (q ** rank - 1), 2 * (q ** rank - q ** k))
    return value.numerator // value.denominator


def axes(degree, rank):
    return [[(1 << b) if i == j else 0
             for j in range(rank) for b in range(degree)]
            for i in range(rank)]


def direct_sum(*blocks):
    sizes = [len(b[0]) for b in blocks]
    out, offset = [], 0
    for block, width in zip(blocks, sizes):
        for row in block:
            out.append([0] * offset + list(row)
                       + [0] * (sum(sizes) - offset - width))
        offset += width
    return out


def subset_grams(field, matrix):
    r, m = len(matrix), len(matrix[0])
    columns = list(zip(*matrix))
    outers = [[field.table[c[i]][c[j]] for i in range(r) for j in range(r)]
              for c in columns]
    grams = [[0] * (r * r)]
    for mask in range(1, 1 << m):
        bit = mask & -mask
        j = bit.bit_length() - 1
        grams.append([a ^ b for a, b in zip(grams[mask ^ bit], outers[j])])
    for mask, gram in enumerate(grams):
        selected = [[matrix[i][j] for j in range(m) if mask >> j & 1]
                    for i in range(r)]
        gram_rows = [gram[i * r:(i + 1) * r] for i in range(r)]
        assert field.rank(selected) == field.rank(gram_rows), (matrix, mask)
    return grams


def projective_vectors(q, rank):
    for first in range(rank):
        for tail in product(range(q), repeat=rank - first - 1):
            yield (0,) * first + (1,) + tail


def hyperplane_controls(field, matrix, grams):
    r, m = len(matrix), len(matrix[0])
    histogram = Counter()
    columns = list(zip(*matrix))
    checked = 0
    for normal in projective_vectors(field.q, r):
        in_plane = 0
        for column in columns:
            dot = 0
            for a, b in zip(normal, column):
                dot ^= field.table[a][b]
            in_plane += dot == 0
        count = 0
        for gram in grams:
            image_inside = True
            for j in range(r):
                dot = 0
                for i in range(r):
                    dot ^= field.table[normal[i]][gram[i * r + j]]
                if dot:
                    image_inside = False
                    break
            count += image_inside
            checked += 1
        assert count == 1 << in_plane
        assert in_plane >= m - field.degree * r
        histogram[in_plane] += 1
    count_hyperplanes = (field.q ** r - 1) // (field.q - 1)
    incidence = (field.q ** (r - 1) - 1) // (field.q - 1)
    assert sum(histogram.values()) == count_hyperplanes
    assert sum(k * v for k, v in histogram.items()) == m * incidence
    return {'hyperplanes': count_hyperplanes,
            'column_count_histogram': dict(sorted(histogram.items())),
            'subset_image_tests': checked}


def binary_rank(vectors):
    basis = {}
    for v in vectors:
        while v:
            pivot = v.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = v
                break
            v ^= basis[pivot]
    return len(basis)


def binary_classification(columns):
    basis, remaining = [], []
    for column in columns:
        if binary_rank(basis + [column]) > len(basis):
            basis.append(column)
        else:
            remaining.append(column)
    coordinates = {0: 0}
    for i, column in enumerate(basis):
        coordinates.update({value ^ column: coefficient | (1 << i)
                            for value, coefficient in list(coordinates.items())})
    used = 0
    for column in remaining:
        support = coordinates[column]
        if not support or support.bit_count() % 2 or support & used:
            return False
        used |= support
    return True


def binary_exhaustion(ambient_rank):
    vectors = list(range(1, 1 << ambient_rank))
    m = len(vectors)
    outers = [sum(((v >> i & 1) * v) << (ambient_rank * i)
                  for i in range(ambient_rank)) for v in vectors]
    grams = [0] * (1 << m)
    bad = bytearray(1 << m)
    ranks = bytearray(1 << m)
    for mask in range(1, 1 << m):
        bit = mask & -mask
        grams[mask] = grams[mask ^ bit] ^ outers[bit.bit_length() - 1]
        columns = [vectors[i] for i in range(m) if mask >> i & 1]
        ranks[mask] = binary_rank(columns)
        rows = [(grams[mask] >> (ambient_rank * i)) & ((1 << ambient_rank) - 1)
                for i in range(ambient_rank)]
        bad[mask] = binary_rank(rows) != ranks[mask]
    # Boolean subset zeta transform: bad[S] records a failing subset of S.
    for i in range(m):
        bit = 1 << i
        for mask in range(1 << m):
            if mask & bit:
                bad[mask] |= bad[mask ^ bit]
    valid_by_rank_size = Counter()
    maxima = {}
    for mask in range(1 << m):
        columns = [vectors[i] for i in range(m) if mask >> i & 1]
        classification = binary_classification(columns)
        assert classification == (not bad[mask]), (ambient_rank, columns)
        if classification:
            r, size = ranks[mask], mask.bit_count()
            valid_by_rank_size[r, size] += 1
            maxima[r] = max(maxima.get(r, 0), size)
    assert all(value == 3 * r // 2 for r, value in maxima.items())
    return {'ambient_rank': ambient_rank, 'all_coefficient_sets': 1 << m,
            'maxima_by_rank': maxima,
            'valid_counts': [{'rank': r, 'columns': size, 'count': count}
                             for (r, size), count in sorted(valid_by_rank_size.items())]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-orthogonal-capacity.json'))
    args = parser.parse_args()
    fields = {f: Field(f, modulus) for f, modulus in
              ((1, 0b11), (2, 0b111), (3, 0b1011), (4, 0b10011))}
    for field in fields.values():
        for a, b, c in product(range(field.q), repeat=3):
            assert field.table[a][b ^ c] == field.table[a][b] ^ field.table[a][c]
            assert field.table[field.table[a][b]][c] == field.table[a][field.table[b][c]]
    binary = [binary_exhaustion(r) for r in range(1, 5)]
    f4_five = [[1, 1, 2, 2, 0], [0, 1, 3, 1, 3]]
    f4_one = [[1, 2]]
    f8_uniform = [[1, 1, 2, 2, 4, 4], [0, 1, 4, 6, 6, 2]]
    cases = [('f2_triangle', 1, [[1, 0, 1], [0, 1, 1]]),
             ('f4_rank1', 2, f4_one), ('f4_rank2', 2, f4_five),
             ('f4_rank3', 2, direct_sum(f4_five, f4_one)),
             ('f4_rank4', 2, direct_sum(f4_five, f4_five)),
             ('f8_rank2_axes', 3, axes(3, 2)),
             ('f8_rank2_uniform', 3, f8_uniform),
             ('f16_rank2_axes', 4, axes(4, 2))]
    witnesses = []
    for name, degree, matrix in cases:
        field = fields[degree]
        r, m = len(matrix), len(matrix[0])
        assert field.rank(matrix) == r
        grams = subset_grams(field, matrix)
        assert len({tuple(gram) for gram in grams}) == 1 << m
        upper = capacity_upper(degree, r)
        assert m == upper
        record = {'name': name, 'field_degree': degree, 'field_modulus': field.modulus,
                  'matrix': matrix, 'rank': r, 'columns': m, 'exact_upper': upper,
                  'direct_subset_gram_checks': len(grams),
                  'hyperplane_identity': hyperplane_controls(field, matrix, grams)}
        if name == 'f8_rank2_uniform':
            assert all(field.rank([[matrix[i][j] for j in (a, b)] for i in range(2)]) == 2
                       for a in range(6) for b in range(a))
            # 2^m sum_S 2^{-rank(S)}, with uniform rank min(2, |S|).
            dimension = sum(Fraction(1 << 6, 1 << min(2, mask.bit_count()))
                            for mask in range(1 << 6))
            assert dimension == 1168
            record['uniform_matroid'] = 'U(2,6)'
            record['full_paired_kernel_dimension'] = int(dimension)
        witnesses.append(record)
    f4_rank6 = direct_sum(f4_five, f4_five, f4_five)
    assert len(f4_rank6[0]) == capacity_upper(2, 6) == 15
    derived = {'name': 'f4_rank6', 'matrix': f4_rank6, 'rank': 6, 'columns': 15,
               'exact_upper': 15, 'validity': 'direct sum of three checked f4_rank2 blocks',
               'direct_subset_gram_checks': 0}
    bounds = []
    for degree in range(1, 17):
        q = 1 << degree
        assert capacity_upper(degree, 2) == (3 if degree == 1 else 5 if degree == 2 else 2 * degree)
        for r in range(1, 33):
            upper = capacity_upper(degree, r)
            assert upper >= degree * r
            if q - 1 >= degree * r:
                assert upper == degree * r
            if r > 1:
                assert upper == subspace_upper(degree, r, r - 1)
            bounds.append({'field_degree': degree, 'rank': r, 'axis_lower': degree * r,
                           'hyperplane_upper': upper})
    report = {'passed': True, 'scope': __doc__, 'binary_exhaustions': binary,
              'directly_checked_witnesses': witnesses, 'derived_witnesses': [derived],
              'integer_bound_controls': bounds,
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'binary_sets_checked': sum(c['all_coefficient_sets'] for c in binary),
                      'witness_subset_gram_checks': sum(c['direct_subset_gram_checks'] for c in witnesses),
                      'hyperplane_subset_image_tests': sum(c['hyperplane_identity']['subset_image_tests'] for c in witnesses),
                      'derived_f4_rank6_capacity': 15, 'output': str(args.output)}))


if __name__ == '__main__':
    main()
