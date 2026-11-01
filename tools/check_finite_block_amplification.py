#!/usr/bin/env python3
"""Exact finite-block amplification controls on complete restricted crossings.

These controls verify the amplification algebra, not a new exponent.
All tensor entries and factorizations are over F4, with no floating rank
or pressure calculation. The input tensors are independently decoded
from the frozen paired-restriction search.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement, permutations, product
import json
from math import factorial, isqrt, prod
from pathlib import Path

import numpy as np


MUL = np.asarray(((0, 0, 0, 0), (0, 1, 2, 3),
                  (0, 2, 3, 1), (0, 3, 1, 2)), dtype=np.uint8)
INV = (0, 1, 3, 2)


def multiply(first, second):
    assert first.shape[1] == second.shape[0]
    result = np.zeros((first.shape[0], second.shape[1]), dtype=np.uint8)
    for k in range(first.shape[1]):
        result ^= MUL[first[:, k, None], second[k, None, :]]
    return result


def factor(matrix):
    """Independent exact RREF factorization A=A[:,pivots] R[:rank,:]."""
    reduced = matrix.copy()
    pivots, row = [], 0
    for column in range(matrix.shape[1]):
        choices = np.flatnonzero(reduced[row:, column])
        if not len(choices):
            continue
        chosen = row + int(choices[0])
        reduced[[row, chosen]] = reduced[[chosen, row]]
        reduced[row] = MUL[INV[int(reduced[row, column])], reduced[row]]
        for other in range(matrix.shape[0]):
            if other != row and reduced[other, column]:
                reduced[other] ^= MUL[reduced[other, column], reduced[row]]
        pivots.append(column)
        row += 1
        if row == matrix.shape[0]:
            break
    first, second = matrix[:, pivots], reduced[:row, :]
    assert np.array_equal(multiply(first, second), matrix)
    return first, second


def labels(index, qx, qy):
    x = index % qx
    index //= qx
    y = index % qy
    index //= qy
    u = index % qy
    v = index // qy
    return x, y, u, v


def index_of(x, y, u, v, qx, qy):
    return x + qx * (y + qy * (u + qy * v))


def flatten(values, qx, qy, cut):
    result = np.zeros((qx * qy, qx * qy), dtype=np.uint8)
    for index, coefficient in enumerate(values):
        x, y, u, v = labels(index, qx, qy)
        row, column = ((x + qx * y, u + qy * v) if cut == 0
                       else (x + qx * u, y + qy * v))
        result[row, column] = coefficient
    return result


def endpoint(values, qx, qy):
    for index, coefficient in enumerate(values):
        x, y, u, v = labels(index, qx, qy)
        assert int(coefficient) == int(x == v and y == u)


def coefficient_atoms(record):
    seed = record['seed']
    digest = sha256(json.dumps(seed['values'], separators=(',', ':')).encode()).hexdigest()
    assert digest == record['tensor_sha256']
    result = []
    for degree in range(seed['degree'] + 1):
        values = np.asarray([((first >> degree) & 1) |
                             (((second >> degree) & 1) << 1)
                             for first, second in seed['values']], dtype=np.uint8)
        if np.any(values):
            result.append(values)
    total = np.zeros(len(seed['values']), dtype=np.uint8)
    for values in result:
        total ^= values
    endpoint(total, seed['q_x'], seed['q_u'])
    return result


def digits(value, base, length):
    result = []
    for _ in range(length):
        result.append(value % base)
        value //= base
    assert value == 0
    return result


def tensor_word(atoms, word, q):
    h, large = len(word), q ** len(word)
    cache = [digits(i, q, h) for i in range(large)]
    values = np.zeros(large ** 4, dtype=np.uint8)
    for index in range(len(values)):
        coordinates = [cache[v] for v in labels(index, large, large)]
        value = 1
        for position, atom in enumerate(word):
            local = index_of(*(coordinate[position] for coordinate in coordinates), q, q)
            value = int(MUL[value, atoms[atom][local]])
        values[index] = value
    return values


def tensor_factors(local_factors, word, q):
    h, large = len(word), q ** len(word)
    cache = [digits(i, q, h) for i in range(large)]
    ranks = [local_factors[atom][0].shape[1] for atom in word]
    size = prod(ranks)
    first = np.zeros((large ** 2, size), dtype=np.uint8)
    second = np.zeros((size, large ** 2), dtype=np.uint8)
    for middle, coordinates in enumerate(product(*(range(rank) for rank in ranks))):
        for outer in range(large ** 2):
            lower, upper = cache[outer % large], cache[outer // large]
            left = right = 1
            for position, atom in enumerate(word):
                local = lower[position] + q * upper[position]
                left = int(MUL[left, local_factors[atom][0][local, coordinates[position]]])
                right = int(MUL[right, local_factors[atom][1][coordinates[position], local]])
            first[outer, middle], second[middle, outer] = left, right
    return first, second


def histogram_query(atoms, composition, external, q, h):
    """Coefficient extraction, independent of the explicit word expansion."""
    large = q ** h
    coordinates = [digits(v, q, h) for v in labels(external, large, large)]
    zero = (0,) * len(atoms)
    state = {zero: 1}
    for position in range(h):
        local = index_of(*(coordinate[position] for coordinate in coordinates), q, q)
        following = {}
        for counts, value in state.items():
            for atom, vector in enumerate(atoms):
                if counts[atom] == composition[atom]:
                    continue
                new = list(counts)
                new[atom] += 1
                key = tuple(new)
                following[key] = following.get(key, 0) ^ int(MUL[value, vector[local]])
        state = following
    return state.get(composition, 0)


def rational_root_sum(products, denominator):
    numerator = 0
    for value in products:
        square = value * denominator ** 2
        root = isqrt(square)
        numerator += root + int(root * root < square)
    return numerator


def grouped_controls(record):
    atoms = coefficient_atoms(record)
    q, h = record['seed']['q_x'], 2
    assert q == record['seed']['q_u']
    factors = [[factor(flatten(atom, q, q, cut)) for atom in atoms] for cut in (0, 1)]
    ranks = [[first.shape[1] for first, _ in side] for side in factors]
    products = [ranks[0][i] * ranks[1][i] for i in range(len(atoms))]
    denominator = 1 << 20
    numerator = rational_root_sum(products, denominator)
    large, total = q ** h, np.zeros(q ** (4 * h), dtype=np.uint8)
    pieces, entries, queries, controls = [], 0, 0, 0
    for multiset in combinations_with_replacement(range(len(atoms)), h):
        words = sorted(set(permutations(multiset)))
        composition = tuple(Counter(multiset).get(i, 0) for i in range(len(atoms)))
        count = factorial(h) // prod(factorial(v) for v in composition)
        assert len(words) == count
        values = np.zeros(large ** 4, dtype=np.uint8)
        for word in words:
            values ^= tensor_word(atoms, word, q)
        total ^= values
        bounds, actual = [], []
        for cut in (0, 1):
            parts = [tensor_factors(factors[cut], word, q) for word in words]
            first = np.concatenate([pair[0] for pair in parts], axis=1)
            second = np.concatenate([pair[1] for pair in parts], axis=0)
            matrix = flatten(values, large, large, cut)
            assert np.array_equal(multiply(first, second), matrix)
            bound = count * prod(ranks[cut][i] ** composition[i] for i in range(len(atoms)))
            assert first.shape[1] == bound
            rank = factor(matrix)[0].shape[1]
            assert rank <= min(bound, large ** 2)
            bounds.append(bound)
            actual.append(rank)
            entries += matrix.size
            controls += 1
        assert bounds[0] * bounds[1] <= (len(atoms) ** 2 * max(products)) ** h
        assert bounds[0] * bounds[1] * denominator ** (2 * h) <= numerator ** (2 * h)
        for external, expected in enumerate(values):
            assert histogram_query(atoms, composition, external, q, h) == int(expected)
            queries += 1
        pieces.append({'composition': composition, 'words': count,
                       'rank_bounds': bounds, 'actual_ranks': actual})
    endpoint(total, large, large)
    return {'input_tensor_sha256': record['tensor_sha256'],
            'local_ranks': ranks, 'rank_products': products,
            'sqrt_sum_upper_bound': [numerator, denominator],
            'grouped_pieces': pieces,
            'counts': {'complete_factorizations': controls, 'factorization_entries': entries,
                       'histogram_queries': queries, 'endpoint_entries': len(total)}}


def pair_swap_controls(record):
    atoms = coefficient_atoms(record)
    qx, qy = record['seed']['q_x'], record['seed']['q_u']
    size = qx * qy
    ranks = [[factor(flatten(atom, qx, qy, cut))[0].shape[1] for atom in atoms]
             for cut in (0, 1)]
    total = np.zeros(size ** 4, dtype=np.uint8)
    checks = entries = 0
    for i, j in product(range(len(atoms)), repeat=2):
        values = np.zeros(size ** 4, dtype=np.uint8)
        for index in range(len(values)):
            x, y, u, v = labels(index, size, size)
            x0, xp = x % qx, x // qx
            yp, y0 = y % qx, y // qx
            up, u0 = u % qx, u // qx
            v0, vp = v % qx, v // qx
            first = atoms[i][index_of(x0, y0, u0, v0, qx, qy)]
            second = atoms[j][index_of(yp, xp, vp, up, qx, qy)]
            values[index] = MUL[first, second]
        total ^= values
        for cut in (0, 1):
            rank = factor(flatten(values, size, size, cut))[0].shape[1]
            assert rank == ranks[cut][i] * ranks[cut][j]
            checks += 1
        entries += len(values)
    endpoint(total, size, size)
    return {'input_tensor_sha256': record['tensor_sha256'], 'dimensions': [qx, qy],
            'local_ranks': ranks, 'paired_atoms': len(atoms) ** 2,
            'counts': {'complete_paired_cut_ranks': checks, 'paired_tensor_entries': entries,
                       'endpoint_entries': len(total)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=Path('evidence/finite-fiber-restriction-search.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    search = json.loads(raw)
    assert search['completed']
    chosen = [next(row for row in search['records'] if row['source_index'] == 0
                   and row['subset_x'] == [0, 2] and row['subset_y'] == second)
              for second in ([0, 2], [0, 1, 2])]
    report = {'scope': 'Exact construction controls only; no finite block below the factoring threshold is claimed.',
              'input_path': str(args.input), 'input_sha256': sha256(raw).hexdigest(),
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'tensor_square': grouped_controls(chosen[0]),
              'rectangular_pair_swap': pair_swap_controls(chosen[1]),
              'passed': True}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True,
                      'tensor_square': report['tensor_square']['counts'],
                      'rectangular_pair_swap': report['rectangular_pair_swap']['counts']}))


if __name__ == '__main__':
    main()
