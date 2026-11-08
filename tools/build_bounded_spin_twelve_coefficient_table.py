#!/usr/bin/env python3
"""Build the exact ten-dimensional twelve-site coefficient multiplication table."""
from __future__ import annotations

from hashlib import sha256
import json
from math import comb
from pathlib import Path
import random
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bounded_spin_crossing import paths
from bounded_spin_six_vertex import monomial_entry

ROOT = Path(__file__).resolve().parents[1]
PROBE = ROOT / 'evidence/bounded-spin-twelve-coefficient-dimension.json'
OUTPUT = ROOT / 'evidence/bounded-spin-twelve-coefficient-table.json'
PRIME = 5
LENGTH = 12
BASIS = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]
SOURCE_NAMES = ('bounded_spin_crossing.py', 'bounded_spin_six_vertex.py')


def original_form(four_paths: tuple) -> list[int]:
    result = [0] * (LENGTH + 1)
    for (_, b, c), value in monomial_entry(four_paths).items():
        for j in range(2 * b + 1):
            result[2 * c + j] = (result[2 * c + j]
                                 + value * comb(2 * b, j) * (-1) ** j) % PRIME
    return result


def rank_mod(matrix: list[list[int]]) -> int:
    a = [[int(x) % PRIME for x in row] for row in matrix]
    rows, columns = len(a), len(a[0])
    r = 0
    for column in range(columns):
        pivot = next((i for i in range(r, rows) if a[i][column]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inverse = pow(a[r][column], -1, PRIME)
        a[r] = [x * inverse % PRIME for x in a[r]]
        for i in range(rows):
            if i != r and a[i][column]:
                factor = a[i][column]
                a[i] = [(x - factor * y) % PRIME for x, y in zip(a[i], a[r])]
        r += 1
        if r == rows:
            break
    return r


def solve_mod(matrix: list[list[int]], rhs: list[int]) -> list[int]:
    a = [[int(x) % PRIME for x in row] for row in matrix]
    b = [int(x) % PRIME for x in rhs]
    rows, columns = len(a), len(a[0])
    augmented = [row + [value] for row, value in zip(a, b)]
    r = 0
    pivots = []
    for column in range(columns):
        pivot = next((i for i in range(r, rows) if augmented[i][column]), None)
        if pivot is None:
            continue
        augmented[r], augmented[pivot] = augmented[pivot], augmented[r]
        inverse = pow(augmented[r][column], -1, PRIME)
        augmented[r] = [x * inverse % PRIME for x in augmented[r]]
        for i in range(rows):
            if i != r and augmented[i][column]:
                factor = augmented[i][column]
                augmented[i] = [(x - factor * y) % PRIME
                                for x, y in zip(augmented[i], augmented[r])]
        pivots.append(column)
        r += 1
        if r == rows:
            break
    for i in range(r, rows):
        if not any(augmented[i][:columns]) and augmented[i][-1]:
            raise ValueError('inconsistent system')
    solution = [0] * columns
    for i, column in enumerate(pivots):
        solution[column] = augmented[i][-1] % PRIME
    return solution


def pair_vectors(ps, row_pair, column_pair):
    q = len(ps)
    row = []
    for a in range(q):
        for b in range(q):
            row.append(original_form((ps[row_pair[0]], ps[row_pair[1]],
                                      ps[a], ps[b])))
    column = []
    for a in range(q):
        for b in range(q):
            column.append(original_form((ps[a], ps[b],
                                         ps[column_pair[0]], ps[column_pair[1]])))
    return row, column


def main() -> None:
    ps = paths(LENGTH)
    q = len(ps)
    index = {tuple(path): i for i, path in enumerate(ps)}
    started = time.time()
    probe = json.loads(PROBE.read_text())
    samples = [tuple(index[tuple(path)] for path in quad)
               for quad in probe['witness_quadruples']]
    if len(samples) != 10:
        raise ValueError('the dimension witness must have ten entries')
    vectors = [original_form(tuple(ps[i] for i in ids)) for ids in samples]
    entry_matrix = [[vector[k] for k in BASIS] for vector in vectors]
    if rank_mod(entry_matrix) != 10:
        raise ValueError('the witness does not separate the coefficient span')

    rows, columns = [], []
    for ids in samples:
        row, column = pair_vectors(ps, (ids[0], ids[1]), (ids[2], ids[3]))
        rows.append(row)
        columns.append(column)
    print(json.dumps({'stage': 'witness rows and columns',
                      'seconds': round(time.time() - started, 1)}), flush=True)

    unit_rhs = [int((ids[0] * q + ids[1]) == (ids[2] * q + ids[3]))
                for ids in samples]
    unit = solve_mod(entry_matrix, unit_rhs)

    table = [[[0] * 10 for _ in range(10)] for _ in range(10)]
    for i in range(10):
        for j in range(10):
            left, right = BASIS[i], BASIS[j]
            rhs = []
            for l in range(10):
                rhs.append(sum(rows[l][s][left] * columns[l][s][right]
                               for s in range(q * q)) % PRIME)
            table[i][j] = solve_mod(entry_matrix, rhs)

    commutative = all(table[i][j] == table[j][i]
                      for i in range(10) for j in range(10))
    associative = True
    for i in range(10):
        for j in range(10):
            for k in range(10):
                left = [0] * 10
                for a in range(10):
                    if table[i][j][a]:
                        factor = table[i][j][a]
                        for b in range(10):
                            left[b] = (left[b] + factor * table[a][k][b]) % PRIME
                right = [0] * 10
                for a in range(10):
                    if table[j][k][a]:
                        factor = table[j][k][a]
                        for b in range(10):
                            right[b] = (right[b] + factor * table[i][a][b]) % PRIME
                if left != right:
                    associative = False
    regular = [[[table[i][j][k] for j in range(10)] for k in range(10)]
               for i in range(10)]
    unit_ok = True
    for i in range(10):
        product = [sum(regular[i][k][j] * unit[j] for j in range(10)) % PRIME
                   for k in range(10)]
        if product != [int(k == i) for k in range(10)]:
            unit_ok = False
    if not (commutative and associative and unit_ok):
        raise ValueError('the computed table failed an algebra identity')

    rng = random.Random(99)
    holdouts = []
    for _ in range(2):
        row_pair = (rng.randrange(q), rng.randrange(q))
        column_pair = (rng.randrange(q), rng.randrange(q))
        row, column = pair_vectors(ps, row_pair, column_pair)
        direct = original_form((ps[row_pair[0]], ps[row_pair[1]],
                                ps[column_pair[0]], ps[column_pair[1]]))
        values = []
        for i in range(10):
            line = []
            for j in range(10):
                value = sum(row[s][BASIS[i]] * column[s][BASIS[j]]
                            for s in range(q * q)) % PRIME
                predicted = sum(table[i][j][k] * direct[BASIS[k]]
                                for k in range(10)) % PRIME
                if value != predicted:
                    raise ValueError('holdout entry does not match the table')
                line.append(value)
            values.append(line)
        holdouts.append({'row_pair': list(row_pair),
                         'column_pair': list(column_pair),
                         'product_values': values})
        print(json.dumps({'stage': 'holdout', 'index': len(holdouts),
                          'seconds': round(time.time() - started, 1)}), flush=True)

    evidence = {
        'passed': True,
        'prime': PRIME,
        'length': LENGTH,
        'basis_indices': BASIS,
        'witness_quadruples': [[list(path) for path in quad]
                               for quad in probe['witness_quadruples']],
        'entry_matrix': entry_matrix,
        'unit': unit,
        'table': table,
        'regular': regular,
        'commutative': commutative,
        'associative': associative,
        'unit_ok': unit_ok,
        'holdouts': holdouts,
        'probe_sha256': sha256(PROBE.read_bytes()).hexdigest(),
        'source_sha256': {
            name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
            for name in SOURCE_NAMES
        },
        'seconds': round(time.time() - started, 1),
        'scope': ('Exact multiplication table of the ten-dimensional twelve-site '
                  'coefficient span in the entry-witness basis, with unit, '
                  'commutativity, associativity, regular representation, and two '
                  'independent holdout checks. No primary decomposition, cut-rank '
                  'bound, or score change is claimed.'),
    }
    OUTPUT.write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps({'passed': True, 'dimension': 10,
                      'unit': unit, 'seconds': evidence['seconds']}))


if __name__ == '__main__':
    main()
