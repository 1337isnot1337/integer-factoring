#!/usr/bin/env python3
"""Construct exact finite duals, or collisions, for six interval jet queries."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from walsh_cauchy_permanents import (
    principal_permanent_polynomials, restricted_interval_features,
)


def feature_matrix(length, precision):
    if length >= 4:
        return restricted_interval_features(length, precision)
    features = [(-1, 0), (-1, 1), (1, 1)]
    selected = [(1 << k) - 1 for k in range(length + 1)]
    values = {center: principal_permanent_polynomials(length, center, precision, selected)
              for center in (-1, 1)}
    rows = [[(values[center][subset] + [0] * precision)[order] for subset in selected]
            for center, order in features]
    return features, rows


def parity_certificate(matrix, precision):
    """Use finite Smith operations and return directly replayable evidence."""
    modulus = 1 << precision
    work = [[value % modulus for value in row] for row in matrix]
    rows, columns = len(work), len(work[0])
    row_operations = [{i: 1} for i in range(rows)]
    transform = [[int(i == j) for j in range(columns)] for i in range(columns)]
    target = [1] * columns
    orders = []
    for index in range(min(rows, columns)):
        best = None
        for row in range(index, rows):
            for column in range(index, columns):
                value = work[row][column]
                if value:
                    order = (value & -value).bit_length() - 1
                    if best is None or order < best[0]:
                        best = order, row, column
        if best is None:
            break
        order, row, column = best
        work[index], work[row] = work[row], work[index]
        row_operations[index], row_operations[row] = row_operations[row], row_operations[index]
        for values in work:
            values[index], values[column] = values[column], values[index]
        for values in transform:
            values[index], values[column] = values[column], values[index]
        target[index], target[column] = target[column], target[index]
        inverse = pow(work[index][index] >> order, -1, modulus)
        work[index] = [value * inverse % modulus for value in work[index]]
        row_operations[index] = {key: value * inverse % modulus
                                 for key, value in row_operations[index].items()}
        for row in range(index + 1, rows):
            factor = work[row][index] >> order
            if factor:
                work[row] = [(a - factor * b) % modulus for a, b in zip(work[row], work[index])]
                for key, value in row_operations[index].items():
                    row_operations[row][key] = (row_operations[row].get(key, 0) - factor * value) % modulus
        for column in range(index + 1, columns):
            factor = work[index][column] >> order
            target[column] = (target[column] - factor * target[index]) % modulus
            for values in transform:
                values[column] = (values[column] - factor * values[index]) % modulus
            work[index][column] = 0
        orders.append(order)

    tail = next((index for index in range(len(orders), columns) if target[index] & 1), None)
    if tail is not None:
        collision = [values[tail] for values in transform]
        assert sum(collision) & 1
        assert all(sum(a * b for a, b in zip(row, collision)) % modulus == 0 for row in matrix)
        return {'certified': False, 'precision': precision, 'nonzero_smith_orders': orders,
                'collision': {str(k): value for k, value in enumerate(collision) if value},
                'target_parity': 1}

    loss = max(order for order, value in zip(orders, target) if value & 1)
    final_modulus = 1 << (loss + 1)
    dual = {}
    for index, (order, value) in enumerate(zip(orders, target)):
        if value & 1:
            factor = 1 << (loss - order)
            for row, coefficient in row_operations[index].items():
                dual[row] = (dual.get(row, 0) + factor * coefficient) % final_modulus
    dual = {row: coefficient for row, coefficient in dual.items() if coefficient}
    assert all(sum(coefficient * matrix[row][column] for row, coefficient in dual.items())
               % final_modulus == 1 << loss for column in range(columns))
    return {'certified': True, 'precision': loss + 1, 'division_loss': loss,
            'nonzero_smith_orders': orders, 'dual': dual}


def generate(lengths):
    records = []
    for length in lengths:
        if length < 2 or length & (length - 1):
            raise ValueError('Packet lengths must be powers of two at least two.')
        precision = length.bit_length() + 2
        features, matrix = feature_matrix(length, precision)
        record = parity_certificate(matrix, precision)
        if record['certified']:
            record['dual'] = [{'center': features[row][0], 'order': features[row][1],
                               'coefficient': coefficient}
                              for row, coefficient in record['dual'].items()]
        record.update({'L': length, 'features': [list(feature) for feature in features]})
        records.append(record)
    root = Path(__file__).resolve().parent.parent
    sources = ['tools/search_interval_permanent_decoder.py', 'tools/walsh_cauchy_permanents.py']
    return {'status': 'exact finite certificates and explicit limited-query collisions',
            'records': records,
            'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest() for name in sources},
            'scope': 'Positive records decode the entire interval principal-permanent span. '
                     'Negative records concern only the specified six functionals. '
                     'No all-size precision theorem or general packet representation is claimed.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lengths', type=int, nargs='+', default=[2, 4, 8, 16, 32, 64, 128])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = generate(args.lengths)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'],
                      'records': [{'L': record['L'], 'certified': record['certified'],
                                   'precision': record['precision']}
                                  for record in result['records']]}, indent=2))


if __name__ == '__main__':
    main()
