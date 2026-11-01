#!/usr/bin/env python3
"""Independent second-coefficient controls and labels for the existing row model."""

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
from random import Random

from bounded_spin_crossing import paths, recursive_vector, path_norm, polynomial_entry

ROOT = Path(__file__).resolve().parents[1]


def transpositions(selected):
    vectors = [recursive_vector(path) for path in selected]
    operators = []
    for i, j in combinations(range(len(selected[0])), 2):
        matrix = []
        for row in vectors:
            values = []
            for path, column in zip(selected, vectors):
                total = 0
                for word, coefficient in column.items():
                    target = word
                    if ((word >> i) ^ (word >> j)) & 1:
                        target ^= (1 << i) | (1 << j)
                    total += row.get(target, 0) * coefficient
                values.append(total * pow(path_norm(path), -1, 5) % 5)
            matrix.append(values)
        operators.append(matrix)
    return operators


def polynomial_product_mod(left, right, modulus):
    assert modulus[-1] == 1
    degree = len(modulus) - 1
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % 5
    for i in range(len(result) - 1, degree - 1, -1):
        coefficient = result[i]
        for j, a in enumerate(modulus):
            result[i - degree + j] = (result[i - degree + j] - coefficient * a) % 5
    return (result + [0] * degree)[:degree]


def evaluate(polynomial, value):
    return sum(c * pow(value, i, 5) for i, c in enumerate(polynomial)) % 5


def derivative(polynomial, value):
    return sum(i*c*pow(value, i-1, 5) for i, c in enumerate(polynomial) if i) % 5


def check():
    rng, controls = Random(20261103), []
    for length, limit in ((4, None), (6, 64), (8, 24)):
        selected = paths(length)
        size, k = len(selected), length // 2
        operators = transpositions(selected)
        for a, b in product(range(size), repeat=2):
            assert sum(t[a][b] for t in operators) % 5 == k * (k - 2) * int(a == b) % 5
        entries = list(product(range(size), repeat=4)) if limit is None else [
            tuple(rng.randrange(size) for _ in range(4)) for _ in range(limit)]
        for a, b, c, d in entries:
            expected = (sum(t[a][c] * t[b][d] for t in operators)
                        + k * (k + 1) * int(a == c and b == d)) * pow(3, -1, 5) % 5
            actual = polynomial_entry(tuple(selected[i] for i in (a, b, c, d)))[2]
            assert actual == expected
        controls.append({'length': length, 'path_count': size,
                         'literal_single_copy_class_sum_controls': size**2,
                         'constructor_C2_entry_controls': len(entries)})

    model_path = ROOT / 'evidence/bounded-spin-rank-model.json'
    certificate = json.loads(model_path.read_text())
    model = certificate['models'][0]
    basis = certificate['basis_indices']
    swap = model['polynomial_basis'][basis.index(0)]
    second = model['polynomial_basis'][basis.index(2)]
    modulus = model['minimal_polynomial']
    assert polynomial_product_mod(swap, swap, modulus) == [1] + [0] * (len(modulus) - 2)
    diagonal = [3*c % 5 for c in polynomial_product_mod(swap, second, modulus)]
    # For length eight, k(k+1)=4*5 vanishes in F5.
    labels = []
    for block in model['blocks']:
        assert block['factor'][1] == 1 and block['exponent'] in (1, 2)
        root = -block['factor'][0] % 5
        sv, tv = evaluate(swap, root), evaluate(diagonal, root)
        td = derivative(diagonal, root)
        assert sv in (1, 4)
        if block['exponent'] == 2:
            assert derivative(swap, root) == 0 and td != 0
        dimension = block['nullities'][-1]
        labels.append({'dimension': dimension, 'replica_swap': sv,
                       'diagonal_transposition_sum': tv,
                       'class_sum_derivative': td if block['exponent'] == 2 else None,
                       'primary_exponent': block['exponent'],
                       'nilpotent_rank': dimension - block['nullities'][1]})
    assert len({(b['replica_swap'], b['diagonal_transposition_sum']) for b in labels}) == len(labels)
    assert sum(b['dimension'] for b in labels) == 169
    return {'scope': 'Uniform C2 formula controls and polynomial labels of the existing eight-site row certificate; no uniform rank bound.',
            'seed': 20261103, 'entry_controls': controls,
            'eight_site_row_primary_labels': labels,
            'joint_generators_cover_row_model': True,
            'constructor_sha256': sha256((ROOT / 'tools/bounded_spin_crossing.py').read_bytes()).hexdigest(),
            'rank_model_sha256': sha256(model_path.read_bytes()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    result = check()
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
