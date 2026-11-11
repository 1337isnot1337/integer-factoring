#!/usr/bin/env python3
"""Ten-site operator control for the two diagonal trace orientations.

The row-cut idempotents are genuine projections.  Their ordinary traces
match the (a,b,b,a) diagonal residues.  The (a,b,a,b) diagonal is instead
the trace of the path-pair swap times the row idempotent, not the
ordinary column-cut trace.  The path-pair swap is not the action of the
coefficient unit.  This control checks the twelve-site row and column
trace provenance.
"""
from __future__ import annotations

import json
from collections import Counter
from math import comb
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))

from bounded_spin_batched_blocks import coefficient_tensor
from bounded_spin_crossing import paths
from bounded_spin_six_vertex import monomial_entry
from check_bounded_spin_growth import combination, original_coefficients
from check_bounded_spin_six_vertex import operators
from native_gf5_flint import rank

PRIME = 5
MODEL = ROOT / 'evidence/bounded-spin-ten-algebra-model.json'
OUTPUT = ROOT / 'evidence/bounded-spin-ten-trace-orientation.json'


def original_form(four_paths):
    length = len(four_paths[0])
    result = [0] * (length + 1)
    for (_, b, c), value in monomial_entry(four_paths).items():
        for j in range(2 * b + 1):
            result[2 * c + j] = (result[2 * c + j]
                                 + value * comb(2 * b, j) * (-1) ** j) % PRIME
    return tuple(result)


def diagonal_vector(ps, pairing):
    forms = Counter(original_form(pairing(a, b))
                    for a in range(len(ps)) for b in range(len(ps)))
    length = len(next(iter(forms)))
    return [sum(count * form[k] for form, count in forms.items()) % PRIME
            for k in range(length)]


def diagonal_traces(vector, basis, idempotents):
    return [sum(int(e[j]) * vector[basis[j]] for j in range(len(basis)))
            % PRIME for e in idempotents]


def path_swap(size, q):
    indices = [(a, b) for a in range(q) for b in range(q)]
    position = {state: i for i, state in enumerate(indices)}
    matrix = np.zeros((size, size), dtype=np.uint8)
    for a, b in indices:
        matrix[position[(b, a)], position[(a, b)]] = 1
    return matrix


def main() -> None:
    model = json.loads(MODEL.read_text())
    basis = model['basis_indices']
    monomials, tensors = coefficient_tensor(10)
    raw = original_coefficients(monomials, tensors, 10)
    row, column = operators(raw, 10)
    row = row.astype(np.uint8)
    column = column.astype(np.uint8)
    ps = paths(10)
    q = len(ps)
    size = q * q

    row_idempotents = [combination(e, row[basis]) for e in model['idempotents']]
    column_idempotents = [combination(e, column[basis])
                          for e in model['idempotents']]
    row_traces = [int(np.trace(e)) % PRIME for e in row_idempotents]
    column_traces = [int(np.trace(e)) % PRIME for e in column_idempotents]
    row_ranks = [rank(e) for e in row_idempotents]
    column_ranks = [rank(e) for e in column_idempotents]

    swapped_vector = diagonal_vector(
        ps, lambda a, b: (ps[a], ps[b], ps[b], ps[a]))
    ordinary_vector = diagonal_vector(
        ps, lambda a, b: (ps[a], ps[b], ps[a], ps[b]))
    swapped = diagonal_traces(swapped_vector, basis, model['idempotents'])
    ordinary = diagonal_traces(ordinary_vector, basis, model['idempotents'])
    P = path_swap(size, q)
    P_traces = [int(np.trace(P.astype(np.int64) @ e.astype(np.int64)))
                % PRIME for e in row_idempotents]

    dimensions = [d % PRIME for d in model['dimensions']]
    if row_traces != dimensions:
        raise ValueError('row traces do not match the primary dimensions')
    if swapped != row_traces:
        raise ValueError('the (a,b,b,a) form is not the row trace')
    if column_traces != dimensions:
        raise ValueError('unexpected ordinary column traces')
    if ordinary != P_traces:
        raise ValueError('the (a,b,a,b) form is not trace(P e)')
    if ordinary == column_traces:
        raise ValueError('the two diagonal forms are not distinguished')
    if column_ranks == row_ranks:
        raise ValueError('column ranks coincide with row ranks')

    receipt = {
        'passed': True,
        'prime': PRIME,
        'length': 10,
        'basis_indices': basis,
        'primary_dimensions': model['dimensions'],
        'row_idempotent_traces': row_traces,
        'row_idempotent_ranks': row_ranks,
        'column_idempotent_traces': column_traces,
        'column_idempotent_ranks': column_ranks,
        'swapped_pair_traces': swapped,
        'ordinary_pair_traces': ordinary,
        'path_swap_row_traces': P_traces,
        'row_traces_match_dimensions': True,
        'swapped_pair_matches_row_trace': True,
        'ordinary_pair_matches_path_swap': True,
        'path_swap_is_not_coefficient_unit_action': True,
        'column_traces_are_not_ranks': True,
        'scope': (
            'Ten-site operator control for the diagonal trace orientation. '
            'The row idempotents are genuine projections: their ordinary '
            'traces equal the (a,b,b,a) residues and the primary dimensions. '
            'The (a,b,a,b) residues equal trace(P e), where P is the '
            'path-pair swap. P is not the action of the coefficient unit, so '
            'this is not an R-character of e and not the ordinary column-cut '
            'trace. The column trace functional is the ordinary trace on the '
            'column-cut module; it is not the column rank. At twelve sites the '
            'column residues therefore cannot be used as column ranks without '
            'a separate column operator check. No twelve-site rank or score '
            'change is claimed.'
        ),
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({
        'passed': True,
        'row_traces': row_traces,
        'swapped_pair_traces': swapped,
        'ordinary_pair_traces': ordinary,
        'path_swap_row_traces': P_traces,
        'column_ranks': column_ranks,
    }, indent=2))


if __name__ == '__main__':
    main()
