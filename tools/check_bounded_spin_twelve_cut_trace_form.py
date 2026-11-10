#!/usr/bin/env python3
"""Cut trace form of the twelve-site coefficient algebra.

The diagonal trace forms at (a,b,b,a) and (a,b,a,b) give the row and column
trace vectors on the cut module.  Combined with the exact coefficient table,
they determine the trace form of every product and the traces of the eight
central idempotents modulo five.  This is a checksum and a structural
invariant, not an exact rank computation.
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from math import comb
from pathlib import Path

from bounded_spin_crossing import paths
from bounded_spin_six_vertex import monomial_entry
from build_bounded_spin_twelve_coefficient_idempotents import nullspace
from native_gf5_flint import rank

ROOT = Path(__file__).resolve().parents[1]
PRIME = 5
LENGTH = 12


def original_form(four_paths):
    n = len(four_paths[0])
    result = [0] * (n + 1)
    for (_, b, c), value in monomial_entry(four_paths).items():
        for j in range(2 * b + 1):
            result[2 * c + j] = (result[2 * c + j]
                                 + value * comb(2 * b, j) * (-1) ** j) % PRIME
    return tuple(result)


def trace_forms():
    ps = paths(LENGTH)
    row, column = Counter(), Counter()
    for a in ps:
        for b in ps:
            row[original_form((a, b, b, a))] += 1
            column[original_form((a, b, a, b))] += 1
    return row, column


def trace_vector(forms, size):
    return [sum(count * form[k] for form, count in forms.items()) % PRIME
            for k in range(size)]


def trace_form(table, vector):
    return [[sum(table[i][j][k] * vector[k] for k in range(len(vector))) % PRIME
             for j in range(len(vector))] for i in range(len(vector))]


def regular_trace_vector(table):
    return [sum(table[i][j][j] for j in range(len(table))) % PRIME
            for i in range(len(table))]


def multiply(table, left, right):
    size = len(left)
    result = [0] * size
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    factor = a * b % PRIME
                    row = table[i][j]
                    for k in range(size):
                        result[k] = (result[k] + factor * row[k]) % PRIME
    return result


def span_basis(vectors):
    span = []
    for vector in vectors:
        if rank(span + [vector]) > len(span):
            span.append(vector)
    return span


def check():
    started = time.time()
    table_data = json.loads((ROOT / 'evidence/bounded-spin-twelve-coefficient-table.json').read_text())
    idempotent_data = json.loads((ROOT / 'evidence/bounded-spin-twelve-coefficient-idempotents.json').read_text())
    table = table_data['table']
    basis = table_data['basis_indices']
    size = len(basis)
    row_forms, column_forms = trace_forms()
    row_vector = trace_vector(row_forms, size)
    column_vector = trace_vector(column_forms, size)
    row_form = trace_form(table, row_vector)
    column_form = trace_form(table, column_vector)
    regular_vector = regular_trace_vector(table)

    assert row_form == [[row_form[j][i] for j in range(size)] for i in range(size)]
    assert column_form == [[column_form[j][i] for j in range(size)] for i in range(size)]

    radical = nullspace(row_form)
    assert rank(row_form) == 9
    assert len(radical) == 1
    radical_vector = radical[0]
    radical_index = next(k for k in range(size) if radical_vector[k] != 0)
    assert all(radical_vector[k] == 0 for k in range(size) if k != radical_index)
    assert radical_index == 8
    radical_label = basis[radical_index]
    assert radical_label == 9
    assert regular_vector[radical_index] != 0

    regular_radical = idempotent_data['trace_form']['radical_basis']
    assert rank(regular_radical + [radical_vector]) == 3

    rt_data = json.loads((ROOT / 'evidence/bounded-spin-twelve-RT-coefficient.json').read_text())
    assert rt_data['passed']
    assert rt_data['basis_indices'] == basis
    unit = [1] + [0] * (size - 1)
    class_sum = rt_data['class_sum_coordinates']
    powers = [unit, class_sum]
    for _ in range(3):
        powers.append(multiply(table, powers[-1], class_sum))
    rt_span = span_basis(powers)
    assert len(rt_span) == rt_data['class_sum_subalgebra_dimension'] == 5
    radical_direction = [0] * size
    radical_direction[radical_index] = 1
    extension_dimension = rank(rt_span + [radical_direction])
    assert extension_dimension == 6
    assert rt_data['extensions_by_one_coefficient'][str(radical_label)] == extension_dimension
    assert [radical_label] not in rt_data['completing_coefficient_pairs']

    idempotents = idempotent_data['primitive_central_idempotents']
    row_traces = [sum(e['idempotent'][j] * row_vector[j] for j in range(size)) % PRIME
                  for e in idempotents]
    column_traces = [sum(e['idempotent'][j] * column_vector[j] for j in range(size)) % PRIME
                     for e in idempotents]

    return {
        'status': 'pass',
        'length': LENGTH,
        'prime': PRIME,
        'basis_indices': basis,
        'trace_vectors': {'row': row_vector, 'column': column_vector},
        'trace_forms': {'row': row_form, 'column': column_form},
        'trace_form_rank': {'row': rank(row_form), 'column': rank(column_form)},
        'trace_radical': {'row': radical, 'column': nullspace(column_form)},
        'cut_radical_vector': radical_vector,
        'cut_radical_basis_position': radical_index,
        'cut_radical_coefficient_index': radical_label,
        'cut_radical_is_coefficient_9': True,
        'regular_trace_vector': regular_vector,
        'regular_trace_radical': regular_radical,
        'cut_radical_outside_regular': True,
        'rt_subalgebra_dimension': len(rt_span),
        'cut_radical_in_rt_subalgebra': rank(rt_span + [radical_direction]) == len(rt_span),
        'cut_radical_extension_dimension': extension_dimension,
        'cut_radical_in_completing_pairs': [radical_label] in rt_data['completing_coefficient_pairs'],
        'idempotent_trace_residues': {'row': row_traces, 'column': column_traces},
        'idempotent_regular_ranks': [e['regular_rank'] for e in idempotents],
        'seconds': round(time.time() - started, 1),
        'scope': ('Cut trace form of the twelve-site central-idempotent algebra. '
                  'The row and column forms have rank 9 with a one-dimensional '
                  'radical spanned by coefficient direction C_9 (basis position 8); '
                  'this direction has nonzero regular trace and is outside the '
                  'regular trace radical. It lies outside the five-dimensional '
                  'coefficient subalgebra <C_0,T> and extends it only to dimension '
                  '6, so it is not a minimal completing direction. The idempotent '
                  'trace residues are ranks modulo five. Exact ranks and any score '
                  'change remain open.'),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text + '\n')
    print(json.dumps({
        'status': result['status'],
        'trace_form_rank': result['trace_form_rank'],
        'cut_radical_vector': result['cut_radical_vector'],
        'idempotent_trace_residues': result['idempotent_trace_residues'],
    }))


if __name__ == '__main__':
    main()
