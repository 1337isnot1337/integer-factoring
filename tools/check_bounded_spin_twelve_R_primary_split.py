#!/usr/bin/env python3
"""R-sign residues of the twelve-site primitive central idempotents.

The row and column cut trace forms are the ordinary and the replica-swap
twisted traces.  Combined, they split the rank residue of each central
idempotent into its R=+ and R=- parts modulo five.  This is a structural
congruence, not an exact rank computation.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from math import comb
from pathlib import Path

from bounded_spin_crossing import paths
from bounded_spin_six_vertex import monomial_entry

ROOT = Path(__file__).resolve().parents[1]
PRIME = 5
INVERSE_TWO = 3
TEN_SIGNS = [1, 1, 1, 1, 1, -1, 1]


def original_form(four_paths):
    n = len(four_paths[0])
    result = [0] * (n + 1)
    for (_, b, c), value in monomial_entry(four_paths).items():
        for j in range(2 * b + 1):
            result[2 * c + j] = (result[2 * c + j]
                                 + value * pow(-1, j) * comb(2 * b, j)) % PRIME
    return tuple(result)


def trace_vectors(length):
    ps = paths(length)
    row, column = Counter(), Counter()
    for a in ps:
        for b in ps:
            row[original_form((a, b, b, a))] += 1
            column[original_form((a, b, a, b))] += 1
    return row, column


def evaluate(idempotents, forms, basis):
    size = len(basis)
    vectors = [sum(count * form[k] for form, count in forms.items()) % PRIME
               for k in range(max(basis) + 1)]
    return [sum(int(e[j]) * vectors[basis[j]] for j in range(size)) % PRIME
            for e in idempotents]


def ten_site_control():
    model = json.loads((ROOT / 'evidence/bounded-spin-ten-algebra-model.json').read_text())
    row_forms, column_forms = trace_vectors(10)
    rows = evaluate(model['idempotents'], row_forms, model['basis_indices'])
    columns = evaluate(model['idempotents'], column_forms, model['basis_indices'])
    expected_rows = [d % PRIME for d in model['dimensions']]
    expected_columns = [s * d % PRIME for s, d in zip(TEN_SIGNS, model['dimensions'])]
    assert rows == expected_rows, (rows, expected_rows)
    assert columns == expected_columns, (columns, expected_columns)
    return {
        'dimensions': model['dimensions'],
        'row_residues': rows,
        'column_residues': columns,
        'R_signs': TEN_SIGNS,
    }


def check():
    control = ten_site_control()
    cut = json.loads((ROOT / 'evidence/bounded-spin-twelve-cut-trace-form.json').read_text())
    separator = json.loads((ROOT / 'evidence/bounded-spin-twelve-separator-split.json').read_text())
    idempotents = json.loads((ROOT / 'evidence/bounded-spin-twelve-coefficient-idempotents.json').read_text())
    assert cut['basis_indices'] == separator['basis_indices'] == idempotents['basis_indices']
    row = cut['idempotent_trace_residues']['row']
    column = cut['idempotent_trace_residues']['column']
    joint = separator['joint_characters_T_S']
    assert len(row) == len(column) == len(joint) == 8
    rows = []
    for i in range(8):
        t, _s = joint[str(i)]
        r_plus = (row[i] + column[i]) * INVERSE_TWO % PRIME
        r_minus = (row[i] - column[i]) * INVERSE_TWO % PRIME
        assert (r_plus + r_minus) % PRIME == row[i]
        assert (r_plus - r_minus) % PRIME == column[i]
        if r_plus == 0 and r_minus == 0:
            status = 'undetermined'
        elif r_plus == 0 or r_minus == 0:
            status = 'pure'
        else:
            status = 'mixed'
        rows.append({
            'idempotent_index': i,
            'T_eigenvalue': t,
            'S_eigenvalue': _s,
            'rank_residue': row[i],
            'R_trace_residue': column[i],
            'R_plus_residue': r_plus,
            'R_minus_residue': r_minus,
            'status': status,
        })
    mixed = [r['idempotent_index'] for r in rows if r['status'] == 'mixed']
    pure = [r['idempotent_index'] for r in rows if r['status'] == 'pure']
    undetermined = [r['idempotent_index'] for r in rows if r['status'] == 'undetermined']
    groups = {}
    for t in (0, 1, 4):
        idx = [r['idempotent_index'] for r in rows if r['T_eigenvalue'] == t]
        groups[str(t)] = {
            'idempotents': idx,
            'rank_residue_sum': sum(row[i] for i in idx) % PRIME,
            'R_plus_residue_sum': sum(r['R_plus_residue'] for r in rows if r['idempotent_index'] in idx) % PRIME,
            'R_minus_residue_sum': sum(r['R_minus_residue'] for r in rows if r['idempotent_index'] in idx) % PRIME,
        }
    return {
        'status': 'pass',
        'prime': PRIME,
        'ten_site_control': control,
        'twelve_site_primary_R_split': rows,
        'mixed_primaries': mixed,
        'pure_primaries_modulo_five': pure,
        'undetermined_primaries': undetermined,
        'T_group_sums': groups,
        'source_sha256': {
            'cut_trace_form': cut.get('source_sha256'),
            'separator_split': separator.get('source_sha256'),
        },
        'scope': ('Modulo-five R-sign split of each twelve-site primitive central '
                  'idempotent image, obtained from the ordinary and replica-swap '
                  'twisted cut trace residues. All seven nontrivial primaries are '
                  'mixed modulo five; idempotent 0 is undetermined because both '
                  'residues vanish. This does not give exact dimensions, exact '
                  'ranks, or a two-cut pairing, and no score change is claimed.'),
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
        'mixed_primaries': result['mixed_primaries'],
        'pure_primaries_modulo_five': result['pure_primaries_modulo_five'],
        'undetermined_primaries': result['undetermined_primaries'],
    }))


if __name__ == '__main__':
    main()
