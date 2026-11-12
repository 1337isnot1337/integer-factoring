#!/usr/bin/env python3
"""Modulo-five separator splitting of the twelve-site cut coarse blocks.

The idempotent trace functional residues from the cut trace form give
nonzero rank residues for the joint T,S primaries, and direct diagonal
matrix entries certify the remaining zero-residue primaries.  A
T-eigenspace group is certified to split in an orientation when at least
two of its joint primaries are nonzero there.  The groups are unions of the
joint R,T coarse blocks; R labels of the primaries are not determined
here.  Exact dimensions and ranks remain open.
"""
from __future__ import annotations

import json
from math import comb
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))

from bounded_spin_crossing import paths
from bounded_spin_six_vertex import monomial_entry

PRIME = 5
RT = ROOT / 'evidence/bounded-spin-twelve-RT-coefficient.json'
CUT = ROOT / 'evidence/bounded-spin-twelve-cut-trace-form.json'
IDEMPOTENTS = ROOT / 'evidence/bounded-spin-twelve-coefficient-idempotents.json'
OUTPUT = ROOT / 'evidence/bounded-spin-twelve-separator-split.json'


def original_form(four_paths):
    length = len(four_paths[0])
    result = [0] * (length + 1)
    for (_, b, c), value in monomial_entry(four_paths).items():
        for j in range(2 * b + 1):
            result[2 * c + j] = (result[2 * c + j]
                                 + value * comb(2 * b, j) * (-1) ** j) % PRIME
    return result


def diagonal_witness(basis, idempotent_data, index, path_indices):
    idempotents = json.loads(IDEMPOTENTS.read_text())['primitive_central_idempotents']
    entry = idempotents[index]
    if entry['index'] != index:
        raise ValueError('unexpected idempotent index')
    vector = entry['idempotent']
    ps = paths(12)
    a, b = ps[path_indices[0]], ps[path_indices[1]]
    row_form = original_form((a, b, b, a))
    column_form = original_form((a, b, a, b))
    row_basis = [row_form[label] for label in basis]
    column_basis = [column_form[label] for label in basis]
    row_entry = sum(vector[k] * row_basis[k] for k in range(len(basis))) % PRIME
    column_entry = sum(vector[k] * column_basis[k] for k in range(len(basis))) % PRIME
    return {
        'idempotent_coefficient_vector': vector,
        'regular_rank': entry['regular_rank'],
        'path_indices': list(path_indices),
        'row_form_on_basis': row_basis,
        'column_form_on_basis': column_basis,
        'row_diagonal_entry': row_entry,
        'column_diagonal_entry': column_entry,
        'nonzero_in_both_orientations': row_entry % PRIME != 0 and column_entry % PRIME != 0,
    }


def idempotent_diagonal_witnesses(basis):
    if basis != [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
        raise ValueError('unexpected coefficient basis')
    witnesses = {
        '0': diagonal_witness(basis, IDEMPOTENTS, 0, (38, 82)),
        '4': diagonal_witness(basis, IDEMPOTENTS, 4, (0, 0)),
        '6': diagonal_witness(basis, IDEMPOTENTS, 6, (0, 0)),
    }
    for index, witness in witnesses.items():
        if not witness['nonzero_in_both_orientations']:
            raise ValueError(f'idempotent {index} diagonal witness vanished')
    return witnesses


def main() -> None:
    rt = json.loads(RT.read_text())
    cut = json.loads(CUT.read_text())
    if not rt.get('passed') or cut.get('status') != 'pass':
        raise ValueError('input certificates are not marked as passed')
    if rt['basis_indices'] != cut['basis_indices']:
        raise ValueError('input basis orientations differ')
    joint = {int(k): tuple(v) for k, v in
             rt['class_sum_separator_joint_characters'].items()}
    residues = cut['idempotent_trace_residues']
    if len(joint) != 8 or len(residues['row']) != 8 or len(residues['column']) != 8:
        raise ValueError('expected eight primaries in both inputs')
    if set(joint) != set(range(8)):
        raise ValueError('unexpected idempotent indices')
    witnesses = idempotent_diagonal_witnesses(rt['basis_indices'])

    groups = {}
    for index, (t_value, s_value) in joint.items():
        groups.setdefault(t_value, []).append(index)
    expected_groups = {0: [0, 6], 1: [1, 2, 5, 7], 4: [3, 4]}
    if groups != expected_groups:
        raise ValueError('unexpected T coarse groups')

    report = {}
    for t_value, indices in groups.items():
        distinct_s = sorted({joint[index][1] for index in indices})
        row_nonzero = [index for index in indices if residues['row'][index] % 5]
        column_nonzero = [index for index in indices
                          if residues['column'][index] % 5]
        for index in indices:
            witness = witnesses.get(str(index))
            if witness is not None and witness['row_diagonal_entry'] % 5:
                row_nonzero = sorted(set(row_nonzero + [index]))
            if witness is not None and witness['column_diagonal_entry'] % 5:
                column_nonzero = sorted(set(column_nonzero + [index]))
        report[str(t_value)] = {
            'idempotents': indices,
            'separator_characters': distinct_s,
            'row_nonzero_primary_residues': row_nonzero,
            'column_nonzero_primary_residues': column_nonzero,
            'row_split_certified': len(row_nonzero) >= 2,
            'column_split_certified': len(column_nonzero) >= 2,
            'row_split_obstructed_by_zero_residue':
                len(row_nonzero) < 2,
            'column_split_obstructed_by_zero_residue':
                len(column_nonzero) < 2,
        }

    if report['1']['row_split_certified'] is not True:
        raise ValueError('T=1 row split was expected')
    if report['1']['column_split_certified'] is not True:
        raise ValueError('T=1 column split was expected')
    for t_value in (0, 4):
        if report[str(t_value)]['row_split_certified'] is not True:
            raise ValueError(f'T={t_value} row split was expected')
        if report[str(t_value)]['column_split_certified'] is not True:
            raise ValueError(f'T={t_value} column split was expected')

    receipt = {
        'passed': True,
        'prime': 5,
        'length': 12,
        'basis_indices': rt['basis_indices'],
        'joint_characters_T_S': {str(k): list(v) for k, v in joint.items()},
        'idempotent_trace_residues': residues,
        'idempotent_diagonal_witnesses': witnesses,
        't_eigenspace_groups': report,
        'certified_split_t_groups': {
            'row': [int(t) for t, row in report.items()
                    if row['row_split_certified']],
            'column': [int(t) for t, row in report.items()
                       if row['column_split_certified']],
        },
        'unresolved_t_groups': {
            'row': [int(t) for t, row in report.items()
                    if not row['row_split_certified']],
            'column': [int(t) for t, row in report.items()
                       if not row['column_split_certified']],
        },
        'scope': (
            'Modulo-five separator splitting of the twelve-site T-eigenspace '
            'groups, using the corrected cut-trace functional residues, the '
            'coefficient-algebra joint T,S characters, and explicit nonzero '
            'diagonal matrix entries of the central idempotent actions. All '
            'three T groups split in both orientations: T=1 by the nonzero '
            'residues of primaries 1 and 7, T=0 by explicit entries of '
            'idempotents 0 and 6, and T=4 by the residue of primary 3 and an '
            'explicit entry of idempotent 4. Because each idempotent action is '
            'a projection, a zero trace residue forces its rank to be a '
            'multiple of five; combined with the witness this gives rank at '
            'least five for idempotents 0, 4, and 6. These T groups are unions '
            'of the joint R,T coarse blocks; R labels of the primaries are not '
            'determined here. This is not an exact dimension or rank '
            'computation and no score change is claimed.'
        ),
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({
        'passed': True,
        'certified_split_t_groups': receipt['certified_split_t_groups'],
        'unresolved_t_groups': receipt['unresolved_t_groups'],
        'idempotent_diagonal_witnesses': witnesses,
        't_eigenspace_groups': report,
    }, indent=2))


if __name__ == '__main__':
    main()
