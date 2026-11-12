#!/usr/bin/env python3
"""Modulo-five separator splitting of the twelve-site cut coarse blocks.

The eleven/eight idempotent trace residues from the cut trace form give
nonzero rank residues for the joint T,S primaries.  A T-eigenspace group is
certified to split in an orientation when at least two of its joint
primaries have nonzero residue there.  The groups are unions of the
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


def idempotent_zero_witnesses(basis):
    if basis[8] != 9:
        raise ValueError('unexpected coefficient-nine basis position')
    idempotents = json.loads(IDEMPOTENTS.read_text())
    zero = idempotents['primitive_central_idempotents'][0]
    if zero['index'] != 0:
        raise ValueError('unexpected first idempotent')
    vector = zero['idempotent']
    if any(value for index, value in enumerate(vector) if index != 8):
        raise ValueError('idempotent zero is not supported on coefficient nine')
    if vector[8] != 4:
        raise ValueError('unexpected coefficient of idempotent zero')
    ps = paths(12)
    row = original_form((ps[38], ps[82], ps[82], ps[38]))[9]
    column = original_form((ps[38], ps[82], ps[38], ps[82]))[9]
    if row != 1 or column != 4:
        raise ValueError('unexpected diagonal witnesses for idempotent zero')
    if vector[8] * row % PRIME == 0 or vector[8] * column % PRIME == 0:
        raise ValueError('idempotent zero diagonal witness vanished')
    return {
        'idempotent_coefficient_vector': vector,
        'basis_position': 8,
        'coefficient': basis[8],
        'row_path_indices': [38, 82],
        'row_degree_nine_coefficient': row,
        'row_diagonal_entry': vector[8] * row % PRIME,
        'column_path_indices': [38, 82],
        'column_degree_nine_coefficient': column,
        'column_diagonal_entry': vector[8] * column % PRIME,
        'nonzero_in_both_orientations': True,
    }


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
    zero_witness = idempotent_zero_witnesses(rt['basis_indices'])

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
        if t_value == 0:
            row_nonzero = sorted(row_nonzero + [0])
            column_nonzero = sorted(column_nonzero + [0])
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
    if report['4']['row_split_certified'] is not False:
        raise ValueError('T=4 row split should not be certified')
    if report['4']['column_split_certified'] is not False:
        raise ValueError('T=4 column split should not be certified')
    if report['0']['row_split_certified'] is not False:
        raise ValueError('T=0 split should remain uncertified')
    if report['0']['column_split_certified'] is not False:
        raise ValueError('T=0 split should remain uncertified')

    receipt = {
        'passed': True,
        'prime': 5,
        'length': 12,
        'basis_indices': rt['basis_indices'],
        'joint_characters_T_S': {str(k): list(v) for k, v in joint.items()},
        'idempotent_trace_residues': residues,
        'idempotent_zero_diagonal_witnesses': zero_witness,
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
            'groups, using the corrected cut-trace idempotent rank residues '
            'and the coefficient-algebra joint T,S characters. T=1 splits '
            'into two nonzero primaries (S=0 and S=1) in both orientations. '
            'T=4 has only one nonzero residue and is not certified to split. '
            'For T=0, the first idempotent has zero trace residue but an '
            'explicit nonzero diagonal entry in each orientation, so its rank '
            'is a positive multiple of five; the second T=0 idempotent still '
            'has zero residue, so only one T=0 component is certified and the '
            'group is not certified to split. These T groups are unions of the '
            'joint R,T coarse blocks; R labels of the primaries are not '
            'determined here. This is not an exact dimension or rank '
            'computation and no score change is claimed.'
        ),
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({
        'passed': True,
        'certified_split_t_groups': receipt['certified_split_t_groups'],
        'unresolved_t_groups': receipt['unresolved_t_groups'],
        'idempotent_zero_diagonal_witnesses': zero_witness,
        't_eigenspace_groups': report,
    }, indent=2))


if __name__ == '__main__':
    main()
