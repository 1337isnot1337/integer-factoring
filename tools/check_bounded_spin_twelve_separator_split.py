#!/usr/bin/env python3
"""Modulo-five separator splitting of the twelve-site cut coarse blocks.

The eleven/eight idempotent trace residues from the cut trace form give
nonzero rank residues for the joint T,S primaries.  A T coarse block is
certified to split in an orientation when at least two of its joint
primaries have nonzero residue there.  Exact dimensions and ranks remain
open.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RT = ROOT / 'evidence/bounded-spin-twelve-RT-coefficient.json'
CUT = ROOT / 'evidence/bounded-spin-twelve-cut-trace-form.json'
OUTPUT = ROOT / 'evidence/bounded-spin-twelve-separator-split.json'


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
    if report['4']['row_split_certified'] is not True:
        raise ValueError('T=4 row split was expected')
    if report['4']['column_split_certified'] is not True:
        raise ValueError('T=4 column split was expected')
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
        'coarse_groups': report,
        'certified_split_groups': {
            'row': [int(t) for t, row in report.items()
                    if row['row_split_certified']],
            'column': [int(t) for t, row in report.items()
                       if row['column_split_certified']],
        },
        'unresolved_groups': {
            'row': [int(t) for t, row in report.items()
                    if not row['row_split_certified']],
            'column': [int(t) for t, row in report.items()
                       if not row['column_split_certified']],
        },
        'scope': (
            'Modulo-five separator splitting of the twelve-site T coarse '
            'blocks, using the cut-trace idempotent rank residues and the '
            'coefficient-algebra joint T,S characters. T=1 and T=4 split in '
            'both orientations; the T=0 group is not certified because the '
            'first idempotent has residue zero in both cuts. This is not an '
            'exact dimension or rank computation and no score change is claimed.'
        ),
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({
        'passed': True,
        'certified_split_groups': receipt['certified_split_groups'],
        'unresolved_groups': receipt['unresolved_groups'],
        'coarse_groups': report,
    }, indent=2))


if __name__ == '__main__':
    main()
