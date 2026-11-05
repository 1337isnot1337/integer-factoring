#!/usr/bin/env python3
"""Track residual pair-row lattice defects across extra exponent windows.

Loads the committed window-4 rank sweep, takes the targets whose HNF
row-lattice index in the full relation kernel is greater than one, and
recomputes them at windows 5 and 6.  The point is to distinguish a
one-window saturation defect from a stable arithmetic defect: if extra pair
exponents do not reach index one, larger windows alone are not the missing
mechanism.  Finite diagnostic only; no uniform rank or saturation theorem.
"""
import argparse
import json

from check_mixed_pair_rank_sweep import probe, summary

W4_PATH = 'evidence/mixed-pair-rank-sweep-w4.json'


def run():
    with open(W4_PATH) as handle:
        w4 = json.load(handle)
    w4_by_target = {record['target']: record for record in w4['results']}
    targets = sorted(
        target for target, record in w4_by_target.items()
        if record.get('lattice_multiple', 1) > 1)
    assert targets
    w5 = probe(max(targets), 5, explicit_targets=targets)
    w5_by_target = {record['target']: record for record in w5['results']}
    residual = sorted(
        target for target, record in w5_by_target.items()
        if record.get('lattice_multiple', 1) > 1)
    if residual:
        w6 = probe(max(residual), 6, explicit_targets=residual)
        w6_by_target = {record['target']: record for record in w6['results']}
    else:
        w6 = None
        w6_by_target = {}
    w6_residual = sorted(
        target for target, record in w6_by_target.items()
        if record.get('lattice_multiple', 1) > 1)
    records = []
    for target in targets:
        record = {
            'target': target,
            'order': w4_by_target[target]['order'],
            'target_index': w4_by_target[target]['index'],
            'w4_lattice_multiple': w4_by_target[target]['lattice_multiple'],
            'w5_lattice_multiple': w5_by_target[target].get(
                'lattice_multiple'),
            'w5_pair_hits': w5_by_target[target]['target_pair_hits'],
        }
        if target in w6_by_target:
            record['w6_lattice_multiple'] = w6_by_target[target].get(
                'lattice_multiple')
            record['w6_pair_hits'] = w6_by_target[target]['target_pair_hits']
        records.append(record)
    return {
        'scope': ('Finite saturation-depth diagnostic on the 40-bit '
                  'order-198850 cycle. Window-4 targets with lattice index '
                  'above one are recomputed at windows 5 and 6. No uniform '
                  'saturation, rank, or complete-bound claim.'),
        'w4_targets': targets,
        'w4_count': len(targets),
        'w5_summary': summary(w5),
        'w5_residual_targets': residual,
        'w5_residual_count': len(residual),
        'w6_summary': None if w6 is None else summary(w6),
        'w6_residual_targets': w6_residual,
        'w6_residual_count': len(w6_residual),
        'records': records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    arguments = parser.parse_args()
    receipt = run()
    text = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    if arguments.output:
        with open(arguments.output, 'w') as handle:
            handle.write(text)
    print(json.dumps({
        'w4_count': receipt['w4_count'],
        'w5_lattice_multiples': receipt['w5_summary']['lattice_multiples'],
        'w5_residual_count': receipt['w5_residual_count'],
        'w5_residual_targets': receipt['w5_residual_targets'],
        'w6_lattice_multiples': (None if receipt['w6_summary'] is None
                                 else receipt['w6_summary']['lattice_multiples']),
        'w6_residual_count': receipt['w6_residual_count'],
        'w6_residual_targets': receipt['w6_residual_targets'],
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
