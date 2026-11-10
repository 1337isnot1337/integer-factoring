#!/usr/bin/env python3
"""Recompute and verify the twelve-site R,T orbit-algebra certificate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import build_bounded_spin_twelve_RT_orbit_algebra as build

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'evidence/bounded-spin-twelve-RT-orbit-algebra.json'


def check():
    stored = json.loads(EVIDENCE.read_text())
    recomputed = build.compute()
    stored_math = {key: value for key, value in stored.items()
                   if key != 'seconds'}
    recomputed_math = {key: value for key, value in recomputed.items()
                       if key != 'seconds'}
    assert stored_math == recomputed_math, 'stored certificate does not reproduce'

    blocks = stored['primitive_central_idempotents']
    assert len(blocks) == 6
    assert sorted(block['regular_rank'] for block in blocks) == [1] * 4 + [2] * 2
    assert sum(1 for block in blocks if block['kind'] == 'scalar') == 4
    assert sum(1 for block in blocks if block['kind'] == 'dual') == 2
    assert stored['trace_form']['rank'] == 6
    assert stored['trace_form']['radical_dimension'] == 2
    assert stored['trace_form']['square_zero'] is True
    assert stored['orbit_dimension'] == 8
    assert stored['full_coefficient_dimension'] == 10
    assert stored['codimension'] == 2
    assert stored['structure'] == 'F5^4 x F5[e]/(e^2) x F5[e]/(e^2)'
    assert stored['full_coefficient_structure'] == (
        'F5^6 x F5[e]/(e^2) x F5[e]/(e^2)')
    assert stored['checks']['regular_rank_sum'] == 8
    return {
        'status': 'pass',
        'structure': stored['structure'],
        'full_coefficient_structure': stored['full_coefficient_structure'],
        'orbit_dimension': stored['orbit_dimension'],
        'full_coefficient_dimension': stored['full_coefficient_dimension'],
        'idempotents': len(blocks),
        'regular_rank_partition': sorted(
            block['regular_rank'] for block in blocks),
        'trace_form': stored['trace_form']['rank'],
        'radical_dimension': stored['trace_form']['radical_dimension'],
        'orbit_evidence_sha256': stored['orbit_evidence_sha256'],
        'scope': stored['scope'],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text + '\n')
    print(text)


if __name__ == '__main__':
    main()
