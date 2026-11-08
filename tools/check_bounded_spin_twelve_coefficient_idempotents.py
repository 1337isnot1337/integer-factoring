#!/usr/bin/env python3
"""Recompute and verify the twelve-site coefficient idempotent certificate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import build_bounded_spin_twelve_coefficient_idempotents as build

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'evidence/bounded-spin-twelve-coefficient-idempotents.json'


def check():
    stored = json.loads(EVIDENCE.read_text())
    recomputed = build.compute()
    stored_math = {k: v for k, v in stored.items() if k != 'seconds'}
    recomputed_math = {k: v for k, v in recomputed.items() if k != 'seconds'}
    assert stored_math == recomputed_math, 'stored certificate does not reproduce'

    blocks = stored['primitive_central_idempotents']
    assert len(blocks) == 8
    assert sorted(b['regular_rank'] for b in blocks) == [1] * 6 + [2] * 2
    assert sum(1 for b in blocks if b['kind'] == 'scalar') == 6
    assert sum(1 for b in blocks if b['kind'] == 'dual') == 2
    assert stored['trace_form']['rank'] == 8
    assert stored['trace_form']['radical_dimension'] == 2
    assert stored['trace_form']['square_zero'] is True
    assert stored['structure'] == 'F5^6 x F5[e]/(e^2) x F5[e]/(e^2)'
    assert stored['checks']['regular_rank_sum'] == 10
    return {
        'status': 'pass',
        'structure': stored['structure'],
        'idempotents': len(blocks),
        'regular_rank_partition': sorted(b['regular_rank'] for b in blocks),
        'trace_form': stored['trace_form']['rank'],
        'radical_dimension': stored['trace_form']['radical_dimension'],
        'table_sha256': stored['table_sha256'],
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
