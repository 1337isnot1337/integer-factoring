#!/usr/bin/env python3
"""Exact T-primary dimensions of the twelve-site cut modules.

The coarse (T,R) joint generalized dimensions are exact module dimensions.
Summing the two R signs gives the exact T-primary dimensions.  These are
not the trace-functional values computed by the T-power checker; in
particular their residues need not agree with the cut trace vectors.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'evidence/bounded-spin-coarse-twelve.json'
EXPECTED = {(0, 1): 1320, (0, 4): 165, (1, 1): 2015,
            (1, 4): 1266, (4, 1): 670, (4, 4): 2485}


def check() -> dict:
    receipt = json.loads(EVIDENCE.read_text())
    joint = {(int(a), int(b)): int(value)
             for a, b, value in receipt['coarse_dimensions']}
    assert joint == EXPECTED
    exact = {t: sum(value for (tt, _), value in joint.items() if tt == t)
             for t in (0, 1, 4)}
    assert exact == {0: 1485, 1: 3281, 4: 3155}
    assert sum(exact.values()) == 89 * 89
    return {
        'passed': True,
        'length': 12,
        'prime': 5,
        'joint_T_R_dimensions': [[a, b, value]
                                 for (a, b), value in sorted(joint.items())],
        'coarse_T_primary_dimensions': {str(k): v for k, v in exact.items()},
        'T_residues_mod_5': {str(k): v % 5 for k, v in exact.items()},
        'sum': sum(exact.values()),
        'scope': ('Exact coarse class-sum T-primary module dimensions obtained '
                  'by summing the two exact R signs of the coarse (T,R) '
                  'generalized eigenspaces. These are dimensions, not cut ranks. '
                  'The coarse class-sum T labels are not yet reconciled with the '
                  'coefficient-table T labels (the residues suggest exchanging '
                  'labels 1 and 4), so do not combine these dimensions with the '
                  'joint T,S characters. Individual central-idempotent ranks and '
                  'the two-cut pairing remain open.'),
    }


if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
