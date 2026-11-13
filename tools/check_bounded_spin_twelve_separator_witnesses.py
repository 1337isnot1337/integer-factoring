#!/usr/bin/env python3
"""Direct nonzero entries for the remaining twelve-site central idempotents.

Idempotent 4 is the second T=4 primary and idempotent 6 is the second T=0
primary.  A nonzero row-cut matrix entry certifies that each cut-module
image is nonzero even though the corresponding trace residue vanishes.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))

from bounded_spin_crossing import paths
from check_bounded_spin_twelve_RT_coefficient import PRIME, original_form

IDEMPOTENTS = ROOT / 'evidence/bounded-spin-twelve-coefficient-idempotents.json'


def entry(coefficients, basis, ps, indices):
    form = original_form(tuple(ps[i] for i in indices))
    return sum(coefficients[k] * form[basis[k]] for k in range(len(basis))) % PRIME


def check() -> dict:
    receipt = json.loads(IDEMPOTENTS.read_text())
    basis = receipt['basis_indices']
    central = {item['index']: item['idempotent']
               for item in receipt['primitive_central_idempotents']}
    ps = paths(12)
    assert len(ps) == 89
    witnesses = {
        4: entry(central[4], basis, ps, (0, 0, 0, 0)),
        6: entry(central[6], basis, ps, (0, 0, 0, 0)),
    }
    assert witnesses == {4: 2, 6: 4}
    return {
        'passed': True,
        'prime': PRIME,
        'length': 12,
        'path_index': 0,
        'path': list(ps[0]),
        'idempotent_entries': {str(k): v for k, v in witnesses.items()},
        'conclusion': ('Idempotents 4 (T=4,S=0) and 6 (T=0,S=4) have nonzero '
                       'row-cut images. Together with the published nonzero '
                       'residues or the idempotent-0 witness, every T group '
                       'has two nonzero certified central components.'),
        'scope': ('Direct diagonal row-cut entry witnesses only. Exact ranks, '
                  'whether the two components exhaust each T-primary block, '
                  'the two-cut pairing, and any score change remain open.'),
    }


if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
