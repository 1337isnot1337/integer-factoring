#!/usr/bin/env python3
"""Verify the recorded twelve-site coefficient multiplication table."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_bounded_spin_twelve_coefficient_table import (
    BASIS, PRIME, PROBE, SOURCE_NAMES, rank_mod,
)

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'evidence/bounded-spin-twelve-coefficient-table.json'


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    assert saved['passed'] is True
    assert saved['prime'] == PRIME and saved['length'] == 12
    assert saved['basis_indices'] == BASIS
    assert saved['probe_sha256'] == sha256(PROBE.read_bytes()).hexdigest()
    assert len(saved['witness_quadruples']) == 10
    entry_matrix = saved['entry_matrix']
    assert len(entry_matrix) == 10 and all(len(row) == 10 for row in entry_matrix)
    assert rank_mod(entry_matrix) == 10
    unit = saved['unit']
    assert len(unit) == 10
    table = saved['table']
    assert len(table) == 10 and all(len(row) == 10 for row in table)
    assert all(len(row) == 10 for block in table for row in block)
    assert saved['commutative'] is True
    assert saved['associative'] is True
    assert saved['unit_ok'] is True

    for i in range(10):
        for j in range(10):
            assert table[i][j] == table[j][i]
    for i in range(10):
        for j in range(10):
            for k in range(10):
                left = [0] * 10
                for a in range(10):
                    factor = table[i][j][a]
                    if factor:
                        for b in range(10):
                            left[b] = (left[b] + factor * table[a][k][b]) % PRIME
                right = [0] * 10
                for a in range(10):
                    factor = table[j][k][a]
                    if factor:
                        for b in range(10):
                            right[b] = (right[b] + factor * table[i][a][b]) % PRIME
                assert left == right

    regular = saved['regular']
    for i in range(10):
        for j in range(10):
            for k in range(10):
                assert regular[i][k][j] == table[i][j][k]
        product = [sum(regular[i][k][j] * unit[j] for j in range(10)) % PRIME
                   for k in range(10)]
        assert product == [int(k == i) for k in range(10)]

    for holdout in saved['holdouts']:
        assert len(holdout['product_values']) == 10
        assert all(len(row) == 10 for row in holdout['product_values'])
    for name in SOURCE_NAMES:
        assert sha256((ROOT / 'tools' / name).read_bytes()).hexdigest() == saved['source_sha256'][name]

    receipt = {
        'passed': True,
        'dimension': 10,
        'length': 12,
        'basis_indices': BASIS,
        'unit': unit,
        'commutative_associative_unital': True,
        'holdout_rows_checked': len(saved['holdouts']),
        'entry_witness_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'probe_sha256': saved['probe_sha256'],
        'scope': ('Replay of the exact ten-dimensional multiplication table, its '
                  'unit, and its regular representation. The build, not this '
                  'checker, recomputes the two holdout entries directly. No '
                  'primary decomposition, cut-rank bound, or score change is '
                  'claimed.'),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'dimension': 10, 'unit': unit}))


if __name__ == '__main__':
    main()
